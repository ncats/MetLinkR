#' Run total harmonization of metabolites
#'
#' @param inputcsv character of pathway to csv containing names of files with metabolite names
#' @param mapping_library_format string, print library in "wide" format, "long" format, or "both"
#'
#' @return a dataframe containing harmonized files across all input files listed in inputcsv
#' @examples
#' \dontrun{
#' finalHarmonized <- runHarmonization(
#'   inputcsv = "HarmInputFiles.csv"
#' )
#' }
#'
#' @importFrom rlang .data
#'
#' @export

harmonizeInputSheets <- function(inputcsv,
                                 mapping_library_format = "both", n_cores = 1) {
  start_time <- Sys.time()
  cluster <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  if(!(mapping_library_format %in% c("long","wide","both"))){
    stop("'mapping library format' must be one of 'long', 'wide', or 'both'")
  }
  ##########################################################################
  ## 1. Read in input files. Output a list of dataframes                  ##
  ##########################################################################
  myinputfiles <- utils::read.csv(inputcsv, header = T)
  list_input_files <- readInputCSVs(inputcsv)
  myinputfiles_list <- as.list(data.frame(t(myinputfiles)))
  myinputfiles_list <- lapply(myinputfiles_list, function(x) {
    out <- t(data.frame(x))
    colnames(out) <- colnames(myinputfiles)
    return(as.data.frame(out))
  })
  message("(1/5) Imported files")

  ##########################################################################
  ## 2. Initial RefMet mappings                                           ##
  ##########################################################################
  mapped_list_input_files <- foreach(i = 1:length(list_input_files)) %dopar% {
      metLinkR:::queryRefMet(
        input_df = list_input_files[[i]],
        filename = myinputfiles_list[[i]]$ShortFileName,
        HMDB_col = myinputfiles_list[[i]]$HMDB,
        metab_col = myinputfiles_list[[i]]$Metabolite_Name,
        CID_col = myinputfiles_list[[i]]$PubChem_CID,
        KEGG_col =myinputfiles_list[[i]]$KEGG,
        LM_col = myinputfiles_list[[i]]$LIPIDMAPS,
        CHEBI_col = myinputfiles_list[[i]]$chebi
      )
  }

  refmet_mapped_ids <- lapply(
    mapped_list_input_files,
    function(x) {
      return(x %>%
               mutate("Origin" = "Original input") %>%
               ## mutate("classFlag" = "RefMet") %>%
               dplyr::filter(`Standardized name` != "-") %>%
               dplyr::group_by(rownum) %>%
               dplyr::filter(priority == min(priority)) %>%
               ## dplyr::ungroup %>%
               as.data.frame() %>%
               dplyr::select(-c(.data$priority)))
    }
  )

  message("(2/5) Performed initial RefMet mapping")
  ################################################################################
  ## 3. Assemble synonyms table from RaMP-DB for missed IDs, perform mass check ##
  ################################################################################
  mapped_rownums <- mapply(function(x,y) {
    return(x %>%
             dplyr::filter(.data$rownum %in% y$rownum) %>%
             dplyr::pull(.data$rownum) %>%
             unique
           )
  }, x = mapped_list_input_files, y = refmet_mapped_ids)

  refmet_unmapped_ids <- mapply(function(x,y){
    return(x[-y,])
  }, x = list_input_files, y = mapped_rownums, SIMPLIFY = FALSE)

  missed_ids <- mapply(function(x, y) {
    if(nrow(x)==0){
      return(NA)
    }else{
      return(extract_identifiers(
        input_df = replaceEmptys(x),
        HMDB_col = y$HMDB,
        metab_col = y$Metabolite_Name,
        CID_col = y$PubChem_CID,
        ramp_prefixes = TRUE
      ))
    }
  }, x = refmet_unmapped_ids, y = myinputfiles_list, SIMPLIFY = FALSE)

  parallel::clusterEvalQ(cluster,db <<- RaMP::RaMP())
  synonym_table_list <- parallel::parLapply(cl=cluster,missed_ids, queryRampSynonyms)
  message("(3/5) Found RaMP synonyms for unmapped inputs")
  ##########################################################################
  ## 4. Re-Query RefMet with synonym table                                ##
  ##########################################################################

    mapped_list_synonyms <- foreach(i = 1:length(list_input_files)) %dopar% {
      metLinkR:::queryRefMet(
        input_df = synonym_table_list[[i]],
        filename = paste0("synonym_table_",myinputfiles_list[[i]]$ShortFileName),
        HMDB_col = NA,
        metab_col = "Synonym",
        CID_col = NA,
        KEGG_col = NA,
        LM_col = NA,
        CHEBI_col = NA
      )
  }
  
  refmet_mapped_synonyms <- list()
  for (i in 1:length(mapped_list_synonyms)) {
    if(length(mapped_list_synonyms[[i]])==1){
      refmet_mapped_synonyms[[i]] <- NA
    }else{
      refmet_mapped_synonyms[[i]] <- mapped_list_synonyms[[i]] %>%
        dplyr::filter(.data$`Standardized name` != "-") %>%
        dplyr::left_join(synonym_table_list[[i]], c("Input name" = "Synonym")) %>%
        dplyr::select(-.data$`Input name`) %>%
        dplyr::rename("Input name" = "Input") %>%
        dplyr::distinct() %>%
        dplyr::mutate("Origin" = "Synonym")
    }
  }

  ## Fix rownums
  refmet_mapped_synonyms <- mapply(function(x,y){
    if(is.logical(y)){
      return(NA)
    }else{
      y = y %>% dplyr::filter(!is.na(.data$`Input name`))
      rownums <- c()
      for(i in 1:nrow(y)){
        matches <- which(
          apply(x, 1,
                function(a) {
                  strip_prefixes(y[i, "Input name"]) %in% a}))
        if(length(matches)!=0){
          rownums <- c(rownums,
                       matches)
        }else{
          rownums <- c(rownums, NA)
        }
      }
      if(nrow(y)==0){
        return(NA)
      }else{
        return(y %>%
                 dplyr::mutate(rownum = rownums))
      }
    }
  },x = list_input_files, y = refmet_mapped_synonyms)
  message("(4/5) Queried RaMP synonyms in RefMet")

  mapped_input_list <- mapply(function(x, y) {
    if(length(y)==1){
      x
    }else{
      rbind(
        x,
        y[, colnames(x)]
      )}
  }, x = refmet_mapped_ids, y = refmet_mapped_synonyms, SIMPLIFY = FALSE)

  ##########################################################################
  ## 4.5 Display mapping rate summary                                     ##
  ##########################################################################
  mapping_rates <-
    calculate_mapping_rates(mapped_input_list, list_input_files, myinputfiles)

  ##########################################################################
  ## 5. Append output to inputs. Merge mapped files and write to csv      ##
  ##########################################################################
  dir.create("metLinkR_output", showWarnings = FALSE)
  appended_inputs <- append_standard_names(mapped_input_list, list_input_files)

  silent <- mapply(function(x,y){
    utils::write.csv(x,
              file = paste0("metLinkR_output/",gsub(".csv","",y),"_metLinkR.csv"),
              row.names = FALSE)
  },x = appended_inputs, y = myinputfiles$FileNames)


  mapping_library <- merge_files(mapped_input_list,myinputfiles)
  multimappings <- find_multimapped_metabolites(mapping_library,myinputfiles)

  if(mapping_library_format=="wide"){
    xlsx::write.xlsx(as.data.frame(mapping_library),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library")
  }else if(mapping_library_format=="long"){
    mapping_library <- pivot_mapping_library(mapping_library)
    xlsx::write.xlsx(as.data.frame(mapping_library),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library")
  }else{
    xlsx::write.xlsx(as.data.frame(mapping_library),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library Wide")
    mapping_library <- pivot_mapping_library(mapping_library)
    xlsx::write.xlsx(as.data.frame(mapping_library),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library Long", append = TRUE)
  }

  ## Write unmapped values 
  missed_mappings <- extract_missing_values(appended_inputs,myinputfiles)
  missed_mappings <- lapply(missed_mappings, function(x){
    return(x %>%
             dplyr::mutate(across(everything(),
                                  as.character)))
  })
  missed_mappings <- dplyr::bind_rows(as.vector(missed_mappings), .id = "Data File")
  xlsx::write.xlsx(missed_mappings,
                   file = paste0("metLinkR_output/mapping_library.xlsx"),
                   row.names=FALSE,sheetName="Missed Metabolites",
                   append = TRUE,showNA=FALSE)

  ## Write multi-mapped metabolites
  xlsx::write.xlsx(multimappings,
                   file = paste0("metLinkR_output/mapping_library.xlsx"),
                   row.names=FALSE,sheetName="MultiMapped Metabolites",
                   append = TRUE,showNA=FALSE)
  ## Write text log
  write_txt_log(start_time,myinputfiles)
  names(mapping_rates[[2]]) = names(mapped_list_input_files) =
    names(mapped_list_synonyms) = myinputfiles$ShortFileName
  
  ## Write PDF report
  write_html_report(mapping_rates,
                   mapped_list_input_files,
                   mapped_list_synonyms)


  print("(5/5) Wrote output files to metLinkR_output/")
  return(mapping_library)
}

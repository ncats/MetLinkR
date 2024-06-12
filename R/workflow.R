#' Run total harmonization of metabolites
#'
#' @param inputcsv character of pathway to csv containing names of files with metabolite names
#' @param outputFileName character of output file
#' @param writecsv boolean for whether or not to write a csv file of the dataframe with harmonized files
#' @param long_mapping_library boolean for whether or not to pivot the mapping library
#'
#' @return a dataframe containing harmonized files across all input files listed in inputcsv
#' @examples
#' \dontrun{
#' finalHarmonized <- runHarmonization(
#'   inputcsv = "HarmInputFiles.csv",
#'   outputFileName = "exampleharmonizedfile"
#' )
#' }
#'
#' @importFrom rlang .data
#'
#' @export

harmonizeInputSheets <- function(inputcsv, outputFileName, writecsv,
                                 long_mapping_library = TRUE) {
  start_time <- Sys.time()
  ##########################################################################
  ## 1. Read in input files. Output a list of dataframes                  ##
  ##########################################################################
  myinputfiles <- utils::read.csv(inputcsv, header = T)
  list_input_files <- readInputCSVs(inputcsv)
  myinputfiles_list <- as.list(data.frame(t(myinputfiles)))
  myinputfiles_list <- parallel::mclapply(myinputfiles_list, function(x) {
    out <- t(data.frame(x))
    colnames(out) <- colnames(myinputfiles)
    return(as.data.frame(out))
  })
  message("(1/5) Imported files")

  ##########################################################################
  ## 2. Initial RefMet mappings                                           ##
  ##########################################################################
  mapped_list_input_files <- parallel::mcmapply(function(x, y) {
    queryRefMet(
      input_df = x,
      filename = y$ShortFileName,
      HMDB_col = y$HMDB,
      metab_col = y$Metabolite_Name,
      CID_col = y$PubChem_CID,
      KEGG_col =y$KEGG,
      LM_col = y$LIPIDMAPS,
      CHEBI_col = y$chebi
    )
  }, x = list_input_files,
  y = myinputfiles_list,
  SIMPLIFY = FALSE)

  refmet_mapped_ids <- parallel::mclapply(
    mapped_list_input_files,
    function(x) {
      return(x %>%
               dplyr::mutate("Origin" = "Raw input") %>%
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

  ## refmet_unmapped_ids <- list()
  ## for (i in 1:length(list_input_files)) {
  ##   refmet_unmapped_ids[[i]] <-
  ##     list_input_files[[i]][which(mapped_list_input_files[[i]]$`Standardized name` == "-"), ]
  ## }
  refmet_unmapped_ids <- mapply(function(x,y){
    return(x[-y,])
  }, x = list_input_files, y = mapped_rownums, SIMPLIFY = FALSE)

  missed_ids <- parallel::mcmapply(function(x, y) {
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

  db <<- RaMP::RaMP()
  synonym_table_list <- parallel::mclapply(missed_ids, queryRampSynonyms)
  message("(3/5) Found RaMP synonyms for unmapped inputs")
  ##########################################################################
  ## 4. Re-Query RefMet with synonym table                                ##
  ##########################################################################

  mapped_list_synonyms <- parallel::mcmapply(function(x, y) {
    queryRefMet(
      input_df = x,
      filename = paste0("synonym_table_", y$ShortFileName),
      HMDB_col = NA,
      metab_col = "Synonym",
      CID_col = NA,
      synonym_search = TRUE
    )
  }, x = synonym_table_list, y = myinputfiles_list, SIMPLIFY = FALSE)

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
  refmet_mapped_synonyms <- parallel::mcmapply(function(x,y){
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
      return(y %>%
               dplyr::mutate(rownum = rownums))
    }
    },x = list_input_files, y = refmet_mapped_synonyms)
  message("(4/5) Queried RaMP synonyms in RefMet")

  mapped_input_list <- parallel::mcmapply(function(x, y) {
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

  if(long_mapping_library){
    mapping_library <- pivot_mapping_library(mapping_library)
  }

  xlsx::write.xlsx(as.data.frame(mapping_library),
                   file = paste0("metLinkR_output/",outputFileName, ".xlsx"),
                   row.names=FALSE,sheetName="Mapping Library")

  ## Write unmapped values to second sheet
  missed_mappings <- extract_missing_values(appended_inputs,myinputfiles)
  missed_mappings <- dplyr::bind_rows(as.vector(missed_mappings), .id = "Data File")
  xlsx::write.xlsx(missed_mappings,
                   file = paste0("metLinkR_output/",outputFileName, ".xlsx"),
                   row.names=FALSE,sheetName="Missed Metabolites",
                   append = TRUE,showNA=FALSE)

  ## Write multi-mapped metabolites to third sheet
  xlsx::write.xlsx(multimappings,
                   file = paste0("metLinkR_output/",outputFileName, ".xlsx"),
                   row.names=FALSE,sheetName="MultiMapped Metabolites",
                   append = TRUE,showNA=FALSE)
  ## Write text log
  write_txt_log(start_time,myinputfiles)

  ## Write PDF report
  write_pdf_report(mapping_rates,
                   mapped_list_input_files,
                   mapped_list_synonyms)


  print("(5/5) Wrote output files to metLinkR_output/")
  return(mapping_library)
}

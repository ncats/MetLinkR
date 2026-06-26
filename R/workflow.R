#' Run total harmonization of metabolites
#'
#' @param inputcsv character of pathway to csv containing names of files with metabolite names
#' @param mapping_library_format string, print library in "wide" format, "long" format, or "both"
#' @param n_cores The number of cores available for parallel computation
#' @param use_ramp_synonyms Query metabolites missed by RefMet in the RaMP database
#' @param remove_parentheses_for_synonym_search Remove any parenthetical text before searching for synonyms in RaMP. May result in loss of structural specificity
#' @param use_metabolon_parsers Strip special characters typical of metabolon inputs
#' @param majority_vote When multiple IDs for a metabolite are provided and they produce different standardized names, choose the highest frequency output
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
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @export

harmonizeInputSheets <- function(inputcsv,
                                 mapping_library_format = "both", n_cores = 1,
                                 use_ramp_synonyms = TRUE,
                                 remove_parentheses_for_synonym_search = TRUE,
                                 use_metabolon_parsers = TRUE,
                                 majority_vote = TRUE,
                                 direct_input = FALSE,
                                 file_list = NA) {
  start_time <- Sys.time()
  cluster <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cluster)
  on.exit(parallel::stopCluster(cluster))
  if(!(mapping_library_format %in% c("long","wide","both"))){
    stop("'mapping library format' must be one of 'long', 'wide', or 'both'")
  }
  ##########################################################################
  ## 1. Read in input files. Output a list of dataframes                  ##
  ##########################################################################

  if(!direct_input){
    myinputfiles <- utils::read.csv(inputcsv, header = T)
    myinputfiles_list <- as.list(data.frame(t(myinputfiles)))
    myinputfiles_list <- lapply(myinputfiles_list, function(x) {
      out <- t(data.frame(x))
      colnames(out) <- colnames(myinputfiles)
      return(as.data.frame(out))
    })
    list_input_files <- readInputCSVs(inputcsv)
  }else{
    myinputfiles <- inputcsv
    myinputfiles <- read.table(text = unlist(myinputfiles), sep =",", header = TRUE, stringsAsFactors = FALSE)
    myinputfiles_list <- as.list(data.frame(t(myinputfiles)))
    myinputfiles_list <- lapply(myinputfiles_list, function(x) {
      out <- t(data.frame(x))
      colnames(out) <- colnames(myinputfiles)
      return(as.data.frame(out))
    })
    list_input_files <- lapply(file_list, function(x){
      return(read.table(text = x, sep =",", header = TRUE, stringsAsFactors = FALSE))
    })
  }

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
  
  initial_disagreements<-lapply(mapped_list_input_files,record_disagreements)

  refmet_mapped_ids <- lapply(mapped_list_input_files, filter_hits,
                              majority_vote = majority_vote)

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

  if(use_ramp_synonyms){
    db <- RaMP::RaMP()
    synonym_table_list <- lapply(
      missed_ids,
      queryRampSynonyms,
      db = db,
      use_metabolon_parsers = use_metabolon_parsers
    )

    if(remove_parentheses_for_synonym_search){
      synonym_table_list <- deparen_names(refmet_unmapped_ids, myinputfiles_list,
                                          synonym_table_list)
    }
  }else{
    synonym_table_list <- replicate(length(list_input_files), NA, simplify = FALSE)
  }
  
  if(use_ramp_synonyms){
    message("(3/5) Found RaMP synonyms for unmapped inputs")
  }else{
    message("(3/5) Skipped RaMP synonym search")
  }

  ##########################################################################
  ## 4. Re-Query RefMet with synonym table                                ##
  ##########################################################################
    mapped_list_synonyms <- lapply(seq_along(list_input_files), function(i) {
      synonym_input <- synonym_table_list[[i]]
      if (!is.data.frame(synonym_input) || nrow(synonym_input) == 0 ||
          !("Synonym" %in% colnames(synonym_input))) {
        return(NA)
      }

      scalarize_col <- function(x) {
        vapply(x, function(value) {
          if (length(value) == 0 || all(is.na(value))) {
            return(NA_character_)
          }
          as.character(value[[1]])
        }, character(1))
      }

      synonym_input <- data.frame(
        Synonym = scalarize_col(as.list(synonym_input$Synonym)),
        Input = if ("Input" %in% colnames(synonym_input)) {
          scalarize_col(as.list(synonym_input$Input))
        } else {
          rep(NA_character_, nrow(synonym_input))
        },
        stringsAsFactors = FALSE
      )

      synonym_input <- synonym_input %>%
        dplyr::filter(!is.na(.data$Synonym) & .data$Synonym != "") %>%
        dplyr::distinct(.data$Synonym, .data$Input, .keep_all = TRUE)
      if (nrow(synonym_input) == 0) {
        return(NA)
      }
      tryCatch(
        metLinkR:::queryRefMet(
        input_df = synonym_input,
        filename = paste0("synonym_table_",myinputfiles_list[[i]]$ShortFileName),
        HMDB_col = NA,
        metab_col = "Synonym",
        CID_col = NA,
        KEGG_col = NA,
        LM_col = NA,
        CHEBI_col = NA
      ),
      error = function(e) {
        utils::write.csv(
          synonym_input,
          file = paste0(
            "metLinkR_output/debug_synonym_input_",
            myinputfiles_list[[i]]$ShortFileName,
            ".csv"
          ),
          row.names = FALSE
        )
        stop(
          paste0(
            "Synonym re-query failed for ",
            myinputfiles_list[[i]]$ShortFileName,
            ". Debug input written to metLinkR_output/debug_synonym_input_",
            myinputfiles_list[[i]]$ShortFileName,
            ".csv. Original error: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      })
  })
  
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
                       matches[1])
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
  },x = list_input_files, y = refmet_mapped_synonyms, SIMPLIFY = FALSE)
  refmet_mapped_synonyms <- lapply(refmet_mapped_synonyms, function(x) {
    if (!is.data.frame(x) || nrow(x) == 0) {
      return(NA)
    }
    x <- x %>%
      dplyr::filter(!is.na(.data$rownum)) %>%
      dplyr::distinct()
    if (nrow(x) == 0) {
      return(NA)
    }
    filter_hits(x, majority_vote = majority_vote)
  })
  message("(4/5) Queried RaMP synonyms in RefMet")

  mapped_input_list <- mapply(function(x, y) {
    if(!is.data.frame(y) | length(y)==1){
      x
    }else{
      rbind(
        x,
        y[, colnames(x)]
      )}
  }, x = refmet_mapped_ids, y = refmet_mapped_synonyms, SIMPLIFY = FALSE)

  ##########################################################################
  ## 5.1 Grab metadata from RaMP                                          ##
  ##########################################################################
  ## if(check_stereochemistry){
  ##   mapped_input_list <- parallel::parLapply(cl=cluster,mapped_input_list,
  ##                                            makeStereochemFlags)
  ## }
  ## if(get_pathway_mappings){
  ##   mapped_input_list <- parallel::parLapply(cl=cluster,mapped_input_list,
  ##                                            getPathwayAndDBData)
  ## }
  ##########################################################################
  ## 5.2 Display mapping rate summary                                     ##
  ##########################################################################
  mapping_rates <-
    calculate_mapping_rates(mapped_input_list, list_input_files, myinputfiles)

  ##########################################################################
  ## 5.3 Append output to inputs. Merge mapped files and write to csv     ##
  ##########################################################################
  dir.create("metLinkR_output", showWarnings = FALSE)
  appended_inputs <- append_standard_names(mapped_input_list, list_input_files)

  silent <- mapply(function(x,y){
    utils::write.csv(x,
              file = paste0("metLinkR_output/",gsub(".csv","",y),"_metLinkR.csv"),
              row.names = FALSE)
  },x = appended_inputs, y = myinputfiles$FileNames)
  
  metadata <- assemble_metadata(mapped_input_list)
  mapping_library <- merge_files(mapped_input_list,myinputfiles)
  multimappings <- find_multimapped_metabolites(mapping_library,myinputfiles)
  mapping_library_long <- pivot_mapping_library(mapping_library)
  
  if(mapping_library_format=="wide"){
    xlsx::write.xlsx(as.data.frame(mapping_library),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library")
  }else if(mapping_library_format=="long"){
    xlsx::write.xlsx(as.data.frame(mapping_library_long),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library")
  }else{
    xlsx::write.xlsx(as.data.frame(mapping_library),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library Wide")
    xlsx::write.xlsx(as.data.frame(mapping_library_long),
                     file = paste0("metLinkR_output/mapping_library.xlsx"),
                     row.names=FALSE,sheetName="Mapping Library Long", append = TRUE)
  }

  ## Write unmapped values 
  missed_mappings <- extract_missing_values(appended_inputs,myinputfiles)
  missed_mappings <- lapply(missed_mappings, function(x){
    return(x %>%
             dplyr::mutate(dplyr::across(dplyr::everything(),
                                  as.character)))
  })
  missed_mappings <- dplyr::bind_rows(as.vector(missed_mappings), .id = "Data File")
  xlsx::write.xlsx(missed_mappings,
                   file = paste0("metLinkR_output/mapping_library.xlsx"),
                   row.names=FALSE,sheetName="Unmapped Metabolites",
                   append = TRUE,showNA=FALSE)

  ## Write disagreements
  xlsx::write.xlsx(do.call(rbind, initial_disagreements),
                   file = paste0("metLinkR_output/mapping_library.xlsx"),
                   row.names=FALSE,sheetName="Ambiguous Metabolites",
                   append = TRUE,showNA=FALSE)
  
  ## Write multi-mapped metabolites
  xlsx::write.xlsx(multimappings,
                   file = paste0("metLinkR_output/mapping_library.xlsx"),
                   row.names=FALSE,sheetName="MultiMapped Metabolites",
                   append = TRUE,showNA=FALSE)
  ## ## Write metadata
  ## xlsx::write.xlsx(metadata,
  ##                  file = paste0("metLinkR_output/mapping_library.xlsx"),
  ##                  row.names=FALSE,sheetName="Metabolite Metadata",
  ##                  append = TRUE,showNA=FALSE)
  
  ## Write text log
  write_txt_log(start_time,myinputfiles)

  ## Write PDF report
  names(mapping_rates[[2]]) = names(mapped_input_list) = myinputfiles$ShortFileName
  write_html_report(mapping_rates,
                   mapped_input_list,
                   mapping_library_long)

  print("(5/5) Wrote output files to metLinkR_output/")
  return(mapping_library)
}

#' Run total harmonization of metabolites
#'
#' @param inputcsv character of pathway to csv containing names of files with metabolite names
#' @param outputFileName character of output file
#' @param writecsv boolean for whether or not to write a csv file of the dataframe with harmonized files
#'
#' @return a dataframe containing harmonized files across all input files listed in inputcsv
#' @examples
#' \dontrun{
#' finalHarmonized <- runHarmonization(
#'   inputcsv = "HarmInputFiles.csv",
#'   outputFileName = "exampleharmonizedfile"
#' )
#' }
#' @export

harmonizeInputSheets <- function(inputcsv, outputFileName, writecsv) {
  ##########################################################################
  ## 1. Read in input files. Output a list of dataframes                  ##
  ##########################################################################
  myinputfiles <- utils::read.csv(inputcsv, header = T)
  list_input_files <- readInputCSVs(inputcsv)
  myinputfiles_list <- as.list(data.frame(t(myinputfiles)))
  myinputfiles_list <- mclapply(myinputfiles_list,function(x){
    out <- t(data.frame(x))
    colnames(out) <- colnames(myinputfiles)
    return(as.data.frame(out))
  })
  print("(1/5) Imported files")

  ##########################################################################
  ## 2. Initial RefMet mappings                                           ##
  ##########################################################################
  mapped_list_input_files <- mcmapply(function(x, y) {
    queryRefMet(
      input_df = x,
      filename = y$ShortFileName,
      HMDB_col = y$HMDB,
      metab_col = y$Metabolite_Name,
      CID_col = y$PubChem_CID
    )
  }, x = list_input_files,
  y = myinputfiles_list,
  SIMPLIFY = FALSE)

  refmet_mapped_ids <- mclapply(
    mapped_list_input_files,
    function(x) {
      return(x %>%
        mutate("Origin" = "Raw input") %>%
        mutate("classFlag" = "Species") %>%
        dplyr::filter(`Standardized name` != "-"))
    }
  )
  print("(2/5) Performed initial RefMet mapping")

  ################################################################################
  ## 3. Assemble synonyms table from RaMP-DB for missed IDs, perform mass check ##
  ################################################################################
  refmet_unmapped_ids <- list()
  for (i in 1:length(list_input_files)) {
    refmet_unmapped_ids[[i]] <-
      list_input_files[[i]][which(mapped_list_input_files[[i]]$`Standardized name` == "-"), ]
  }
  missed_ids <- mcmapply(function(x, y) {
    extract_identifiers(
      input_df = replaceEmptys(x),
      HMDB_col = y$HMDB,
      metab_col = y$Metabolite_Name,
      CID_col = y$PubChem_CID,
      ramp_prefixes = TRUE
    )
  }, x = refmet_unmapped_ids, y = myinputfiles_list, SIMPLIFY = FALSE)

  db <<- RaMP()
  synonym_table_list <- mclapply(missed_ids, queryRampSynonyms)
  print("(3/5) Found RaMP synonyms for unmapped inputs")
  ##########################################################################
  ## 4. Re-Query RefMet with synonym table                                ##
  ##########################################################################
  mapped_list_synonyms <- mcmapply(function(x, y) {
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
    refmet_mapped_synonyms[[i]] <- mapped_list_synonyms[[i]] %>%
      dplyr::filter(`Standardized name` != "-") %>%
      dplyr::left_join(synonym_table_list[[i]], c("Input name" = "Synonym")) %>%
      dplyr::select(-`Input name`) %>%
      dplyr::rename("Input name" = "Input") %>%
      dplyr::distinct() %>%
      mutate("Origin" = "Synonym")
  }
  print("(4/5) Queried RaMP synonyms in RefMet")

  mapped_input_list <- mcmapply(function(x, y) {
    rbind(
      x,
      y[, colnames(x)]
    )
  }, x = refmet_mapped_ids, y = refmet_mapped_synonyms, SIMPLIFY = FALSE)

  ##########################################################################
  ## 4.5 Display mapping rate summary                                     ##
  ##########################################################################
  mapping_rates <-
    calculate_mapping_rates(mapped_input_list, list_input_files, myinputfiles)

  ##########################################################################
  ## 5. Merge mapped files and write to csv                               ##
  ##########################################################################
  mapping_library <- merge_files(mapped_input_list,myinputfiles)
  
  if (writecsv) {
    utils::write.csv(mapping_library, file = paste0(outputFileName, ".csv"),
                     row.names=FALSE,na="-")
  }
  
  print(paste0("(5/5) Wrote output file: ", outputFileName,".csv"))
  return(mapping_library)
}

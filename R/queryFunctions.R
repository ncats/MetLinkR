#' Query RefMet for standardized names. Code from Eoin Fahy, received on 5/8/2023
#'
#' @param input_df dataframe of metabolites
#' @param filename name of file
#' @param HMDB_col HMDB column name
#' @param CID_col CID column name
#' @param KEGG_col KEGG column name
#' @param LM_col LM column name
#' @param CHEBI_col CHEBI column name
#' @param metab_col metabolite column name
#' @param synonym_search boolean for whether to search synonyms
#' #' @return a dataframe of RefMet query results: input/standardized name, classes, molecular formula, mass
#' @examples
#' \dontrun{
#' queryResults <- queryRefMet(mets = "2'-Deoxyuridine")
#' }
#' @export
queryRefMet <- function(input_df, filename, HMDB_col, CID_col, KEGG_col = NA,
                        LM_col = NA, CHEBI_col = NA,
                        metab_col,synonym_search=FALSE) {
  if(length(input_df)==1){
    return(NA)
  }
  input_df <- input_df %>% replaceEmptys()
  id_df <- extract_identifiers(
    input_df = input_df, HMDB_col = HMDB_col, CID_col = CID_col,
    KEGG_col = KEGG_col, LM_col, CHEBI_col, metab_col
  )
  colnames(id_df) <- c("rownum","ID","priority","origin")
  id_vector <- id_df$ID
  if(synonym_search){
    id_vector = unique(id_vector)
  }
  id_vector <- paste0(id_vector, collapse = "\n")
  res <- httr::RETRY("POST", "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_min.php",
                     body = list(metabolite_name = id_vector),
                     encode = "form",
                     times = 10
                     )
  x <- httr::content(res)
  y <- strsplit(x, "\n")
  df <- data.frame(ncol = 7)

  for (i in 1:length(y[[1]])) {
    if (nchar(y[[1]][i]) > 1) {
      z <- strsplit(y[[1]][i], "\t")
      for (j in 1:length(z[[1]])) {
        df[i, j] <- z[[1]][j]
      }
    }
  }

  df1 <- df[rowSums(is.na(df)) != ncol(df), ]
  colnames(df1) <- df1[1, ]
  df1 <- df1[-c(1), ]
  df1 <- df1 %>%
    dplyr::left_join(id_df, by = c("Input name" = "ID"))
  return(df1)
}


##' @title Query Ramp Synonyms
##' @param ids a string of metabolites separated by newline
##'
##' @return a dataframe of all RaMP synonyms found for inputs,
##' with class/species flags
##' @author Andrew Patt
queryRampSynonyms <- function(ids){
  RaMP_prefixes <-
    RaMP::getPrefixesFromAnalytes(db = db, analyteType="metabolite")[[2]]
  RaMP_prefixes <- strsplit(RaMP_prefixes,", ")[[1]]

  list_ids <- ids[grepl(paste0(RaMP_prefixes, collapse = "|"),ids)]

  list_ids <- sapply(list_ids,shQuote)
  list_ids <- paste(list_ids,collapse = ",")

  list_names <- parse_names(ids, RaMP_prefixes)

  queryId <- paste0(
    "SELECT DISTINCT source.rampId,source.sourceId,source.commonName,
     chem_props.mol_formula, chem_props.mw
     FROM source
     LEFT JOIN chem_props ON source.rampId = chem_props.ramp_id
     WHERE source.sourceId in (",list_ids,")
     OR source.commonName in (",list_names,");")
  resRampId <- RaMP::runQuery(queryId,db=db)

  if(nrow(resRampId)==0){
    dbID_synonyms = NA
  }else{
    ## Mass check
    ## checkValid_ids <- massCheck(resRampId)
    ## if("rampId" %in% rownames(checkValid_ids)){
    ##   checkValid_ids <- t(checkValid_ids) %>% as.data.frame
    ## }

    resRampIdStr <- sapply(resRampId,shQuote)
    resRampIdStr <- paste(resRampIdStr,collapse = ",")
    querywRamp <- paste0(
      "SELECT DISTINCT rampId, Synonym FROM
       analyteSynonym WHERE rampId in (",resRampIdStr, ")")
    dbID_synonyms <- RaMP::runQuery(querywRamp,db=db)
    dbID_synonyms <- dbID_synonyms %>%
      ## dplyr::left_join(checkValid_ids, by ="rampId") %>%
      ## dplyr::filter(classFlag != "Invalid") %>%
      dplyr::left_join(resRampId, by = "rampId") %>%
      dplyr::select(`Synonym`,
                    ## `classFlag`,
                    `sourceId`,`commonName`)
  }

  ## Have to re-map original inputs to Synonyms
  if(length(nrow(dbID_synonyms))!=0){
    input_vector <- c()
    for(i in 1:nrow(dbID_synonyms)){
      if(grepl(dbID_synonyms$sourceId[i],list_ids)){
        input_vector <- c(input_vector,dbID_synonyms$sourceId[i])
      }else{
        input_vector <- c(input_vector,dbID_synonyms$commonName[i])
      }
  }
    dbID_synonyms <- dbID_synonyms %>%
      dplyr::mutate(Input = input_vector) %>%
      dplyr::select(`Synonym`,
                    ##`classFlag`,
                    `Input`)
    return(dbID_synonyms)
  }else{
    return(NA)
  }
}

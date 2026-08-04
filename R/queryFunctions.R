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
  if (length(input_df) == 1 || !is.data.frame(input_df) || nrow(input_df) == 0) {
    return(NA)
  }
  if (!is.na(metab_col)) {
    if (!(metab_col %in% colnames(input_df))) {
      return(NA)
    }
    input_df <- input_df[!is.na(input_df[[metab_col]]) & nzchar(input_df[[metab_col]]), , drop = FALSE]
    if (nrow(input_df) == 0) {
      return(NA)
    }
  }
  input_df <- input_df %>% replaceEmptys()
  id_df <- extract_identifiers(
    input_df = input_df, HMDB_col = HMDB_col, CID_col = CID_col,
    KEGG_col = KEGG_col, LM_col, CHEBI_col, metab_col
  )

  colnames(id_df) <- c("rownum","ID","priority","origin")
  
  
  id_vector <- unique(id_df$ID)
  if(synonym_search){
    id_vector = unique(id_vector)
  }
  refmet_batch_size <- 500
  x <- seq_along(id_vector)
  id_vector_list <- split(id_vector, ceiling(x / refmet_batch_size))
  df_list <- lapply(id_vector_list, send_id_vector_to_RefMet)
  df <- do.call(rbind, df_list)
  
  df1 <- df[rowSums(is.na(df)) != ncol(df), ]
  colnames(df1) <- df1[1, ]
  df1 <- df1[-c(1), ]
  df1 <- df1 %>%
    dplyr::left_join(id_df, by = c("Input name" = "ID"))
  return(df1)
}


##' @title send_id_vector_to_RefMet
##' @param id_vector ids to send to refMet (vector) 
##' @return Dataframe of RefMet results
##' @author Patt
send_id_vector_to_RefMet <- function(id_vector){
  id_vector <- paste0(id_vector, collapse = "\n")
  res <- httr::RETRY(
    "POST",
    "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_min.php",
    body = list(metabolite_name = id_vector),
    encode = "form",
    times = 10,
    terminate_on = c(400, 401, 403, 404),
    httr::timeout(60)
  )
  if (httr::http_error(res)) {
    stop(
      sprintf(
        "RefMet request failed with HTTP status %s",
        httr::status_code(res)
      ),
      call. = FALSE
    )
  }
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
  return(df)
}

##' @title refmet_metadata
##' @param metabolite RefMet metabolite name
##' @return metadata
##' @author Eoin Fahy
refmet_metadata <- function(metabolite){
  metabolite_enc <- URLencode(metabolite)
  myurl <-paste0("https://www.metabolomicsworkbench.org/data/metstatR.php?REFMET_NAME=",metabolite_enc)
  h <- curl::new_handle()
  req <- curl::curl_fetch_memory(myurl, handle = h)
  metadata <- read.table(text = rawToChar(req$content), header = TRUE, na.strings = "-", stringsAsFactors = FALSE, quote = "", comment.char = "", sep="\t");
  metadata
}


##' @title Query Ramp Synonyms
##' @param ids a string of metabolites separated by newline
##' @param db connection to RaMP sqlite database
##' @param use_metabolon_parsers Strip special characters typical of metabolon inputs
##' @return a dataframe of all RaMP synonyms found for inputs,
##' with class/species flags
##' @importFrom rlang .data
##' @author Andrew Patt
queryRampSynonyms <- function(ids, db = RaMP::RaMP(), use_metabolon_parsers = TRUE){
  RaMP_prefixes <-
    RaMP::getPrefixesFromAnalytes(db = db, analyteType="metabolite")[[2]]
  RaMP_prefixes <- strsplit(RaMP_prefixes,", ")[[1]]

  list_ids <- ids[grepl(paste0(RaMP_prefixes, collapse = "|"),ids)]

  list_ids <- sapply(list_ids,shQuote)
  list_ids <- paste(list_ids,collapse = ",")

  if(use_metabolon_parsers){
    list_names <- parse_metabolon_names(ids, RaMP_prefixes)
  }else{
    list_names <- sapply(ids,shQuote)
    list_names <- paste(list_names,collapse = ",")
  }

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

    resRampIdStr <- sapply(unique(resRampId$rampId), shQuote)
    resRampIdStr <- paste(resRampIdStr,collapse = ",")
    querywRamp <- paste0(
      "SELECT DISTINCT rampId, Synonym FROM
       analyteSynonym WHERE rampId in (",resRampIdStr, ")")
    dbID_synonyms <- RaMP::runQuery(querywRamp,db=db)
    dbID_synonyms <- dbID_synonyms %>%
      ## dplyr::left_join(checkValid_ids, by ="rampId") %>%
      ## dplyr::filter(classFlag != "Invalid") %>%
      dplyr::left_join(resRampId, by = "rampId") %>%
      dplyr::select(.data$`Synonym`,
                    ## `classFlag`,
                    .data$`sourceId`,.data$`commonName`)
  }

  ## Have to re-map original inputs to Synonyms
  if(length(nrow(dbID_synonyms))!=0){
    dbID_synonyms <- dbID_synonyms %>%
      dplyr::mutate(
        Input = dplyr::if_else(
          !is.na(.data$sourceId) & .data$sourceId != "" &
            .data$sourceId %in% ids,
          .data$sourceId,
          .data$commonName,
          missing = .data$commonName
        )
      ) %>%
      dplyr::select(.data$`Synonym`,
                    ##`classFlag`,
                    .data$`Input`) %>%
      dplyr::filter(!is.na(.data$`Synonym`) & !is.na(.data$`Input`)) %>%
      dplyr::filter(.data$`Synonym` != "" & .data$`Input` != "") %>%
      dplyr::distinct()
    return(dbID_synonyms)
  }else{
    return(NA)
  }
}


##' @title getChemProps
##' @param db db connection
##' @param nameList input vector of metabolite names
##' @return 
##' @author Patt
getChemPropsFromName <- function(db, nameList) {
  strNameList <- formatListAsString(nameList)

  query <- paste0("select * from chem_props where ramp_id in (
    select rampId
        from analytesynonym
    where Synonym in (", strNameList, ")
    union
    select rampId
        from source
    where commonName in (", strNameList, "))")

  return (RaMP::runQuery(db = db, sql = query))

}

##' @importFrom rlang .data
makeStereochemFlags <- function(metab_df, db = RaMP::RaMP()){
  metab_df <- metab_df %>%
    dplyr::mutate(RaMP_input = ifelse(.data$origin=="common name",
           .data$`Input name`,
           paste(.data$origin, .data$`Input name`,sep=":")))

  ## Query names
  metab_names <- metab_df %>%
    dplyr::filter(.data$origin=="common name") %>%
    dplyr::select(.data$`RaMP_input`)
  if(nrow(metab_names)==0){
    inchis_name <- as.data.frame(matrix(ncol=2,nrow=0)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    colnames(inchis_name) = c("common_name", "inchi_name")
    ambiguous_IDs <- NULL
    unambiguous_IDs <- NULL
  }else{
    output_list <- sapply(metab_names,shQuote)
    strNameList <- paste(output_list,collapse = ",")
    
    query <- paste0("select * from chem_props 
    left join analytesynonym
    on chem_props.ramp_id = analytesynonym.rampId
    and Synonym in (", strNameList, ")
    left join source
    on source.rampId = analytesynonym.rampId
    and commonName in (", strNameList, ")")
    
    inchis_name  <- RaMP::runQuery(db = db, sql = query) %>%
      dplyr::select(c("common_name","Synonym", "inchi")) %>%
      dplyr::rename("inchi_name"="inchi") %>%
      tidyr::pivot_longer(-.data$inchi_name) %>%
      dplyr::select(-c("name")) %>%
      dplyr::rename("common_name"="value") %>%
      unique
    
    unambiguous_names <- names(which(table(inchis_name$common_name)==1))
    ambiguous_names <- names(which(table(inchis_name$common_name)>1))
    inchis_name <- inchis_name %>%
      dplyr::filter(.data$common_name %in% unambiguous_names) %>%
      dplyr::filter(.data$common_name %in% metab_df$RaMP_input)
  }
  ## Query IDs
  metab_ids <- metab_df %>%
    dplyr::filter(.data$origin!="common name") %>%
    dplyr::select(.data$`RaMP_input`)
  if(nrow(metab_ids)==0){
    inchis_ID <- as.data.frame(matrix(ncol=2,nrow=0)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    colnames(inchis_ID) = c("sourceId", "inchi_ID")
    ambiguous_IDs <- NULL
    unambiguous_IDs <- NULL
  }else{
    output_list <- sapply(metab_ids,shQuote)
    strNameList <- paste(output_list,collapse = ",")
    
    query <- paste0("select * from chem_props 
    left join source
    on ramp_id = rampId
    and sourceId in (", strNameList, ")")
    inchis_ID  <- RaMP::runQuery(db = db, sql = query) %>%
      dplyr::select(c("sourceId","inchi"))  %>%
      dplyr::rename("inchi_ID"="inchi") %>%
      unique
    unambiguous_IDs <- names(which(table(inchis_ID$sourceId)==1))
    ambiguous_IDs <- names(which(table(inchis_ID$sourceId)>1))
    inchis_ID <- inchis_ID %>%
      dplyr::filter(.data$sourceId %in% unambiguous_IDs)
  }

  ## Merge and output
  metab_df <- metab_df %>%
    dplyr::left_join(inchis_name, by = c("RaMP_input" = "common_name")) %>%
    dplyr::left_join(inchis_ID, by = c("RaMP_input" = "sourceId")) %>%
    dplyr::mutate("inchi" = ifelse(is.na(.data$inchi_ID),
                                   .data$inchi_name,
                                   .data$inchi_ID)) %>%   
    dplyr::mutate("Stereochem" = ifelse(is.na(.data$inchi),"UNKNOWN",
                                        ifelse(.data$RaMP_input %in% ambiguous_IDs |
                                                 .data$RaMP_input %in% ambiguous_names,
                                               "MULTIPLE INCHIS",
                                               ifelse(grepl("b|t|m|s", .data$inchi),
                                                      "TRUE","FALSE")))) %>%
    dplyr::select(-c("inchi_name","inchi_ID","RaMP_input"))
  return(metab_df)
}

##' @importFrom rlang .data
getPathwayAndDBData <- function(metab_df, db = RaMP::RaMP()){
  metab_names <- metab_df %>% dplyr::pull(.data$`Standardized name`)

  pathways <- RaMP::getPathwayFromAnalyte(analytes = metab_names, 
                                          namesOrIds = "names") %>%
    dplyr::select(-.data$sourceIds) %>%
    dplyr::select(-.data$pathwaySource) %>%
    dplyr::group_by(.data$commonName) %>%
    ## dplyr::mutate(pathwayNames = paste0(unique(pathwayName), collapse = "||")) %>%
    dplyr::select(-.data$pathwayName) %>%
    dplyr::mutate(pathwayIds = paste0(unique(.data$pathwayId), collapse = "||")) %>%
    dplyr::select(-.data$pathwayId) %>%
    unique
  
  metab_df <- metab_df %>%
    dplyr::mutate(lowname = tolower(.data$`Standardized name`)) %>%
    dplyr::left_join(pathways, by = c("lowname" = "commonName")) %>%
    dplyr::select(-.data$lowname)
  return(metab_df)
}

#' Replace empty cells with NA
#'
#' @param input_df a dataframe with empty cells
#'
#' @return a modified dataframe where all empty cells are replaced with NA
#' @examples
#' \dontrun{
#' updatedf <- replaceEmptys(input_df = data.frame(A = c("", "B"), B = c("C", "E")))
#' }
#'
replaceEmptys <- function(input_df) {
  for (i in 1:nrow(input_df)) {
    for (j in 1:ncol(input_df)) {
      if (!is.na(input_df[i, j]) && inherits(input_df[i, j], "character") && nchar(input_df[i, j]) == 0) {
        input_df[i, j] <- NA
      }
    }
  }
  return(input_df)
}

##' @title Extract identifiers
##'
##' @param input_df a dataframe with metabolite names, ids (HMDB, CID), information
##' @param HMDB_col name of column containing HMDB IDs
##' @param CID_col name of column containing PubChem IDs
##' @param KEGG_col name of column containing KEGG IDs
##' @param LM_col name of column containing LipidMaps IDs
##' @param CHEBI_col name of column containing ChEBI IDs
##' @param ramp_prefixes boolean for whether to add RaMP prefixes to IDs
##' @param metab_col name of column containing metabolite names
##'
##' @return one id per metabolite based on preferred ID types
##' @author Andrew Patt
extract_identifiers <- function(input_df, HMDB_col, CID_col,
                                KEGG_col = NA,
                                LM_col = NA, CHEBI_col = NA,
                                metab_col,
                                ramp_prefixes = FALSE
                                ) {
  id_vector <- c()
  temp_vector <- c()
  origin_vector <- c()
  id_df <- as.data.frame(matrix(nrow=0,ncol=4))
  colnames(id_df) <- c("rownum", "ID", "priority","origin")
  for (x in 1:nrow(input_df)) {
    if (ramp_prefixes) {
      if (!is.na(HMDB_col)) {
        if(is.na(input_df[x, HMDB_col])){
          temp_vector <- c(temp_vector,NA)
          origin_vector <- c(origin_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("hmdb:", input_df[x, HMDB_col]))
          origin_vector <- c(origin_vector, "hmdb")
        }
      }
      if (!is.na(KEGG_col)) {
        if(is.na(input_df[x, KEGG_col])){
          temp_vector <- c(temp_vector,NA)
          origin_vector <- c(origin_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("kegg:", input_df[x, KEGG_col]))
          origin_vector <- c(origin_vector, "kegg")
        }
      }
      if (!is.na(LM_col)) {
        if(is.na(input_df[x, LM_col])){
          temp_vector <- c(temp_vector,NA)
          origin_vector <- c(origin_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("LIPIDMAPS:", input_df[x, LM_col]))
          origin_vector <- c(origin_vector, "LIPIDMAPS")
        }
      }
      if (!is.na(CHEBI_col)) {
        if(is.na(input_df[x, CHEBI_col])){
          temp_vector <- c(temp_vector,NA)
          origin_vector <- c(origin_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("chebi:", input_df[x, CHEBI_col]))
          origin_vector <- c(origin_vector, "chebi")
        }
      }
      if (!is.na(metab_col)) {
        temp_vector <- c(temp_vector, input_df[x, metab_col])
        origin_vector <- c(origin_vector, "common name")
      }
      if (!is.na(CID_col)) {
        if(is.na(input_df[x, CID_col])){
          temp_vector <- c(temp_vector,NA)
          origin_vector <- c(origin_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("CAS:", input_df[x, CID_col]))
          origin_vector <- c(origin_vector,"CAS")
        }
      }
    } else {
      if (!is.na(HMDB_col)) {
        temp_vector <- c(temp_vector, input_df[x, HMDB_col])
        origin_vector <- c(origin_vector, "hmdb")
      }
      if (!is.na(KEGG_col)) {
        temp_vector <- c(temp_vector, input_df[x, KEGG_col])
        origin_vector <- c(origin_vector, "kegg")
      }
      if (!is.na(LM_col)) {
        temp_vector <- c(temp_vector, input_df[x, LM_col])
        origin_vector <- c(origin_vector, "LIPIDMAPS")
      }
      if (!is.na(CHEBI_col)) {
        temp_vector <- c(temp_vector, input_df[x, CHEBI_col])
        origin_vector <- c(origin_vector, "chebi")
      }
      if (!is.na(metab_col)) {
        temp_vector <- c(temp_vector, input_df[x, metab_col])
        origin_vector <- c(origin_vector, "common name")
      }
      if (!is.na(CID_col)) {
        temp_vector <- c(temp_vector, input_df[x, CID_col])
        origin_vector <- c(origin_vector, "CAS")
      }
    }

    ## Pick the highest priority ID or return NA
    if(ramp_prefixes){
      if (length(temp_vector) == 0 | all(is.na(temp_vector))) {
        id_vector <- c(id_vector, NA)
      } else {
        ## Use first ID if multiple IDs per cell are specified
        if (grepl(";", stats::na.omit(temp_vector)[1])) {
          id_vector <- c(id_vector, strsplit(stats::na.omit(temp_vector)[1], ";")[[1]][1])
        } else {
          id_vector <- c(id_vector, stats::na.omit(temp_vector)[1])
        }
      }
    }else{
      if(any(grepl(";",temp_vector))){
        multi_index <- which(grepl(";",temp_vector))
        origin_vector <- append(origin_vector, origin_vector[multi_index],
                                after = multi_index)
        temp_vector <- unlist(strsplit(temp_vector,";"))
      }
      temp_df <- data.frame(x, temp_vector,1:length(temp_vector),origin_vector)
      id_df <- rbind(id_df,
                     temp_df)

    }
    temp_vector <- c()
    origin_vector <- c()
  }
  return(ifelse(ramp_prefixes,return(id_vector),return(id_df)))
}

##' @title Mass Check
##' @importFrom rlang .data
##' @param synonym_DF dataframe of synonyms
##'
##' @return dataframe with three possible values per rampId: species, class or invalid
##' @author Iris Pang, Andrew Patt
massCheck <- function(synonym_DF) {
  ID_flags <- sapply(unique(synonym_DF$rampId), function(x){
    sliceDF <- synonym_DF %>%
      dplyr::filter(.data$rampId==x)
    if (length(unique(sliceDF$mol_formula)) == 1) {
      return(data.frame(x,"Species"))
    }else{
      if(any(!is.na(sliceDF$mw))){
        minVal <- as.numeric(min(sliceDF$mw, na.rm = TRUE))
        if (all(stats::na.omit(sliceDF$mw) <= minVal * 1.1)) {
          return(data.frame(x,"Class"))
        }else{
          return(data.frame(x,"Invalid"))
        }
      }else{
        return(data.frame(x,"Invalid"))
      }
    }
  })
  out <- t(ID_flags)
  out <- data.frame(out)
  colnames(out) <- c("rampId","classFlag")
  out <- apply(out,2,unlist) %>% as.data.frame

  return(out)
}

##' @title Calculate Mapping Rates
##'
##' @param mapped_input_list list of mapped input files
##' @param myinputfiles list of input files
##' @param list_input_files list of input files
##'
##' @return vector of mapping rates for each input file
##' @author Andrew Patt
calculate_mapping_rates <- function(mapped_input_list, list_input_files,
                                    myinputfiles){
  mapping_rates = mapply(function(x,y) length(unique(x$`Input name`))/nrow(y),
                         x = mapped_input_list,
                         y = list_input_files,
                         SIMPLIFY = FALSE)
  global_mapping_rate =
    sum(sapply(mapped_input_list,
               function(x)
                 return(length(unique(x$`Input name`)))))/
    sum(sapply(list_input_files, function(x)
      return(nrow(x))))

  mapping_rates_str <- c()
  for(i in 1:length(mapping_rates)){
    mapping_rates_str <- c(mapping_rates_str,
                           paste0(myinputfiles$ShortFileName[i], ": ",
                                  round(mapping_rates[[i]],3) * 100,"%\n"))
  }
  mapping_rate_out <- paste0(
    crayon::bold("MetLinkR achieved the following mapping rates:\n"),
    paste0(mapping_rates_str, collapse=""),
    crayon::bold("Global mapping rate: ",
                 round(global_mapping_rate,3) * 100, "%\n"))
  cat(mapping_rate_out)
  return(list(global_mapping_rate,mapping_rates))
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n){
  substr(x, 1, n)
}

strip_prefixes <- function(id){
  return(gsub("hmdb:|kegg:|LIPIDMAPS:|pubchem:|CAS:","",id))
}

formatListAsString <- function(idList) {
  output_list <- sapply(idList,shQuote)
  output_list <- paste(output_list,collapse = ",")
  return (output_list)
}

record_disagreements <- function(refmet_df){
  disagreements <- refmet_df %>%
    makeStereochemFlags %>%
    dplyr::filter(.data$`Standardized name` != "-") %>%
    dplyr::group_by(.data$rownum) %>%
    dplyr::filter(dplyr::n_distinct(.data$`Standardized name`)>1) %>%
    dplyr::mutate(input_group_label = dplyr::cur_group_id()) %>%
    ## Put group label as first column
    dplyr::mutate(all_isomers = ifelse(dplyr::n_distinct(.data$`Formula`)>1,
                                       FALSE, TRUE)) %>%
    as.data.frame %>%
    dplyr::select(-c("priority", "rownum"))
}

filter_hits <- function(x, majority_vote){
  if(majority_vote){
    x <- x %>%
      dplyr::mutate("Origin" = "Original input") %>%
      dplyr::filter(.data$`Standardized name` != "-")

    out <- x[0,]
    for(i in unique(x$rownum)){
      temp_df <- x %>%
        dplyr::filter(.data$rownum==i)
      if(nrow(temp_df)<3){
        if(any(!is.na(temp_df$priority))){
          out <- rbind(out,
                       temp_df %>%
                         dplyr::filter(.data$priority==min(.data$priority)))
        }
      }else{
        freqs = sort(table(temp_df$`Standardized name`), decreasing = TRUE)
        if(length(which(freqs==max(freqs)))==1){
          consensus = names(which(freqs==max(freqs)))
          temp_df = temp_df %>%
            dplyr::filter(.data$`Standardized name` == consensus) %>%
            dplyr::filter(.data$`priority` == min(.data$`priority`))
          out <- rbind(out,
                       temp_df)
        }else{
          out <- rbind(out,
                       temp_df %>%
                         dplyr::filter(.data$priority==min(temp_df$priority)))
        }
      }
    }
    return(out)
  }else{
    x %>%
      dplyr::mutate("Origin" = "Original input") %>%
      dplyr::filter(.data$`Standardized name` != "-") %>%
      dplyr::group_by(.data$rownum) %>%
      dplyr::filter(.data$priority == min(.data$priority)) %>%
      as.data.frame() %>%
      dplyr::select(-c(.data$priority))
  }
}

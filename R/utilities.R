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

##' @param input_df a dataframe with metabolite names, ids (HMDB, CID), information
##' @param HMDB_col name of column containing HMDB IDs
##' @param CID_col name of column containing PubChem IDs
##' @param KEGG_col name of column containing KEGG IDs
##' @param LM_col name of column containing LipidMaps IDs
##' @param CHEBI_col
##' @param metab_col
##' @return one id per metabolite based on preferred ID types
##' @author Andrew Patt
extract_identifiers <- function(input_df, HMDB_col, CID_col,
                                KEGG_col = NA,
                                LM_col = NA, CHEBI_col = NA,
                                metab_col,
                                ramp_prefixes = FALSE) {
  id_vector <- c()
  temp_vector <- c()
  id_df <- as.data.frame(matrix(nrow=0,ncol=3))
  colnames(id_df) <- c("rownum", "ID", "priority")

  for (x in 1:nrow(input_df)) {
    if (ramp_prefixes) {
      if (!is.na(HMDB_col)) {
        if(is.na(input_df[x, HMDB_col])){
          temp_vector <- c(temp_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("hmdb:", input_df[x, HMDB_col]))
        }
      }
      if (!is.na(KEGG_col)) {
        if(is.na(input_df[x, KEGG_col])){
          temp_vector <- c(temp_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("kegg:", input_df[x, KEGG_col]))
        }
      }
      if (!is.na(LM_col)) {
        if(is.na(input_df[x, LM_col])){
          temp_vector <- c(temp_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("LIPIDMAPS:", input_df[x, LM_col]))
        }
      }
      if (!is.na(CHEBI_col)) {
        if(is.na(input_df[x, CHEBI_col])){
          temp_vector <- c(temp_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("chebi:", input_df[x, CHEBI_col]))
        }
      }
      if (!is.na(metab_col)) {
        temp_vector <- c(temp_vector, input_df[x, metab_col])
      }
      if (!is.na(CID_col)) {
        if(is.na(input_df[x, CID_col])){
          temp_vector <- c(temp_vector,NA)
        }else{
          temp_vector <- c(temp_vector, paste0("CAS:", input_df[x, CID_col]))
        }
      }
    } else {
      if (!is.na(HMDB_col)) {
        temp_vector <- c(temp_vector, input_df[x, HMDB_col])
      }
      if (!is.na(KEGG_col)) {
        temp_vector <- c(temp_vector, input_df[x, KEGG_col])
      }
      if (!is.na(LM_col)) {
        temp_vector <- c(temp_vector, input_df[x, LM_col])
      }
      if (!is.na(CHEBI_col)) {
        temp_vector <- c(temp_vector, input_df[x, CHEBI_col])
      }
      if (!is.na(metab_col)) {
        temp_vector <- c(temp_vector, input_df[x, metab_col])
      }
      if (!is.na(CID_col)) {
        temp_vector <- c(temp_vector, input_df[x, CID_col])
      }
    }

    ## Pick the highest priority ID or return NA
    if(ramp_prefixes){
      if (length(temp_vector) == 0 | all(is.na(temp_vector))) {
        id_vector <- c(id_vector, NA)
      } else {
        ## Use first ID if multiple IDs per cell are specified
        if (grepl(";", na.omit(temp_vector)[1])) {
          id_vector <- c(id_vector, strsplit(na.omit(temp_vector)[1], ";")[[1]][1])
        } else {
          id_vector <- c(id_vector, na.omit(temp_vector)[1])
        }
      }
    }else{
      if(any(grepl(";",temp_vector))){
        temp_vector <- unlist(strsplit(temp_vector,";"))
      }
      temp_df <- data.frame(x, temp_vector,1:length(temp_vector))
      id_df <- rbind(id_df,
                     temp_df)
      
    }
    temp_vector <- c()
  }
  return(ifelse(ramp_prefixes,return(id_vector),return(id_df)))
}

##' @param rampId_DF RaMP IDs for mass check
##' @return dataframe with three possible values per rampId: species, class or invalid
##' @author Iris Pang, Andrew Patt
massCheck <- function(synonym_DF) {
  ID_flags <- sapply(unique(synonym_DF$rampId), function(x){
    sliceDF <- synonym_DF %>%
      dplyr::filter(rampId==x)
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

##' @param mapped_input_list 
##' @param list_input_files 
##' @return 
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
  return(mapping_rates)
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

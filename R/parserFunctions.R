#' @title Read Input CSVs
#' @param inputcsv character of pathway to csv containing names of files with metabolite names
#' @return a list of input files
readInputCSVs <- function(inputcsv){
  myinputfiles <- utils::read.csv(inputcsv, header=T)

  list_input_files <- list()
  for (i in 1:nrow(myinputfiles)) {
    temp <- utils::read.csv(myinputfiles$FileNames[i], header=T)
    list_input_files[[i]] <- temp
  }
  names(list_input_files) <- myinputfiles$ShortFileName
  return(list_input_files)
}


##' @title Parse Names
##' @param ids A vector of metabolite ids
##' @param RaMP_prefixes A vector of RaMP prefixes
##'
##' @return A vector of metabolite names
##' @author Patt
parse_metabolon_names <- function(ids, RaMP_prefixes){
  ## Remove saturation levels and MS1 ID status from Metabolon-style names
  list_names <- ids[!grepl(paste0(RaMP_prefixes, collapse = "|"),ids)]
  list_names <- sapply(list_names, function(x){
    if(grepl("\\*",substrRight(x,1))){
      x <- substrLeft(x,nchar(x)-1)
    }
    if(grepl("(\\*)",substrRight(x,3))){
      x <- substrLeft(x, nchar(x)-3)
    }
    return(x)
  })

  ## Remove unnamed Metabolon metabolites to speed up query
  list_names <- list_names[!grepl("X - ",substrLeft(list_names,4))]

  list_names <- sapply(list_names,shQuote)
  list_names <- paste(list_names,collapse = ",")
  
  return(list_names)
}

##' @importFrom rlang .data
append_standard_names <- function(mapped_input_list, list_input_files){
  out <- mapply(function(x, y) {
    y <- y %>%
      dplyr::mutate("rownum" = 1:nrow(y))
    x <- x %>%
      dplyr::distinct(.data$rownum, .keep_all = TRUE) %>%
      dplyr::select(c(.data$`Standardized name`,.data$`rownum`))
    out <- y %>%
      dplyr::left_join(x,by="rownum") %>%
      dplyr::select(-.data$`rownum`)
    return(out)
  }, x = mapped_input_list, y = list_input_files, SIMPLIFY = FALSE)
  return(out)
}

## utils::globalVariables(".")

##' @title Merge Files
##' @param mapped_input_list Mapped input list
##' @param myinputfiles Input files
##' @return Merged files
##' @importFrom rlang :=
##' @importFrom rlang .data
merge_files <- function(mapped_input_list, myinputfiles){
  prefixes <- RaMP::getPrefixesFromAnalytes(analyteType="metabolite")$idTypes
  prefixes <- strsplit(prefixes, ", ")[[1]]
  prefixes <- paste(prefixes, collapse = ":|")

  standard_name_list <-
    sapply(mapped_input_list,
           function(x){return(x$`Standardized name`)}) %>%
    unlist %>% unique %>% sort
  merged_df <- data.frame(`Harmonized name` = standard_name_list)
  colnames(merged_df) = "Harmonized name"
  for(i in 1:length(mapped_input_list)){
    mapped_input_list[[i]]$`Input name` <-
      sapply(mapped_input_list[[i]]$`Input name`, function(x){
        x <- gsub(prefixes,"",x)
        return(x)
    })
    merged_df <- merged_df %>%
      dplyr::left_join(mapped_input_list[[i]] %>%
                         dplyr::select("Input name", "Standardized name", "Origin" ##,
                                       ##classFlag
                                       ),
                       by = c("Harmonized name" = "Standardized name")) %>%
      dplyr::rename(!!paste0("Input name (",myinputfiles$ShortFileName[i],")"):="Input name") %>%
      dplyr::rename(!!paste0("Origin (",myinputfiles$ShortFileName[i],")"):="Origin") ## %>%
      ## dplyr::rename(!!paste0("Class flag (",myinputfiles$ShortFileName[i],")"):="classFlag")
  }
  out <- merged_df %>%
    dplyr::group_by(.data$`Harmonized name`) %>%
    dplyr::summarize_all(~ paste(unique(.), collapse = ";", sep = "")) %>%
    dplyr::ungroup() %>%
    replace(.=="NA","-")
  return(out)
}

deparen_names <- function(refmet_unmapped_ids, myinputfiles_list,
                          synonym_table_list){
  deparen_names <- mapply(function(x, y) {
    if(nrow(x)==0){
      return(NULL)
    }else{
      deparen_df <- data.frame(
        Synonym = gsub("\\s*\\([^\\)]+\\)","",
                       x[,y$Metabolite_Name]),
        Input = x[,y$Metabolite_Name]
      )
      return(deparen_df)
    }
  }, x = refmet_unmapped_ids, y = myinputfiles_list, SIMPLIFY = FALSE)
  
  synonym_table_list <- mapply(function(x, y) {
    if(is.null(x)){
      return(y)
    }else{
      return(rbind(x,y))
    }
  }, x = synonym_table_list, y = deparen_names, SIMPLIFY = FALSE)
}

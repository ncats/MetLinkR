# Jun 2023
# Andy Patt & Fred Lee
# Functions required to run HarmonizedFileforShare.R

############################################################################
## Read in database                                                       ##
############################################################################

## database <- jsonlite::fromJSON("refmet_db.json")
## database <- do.call(rbind, database)
# library(DBI)
# library(RaMP)
# library(parallel)

############################################################################
## Helper functions                                                       ##
############################################################################

#' Replace empty cells with NA
#'
#' @param input_df a dataframe with empty cells
#'
#' @return a modified dataframe where all empty cells are replaced with NA
#' @examples
#' \dontrun{
#' updatedf <- replaceEmptys(input_df = data.frame(A = c("", "B"), B = c("C", "E")))
#' }
#' @export
#'

replaceEmptys <- function(input_df){
  for (i in 1:nrow(input_df)) {
    for (j in 1:ncol(input_df)) {
      if (!is.na(input_df[i, j]) && inherits(input_df[i, j],"character") && nchar(input_df[i, j]) == 0) {
        input_df[i, j] <- NA
      }
    }
  }
  return(input_df)
}

#' Query RefMet for standardized names
#'
#' @param mets a string of metabolites separated by newline
#'
#' @return a dataframe of RefMet query results: input/standardized name, classes, molecular formula, mass
#' @examples
#' \dontrun{
#' queryResults <- queryRefMet(mets = "2'-Deoxyuridine")
#' }
#' @export
#'

queryRefMet <- function(mets){
  res <- httr::RETRY("POST", "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_min.php",
                     body = list(metabolite_name =mets),
                     encode = "form",
                     times = 5)
  x<-httr::content(res)
  y<-strsplit(x,"\n")
  df<-data.frame(ncol=7)

  for (i in 1:length(y[[1]])){
    if(nchar(y[[1]][i])>1){
      z<-strsplit(y[[1]][i],"\t")
      for (j in 1:length(z[[1]])){
        df[i,j]<-z[[1]][j]
      }
    }
  }

  df1<-df[rowSums(is.na(df)) != ncol(df), ]
  colnames(df1)=df1[1,]
  df1<-df1[-c(1),]
  return(df1)
}



############################################################################
## Main function                                                          ##
############################################################################

#' Mapping metabolites from data to refmet for standardized names
#'
#' @param input_df a dataframe with metabolite names, ids (HMDB, CID), information
#' @param filename name of file with data in input_df
#' @param HMDB_col name of column containing HMDB IDs
#' @param metab_col name of column containing metabolite names
#' @param CID_col name of column containing PubChem IDs
#'
#' @return list of dataframe of metabolite and information from RefMet and list of metabolites unharmonizable
#' @examples
#' \dontrun{
#' metFiles <- mapMetabolitesToRefMet(input_df = metabolite_information,
#' filename = "Original File", HMDB_col = "HMDB_ids",
#' metab_col = "metab_name", CID_col = "pubchem")
#' }

mapMetabolitesToRefMet <- function(input_df, filename, HMDB_col,
                                   metab_col, CID_col){

  #? include this?
  #source("Credentials.R")
  # pkg.globals <- RaMP::setConnectionToRaMP(
  #   dbname = "ramp", username = "root", conpass = "",
  #   host = "localhost"
  # )

  if(!exists("setConnectionToRaMP")){
    stop("Please install and connect to RaMP")
  }

  if (file.exists(paste0(filename, ".rds"))){
    filteredMetMappings <- readRDS(paste0(filename, ".rds"))
    return(filteredMetMappings)
  }

  ## Replace empty cells with NA
  #' @importFrom magrittr %>%
  input_df <- input_df %>% replaceEmptys

  ## Rows with multiple HMDB IDs listed (mixture)
  mixture_indices <- c(which(grepl(",", input_df[, HMDB_col]))) #identify indices of spec col w commas
  ## Initialize output df
  output_df <- stats::setNames(data.frame(matrix(nrow = nrow(input_df), ncol = 4)),
                               c(paste0("file_", filename), "input", "inputtype", "refmet_name"))

  id_vector <- apply(input_df, 1, function(x){ #iterate over rows
    if(!is.na(x[HMDB_col])){
      return(x[HMDB_col]) #returns value based on available ID
    }else if(!is.na(x[metab_col])){
      return(x[metab_col])
    }else{
      return(x[CID_col])
    }
  })

  id_vector <- paste0(id_vector,collapse = "\n")

  refMetMappings <- queryRefMet(id_vector)

  noRefMet = refMetMappings %>% dplyr::filter("Standardized name"=="-"|"Standardized name" == "")

  print(paste0("RefMet found ", round((1-nrow(noRefMet)/nrow(input_df))*100,1), "% of input metabolites outright"))
  ## If there's an HMDB, report directly. If not, map HMDB, Name, then PubChem to RaMP.
  ## If not, return "unmappable"


  comNames <- noRefMet[1]
  comNamesList <- apply(comNames, 1, function(row) toString(row))
  comNamesList <- comNamesList[comNamesList != "NA"]

  idsOnly <- comNamesList[grepl("HMDB", comNamesList)]

  nameOnly <- comNamesList[!grepl("HMDB", comNamesList)]

  con <- RaMP::connectToRaMP()

  idStore <- data.frame(nrows = length(idsOnly), ncol = 2)
  colnames(idStore) <- c("commonName", "sourceId")

  nameStore <- data.frame(nrows = length(nameOnly), ncol = 2)
  colnames(nameStore) <- c("commonName", "sourceId")

  origToSyn <- data.frame()
  ## Change this for() loop into an apply() function.
  ## Have the apply function output a dataframe of transformed IDs, then do batch query in RaMP
  ## We can study getPathwayFromAnalyte() from the RaMP package as an example batch query
  if (length(idsOnly) != 0){

    rawId_check <- function(x){

      #if multiple hmdbs provided in one input instance
      if (grepl(",", x)){
        new_ids <- strsplit(x, ",")
        x <- unlist(new_ids)
        x <- lapply(x, rawId_check)
        return(x)
      }

      #remove asterisks if any
      while(nchar(x) >= 4 & substr(x, nchar(x), nchar(x)) == "*"){
        x <- substr(x, 1, nchar(x) - 1)
      }

      #configure hmdbs with the same hmdb: tag and add zeros
      if (substr(x, 1, 5) != "hmdb:"){
        if (nchar(substr(x, 5, nchar(x))) != 7){
          zerosToAdd <- 7 - nchar(substr(x, 5, nchar(x)))
          x <- gsub("HMDB", paste0("hmdb:HMDB", strrep("0", times = zerosToAdd)), x)
        } else if (nchar(substr(x, 5, nchar(x))) == 7) {
          x <- gsub("HMDB", "hmdb:HMDB", x)
        }

      }

      return (x)
    }

    result <- as.data.frame(lapply(idsOnly, rawId_check))

    if (any(grepl(",", idsOnly))) {
      whereComma <- grep(",", idsOnly)
      for (i in rev(whereComma)) {
        replic <- length(gregexpr(",", idsOnly[i])[[1]])
        newEntries <- as.character(rep(idsOnly[i], replic))
        idsOnly <- c(idsOnly[1:i], newEntries, idsOnly[(i + 1):length(idsOnly)])
      }
    }
    origToConverti <- as.data.frame(cbind(idsOnly, unlist(result)))


    colnames(origToConverti) <- c("idsOnly", "converted")

    sql <- paste0("SELECT DISTINCT sourceId, commonName FROM source WHERE sourceId in (", toString(apply(result, 1, shQuote)), ")")

    foundId <- RMariaDB::dbGetQuery(con, sql)

    indices <- match(foundId$sourceId, origToConverti$converted)
    origNames <- origToConverti$idsOnly[indices]
    toUpdate <- as.data.frame(cbind(origNames, foundId$commonName, rep(FALSE, nrow(foundId))))

    origToSyn <- rbind(origToSyn, toUpdate)
    names(origToSyn) <- c("OrigName", "Synonym", "synFlag")

    identifiedIds <- foundId$sourceId
    identifiedIds <- identifiedIds[!duplicated(identifiedIds)]

    stillNeedIds <- result[!unlist(result) %in% identifiedIds]

    origToSyn <- origToSyn[!(origToSyn$OrigName %in% stillNeedIds), ]
    origToSyn <- origToSyn[!(origToSyn$"Synonym" %in% stillNeedIds), ]

  }
  else {
    foundId <- data.frame()
    stillNeedIds <- data.frame()
  }

  if (length(nameOnly) != 0){

    rawName_check <- function(x){

      #remove asterisks at the end of inputs
      while(!is.na(x) & nchar(x) >= 4 & substr(x, nchar(x), nchar(x)) == "*"){
        x <- substr(x, 1, nchar(x) - 1)
      }
      #remove patterns (1), (2),... (9) from inputs
      if (!is.na(x) & nchar(x) >= 4 & substr(x, nchar(x) - 2, nchar(x) -2) == "("){
        x <- substr(x, 1, nchar(x) - 4)
      }

      return (x)
    }

    result <- as.data.frame(lapply(nameOnly,rawName_check))
    origToConvertn <- as.data.frame(cbind(nameOnly, unlist(result)))
    colnames(origToConvertn) <- c("nameOnly", "converted")

    sql <- paste0("SELECT DISTINCT commonName, sourceId FROM source WHERE commonName in (", toString(apply(result, 1, shQuote)), ")")

    foundName <- RMariaDB::dbGetQuery(con, sql)

    indices <- match(tolower(foundName$commonName), tolower(origToConvertn$converted))
    origNames <- origToConvertn$nameOnly[indices]
    toUpdate <- as.data.frame(cbind(origNames, foundName$sourceId, rep(FALSE, nrow(foundName))))
    names(toUpdate) <- c("OrigName", "Synonym", "synFlag")

    origToSyn <- rbind(origToSyn, toUpdate)

    identifiedNames <- foundName$commonName
    identifiedNames <- identifiedNames[!duplicated(identifiedNames)]

    stillNeedNames <- result[!tolower(unlist(result)) %in% tolower(identifiedNames)]
    stillNeedNames <- stillNeedNames[!duplicated(stillNeedNames)]



    origToSyn <- origToSyn[!(trimws(origToSyn$OrigName) %in% trimws(stillNeedNames)), ]
    origToSyn <- origToSyn[!(trimws(origToSyn$"Synonym") %in% trimws(stillNeedNames)), ]

  }
  else {
    foundName <- data.frame()
    stillNeedNames <- data.frame()
  }

  RMariaDB::dbDisconnect(con)


  con <- RaMP::connectToRaMP()

  specCheck <- function(rampId){
    checkIds <- data.frame()
    query <- paste0("select distinct ramp_id, mol_formula, mw from chem_props where ramp_id in (", toString(sapply(rampId, shQuote)), ")")
    molCheck <- RMariaDB::dbGetQuery(con, query)

    tryCatch(
      {
        if(length(unique(molCheck$mol_formula)) == 1){
          checkIds <- cbind(rampId, rep(TRUE, nrow(rampId)))
        } else {
          if (length(unique(molCheck$"ramp_id")) == 1){
            checkIds <- cbind(rampId, rep(FALSE, nrow(rampId)))
          } else {
            tempmolCheck <- stats::na.omit(molCheck)
            minVal <- as.numeric(min(tempmolCheck$mw))
            if(all(tempmolCheck$mw <= minVal*1.1)){
              browser()
              checkIds <- cbind(rampId, rep(TRUE, nrow(rampId)))
            }
            molFormulas <- molCheck %>% dplyr::group_by("ramp_id") %>% dplyr::summarize_all(~ paste(unique(.), collapse = '||', sep = ""))
            if(any(!grepl("\\|\\|", molFormulas$mol_formula))){
              uniqueEntry <- molFormulas[!grepl("\\|\\|", molFormulas$mol_formula),]
              if(length(unique(uniqueEntry$mol_formula)) != 1){
                checkIds <- cbind(rampId, rep(FALSE, nrow(rampId)))
              } else {
                validRampIds <- uniqueEntry$"ramp_id"
                checkIds <- data.frame(validRampIds, rep(TRUE, length(validRampIds)))
                notValid <- molFormulas$"ramp_id"[!molFormulas$"ramp_id" %in% validRampIds]
                notValidIds <- data.frame(notValid, rep(FALSE, length(notValidIds)))
                checkIds <- rbind(checkIds, notValidIds)
              }
            }
          }
        }
      }, error = function(err){
        browser()
      }
    )

    colnames(checkIds) <- c("rampId", "valid")
    return (checkIds)
  }

  RMariaDB::dbDisconnect(con)

  #getting synonyms for queries
  con <- RaMP::connectToRaMP()

  synNamefromId <- character(0)
  checkValid <- data.frame()
  classCheck <- data.frame()

  allFound <- rbind(foundId, foundName)
  #Get synonyms for RaMP-mapped HMDB ids
  for (id in foundId$sourceId){
    queryId <- paste0("SELECT DISTINCT rampId FROM source WHERE sourceId = \"", id, "\"")
    resRampId <- RMariaDB::dbGetQuery(con, queryId)

    checkValid <- specCheck(resRampId)

    querywRamp <- paste0("SELECT DISTINCT rampId, Synonym FROM analyteSynonym WHERE rampId in (", toString(sapply(resRampId, shQuote)), ")")
    resRampIdi <- RMariaDB::dbGetQuery(con, querywRamp)

    int <- dplyr::full_join(resRampIdi, checkValid, by = "rampId", relationship = "many-to-many")
    classCheck <- rbind(classCheck, int)

    resRampIdi <- as.data.frame(resRampIdi$"Synonym")

    origId <- origToConverti$idsOnly[which(origToConverti$converted == id)]
    toAdd <-do.call(rbind, sapply(origId, function(id) {
      append <- list()
      newRows <- apply(resRampIdi, 1, function(syn) {
        c(id, syn, FALSE)
      })
      for (i in 1:nrow(resRampIdi)) {
        append <- c(append, list(newRows[, i]))
      }
      append
    }))

    toAdd <- as.data.frame(toAdd)
    colnames(toAdd) <- colnames(origToSyn)

    origToSyn <- rbind(origToSyn, toAdd)

    origToSyn <- origToSyn[!duplicated(origToSyn),]

    synNamefromId <- c(synNamefromId, resRampIdi)

  }

  ## Get synonyms for RaMP-mapped common names
  for (id in foundName$commonName){


    queryName <- paste0("SELECT DISTINCT rampId FROM source WHERE commonName =\"", id, "\"")
    resRampId <- RMariaDB::dbGetQuery(con, queryName)

    checkValid <- specCheck(resRampId)

    # resRampId <- checkValid$rampId[which(checkValid$valid)]
    # int <- data.frame(checkValid$rampId[which(!checkValid$valid)], rep(id, length(checkValid$rampId[which(!checkValid$valid)])))
    # classLvl <- rbind(classLvl, int)
    querywRampn <- paste0("SELECT DISTINCT rampId, Synonym FROM analyteSynonym WHERE rampId in (", toString(sapply(resRampId, shQuote)), ")")
    resRampIdn <- RMariaDB::dbGetQuery(con, querywRampn)
    resRampIdn <- resRampIdn[resRampIdn$"Synonym" != id, ]
    resRampIdn <- resRampIdn %>% stats::na.omit

    if (nrow(resRampIdn) != 0){
      int <- dplyr::full_join(resRampIdn, checkValid, by = "rampId", relationship = "many-to-many")
      classCheck <- rbind(classCheck, int)
    }

    resRampIdn <- as.data.frame(resRampIdn$"Synonym")

    origName <- origToConvertn$nameOnly[which(tolower(origToConvertn$converted) == tolower(id))]

    if (nrow(resRampIdn) > 0){

      toAdd <- do.call(rbind, lapply(origName, function(name) {
        newRows <- do.call(rbind, lapply(resRampIdn, function(syn) {
          data.frame(name = name, syn = syn, valid = FALSE)
        }))
        newRows
      }))

      toAdd <- as.data.frame(toAdd)

      colnames(toAdd) <- colnames(origToSyn)

      origToSyn <- rbind(origToSyn, toAdd)

      origToSyn <- origToSyn[!duplicated(origToSyn),]
    }
    synNamefromId <- c(synNamefromId, resRampIdn)

  }
  fullStillNeed <- c(stillNeedIds, stillNeedNames)
  fullStillNeed <- fullStillNeed[!duplicated(fullStillNeed)]

  ## Get all synonyms from unmapped names
  for (mets in fullStillNeed){


    querymet <- paste0("SELECT DISTINCT rampId FROM analyteSynonym where Synonym = \"", mets, "\"")
    resmet <- RMariaDB::dbGetQuery(con, querymet)

    if (length(row.names(resmet)) != 0){

      checkValid <- specCheck(resmet)

      queryForSyn <- paste0("SELECT DISTINCT rampId, Synonym FROM analyteSynonym WHERE rampId in (", toString(sapply(resmet, shQuote)), ")")
      resSyn <- RMariaDB::dbGetQuery(con, queryForSyn)
      resSyn <- resSyn[resSyn$"Synonym" != mets,]

      int <- dplyr::full_join(resSyn, checkValid, by = "rampId", relationship = "many-to-many")
      classCheck <- rbind(classCheck, int)

      resSyn <- as.data.frame(resSyn$"Synonym")
      resSyn <- data.frame(mets = unlist(resSyn))


      origMetsCol <- rep(mets, nrow(resSyn))
      flagCol <- rep(TRUE, nrow(resSyn))
      toAttach <- cbind(origMetsCol, resSyn, flagCol)
      names(toAttach) <- c("OrigName", "Synonym", "synFlag")

      origToSyn <- rbind(origToSyn, toAttach)

      synNamefromId <- c(synNamefromId, resSyn)
      fullStillNeed <- fullStillNeed[-which(fullStillNeed == mets)]

    }
  }

  RMariaDB::dbDisconnect(con)

  filteredMetMappings <- refMetMappings %>% dplyr::filter("Standardized name" != "-")
  `Original Input` <- filteredMetMappings$`Input name`
  `Synonym Direct Use` <- rep(FALSE, nrow(filteredMetMappings))
  `Class or Species Based` <- rep("NA", nrow(filteredMetMappings))
  filteredMetMappings <- cbind(`Original Input`, filteredMetMappings, `Synonym Direct Use`, `Class or Species Based`)
  ## Follow up refMet query for missed metabolites
  for (grp in 1:length(synNamefromId)){
    inspectSynList <- unlist(synNamefromId[grp])
    inspectSynList <- paste0(inspectSynList, collapse = "\n")

    newrefMetMappings <- queryRefMet(inspectSynList)
    newrefMetMappings <- newrefMetMappings %>% dplyr::filter ("Standardized name" != "-")


    if(nrow(newrefMetMappings) != 0){
      #browser()
      #`Original Input` <- origToSyn$OrigName[origToSyn$Synonym %in% newrefMetMappings$`Input name`]
      classChecksub <- subset(classCheck, select = c("Synonym", "valid"))
      int <- classChecksub[(classChecksub$"Synonym" %in% newrefMetMappings$`Input name`),]
      int$"classSpec" <- ifelse(int$"valid", "species", "class")

      classChecksub <- subset(int, select = c("Synonym", "classSpec"))

      newrefMetMappings <- dplyr::full_join(newrefMetMappings, classChecksub, by = c(`Input name` = "Synonym"), relationship = "many-to-many")

      newrefMetMappings <- dplyr::full_join(origToSyn %>% dplyr::filter("Synonym" %in% newrefMetMappings$`Input name`), newrefMetMappings, by = c("Synonym" = "Input name"), relationship = "many-to-many")
      newrefMetMappings <- newrefMetMappings[, c("OrigName", "Synonym", "Standardized name", "Formula", "Exact mass", "Super class", "Main class", " Sub class", "synFlag", "classSpec")]
      colnames(newrefMetMappings) <- colnames(filteredMetMappings)
      newrefMetMappings <- newrefMetMappings %>%
        dplyr::distinct()
    }


    if (length(row.names(newrefMetMappings)) == 0){
      fullStillNeed <- c(fullStillNeed, toString(unlist(synNamefromId[grp])[1]))
    }

    filteredMetMappings <- rbind(filteredMetMappings, newrefMetMappings)
  }
  fullStillNeed <- fullStillNeed[!duplicated(fullStillNeed)]
  filteredMetMappings <- filteredMetMappings[!duplicated(filteredMetMappings),]

  multimappings <- nrow(filteredMetMappings)-length(unique(filteredMetMappings$`Input name`))

  if(nchar(CID_col) != 0){
    pubChemCID <- input_df[unlist(input_df[,metab_col]) %in% unlist(fullStillNeed),]
    stillNeedPC <- fullStillNeed[!(fullStillNeed %in% pubChemCID[,metab_col])]

    refToOrig <- data.frame(stillNeed = NULL, origin = NULL)
    hmdbIds <- NULL
    if(any(grepl("HMDB", stillNeedPC))){
      hmdbIds <- gsub("hmdb:", "", stillNeedPC[grepl("HMDB", stillNeedPC)])
      refToOrig <- rbind(refToOrig, data.frame(stillNeed = hmdbIds, origin = unlist(stillNeedPC[which(grepl("HMDB", stillNeedPC))])))
      colnames(refToOrig) <- c("stillNeed", "origin")
      stillNeedPC <- stillNeedPC[-which(grepl("HMDB", stillNeedPC))]

    }

    originalInputs <-  origToSyn$OrigName[origToSyn$"Synonym" %in% stillNeedPC]

    refToOrig <- rbind(refToOrig, data.frame(stillNeed = origToSyn$"Synonym"[origToSyn$"Synonym" %in% stillNeedPC] , origin = originalInputs))
    stillNeedPC <- subset(stillNeedPC, !(stillNeedPC %in% origToSyn$"Synonym"))

    nameCheck <- origToConvertn[origToConvertn$converted %in% stillNeedPC,]

    originalInputs <- c(originalInputs, nameCheck$nameOnly[nameCheck$converted %in% stillNeedPC])
    hmdbOrigInputs <- c(hmdbIds,originalInputs[grepl("HMDB", originalInputs)])
    metOrigInputs <- originalInputs[!grepl("HMDB",originalInputs)]

    refToOrig <- rbind(refToOrig, data.frame(stillNeed = nameCheck$converted, origin = nameCheck$nameOnly))

    hmdbCols <- names(input_df)[which(grepl("HMDB|hmdb", names(input_df)))][1]
    toAddhmdb <- input_df[unlist(input_df[,hmdbCols]) %in% unlist(hmdbOrigInputs),]
    toAddmet <- input_df[unlist(input_df[,metab_col]) %in% unlist(metOrigInputs),]

    pubChemCID <- rbind(pubChemCID, toAddhmdb)
    pubChemCID <- rbind(pubChemCID, toAddmet)

    if(any(grepl("pubchem:", pubChemCID[,CID_col]))){
      pubChemCID[,CID_col] <- gsub("pubchem:", "", pubChemCID[,CID_col])
    }
    cid <- paste0(pubChemCID[,CID_col], collapse = "\n")

    refMetResults <- queryRefMet(cid)
    refMetResults <- refMetResults %>% dplyr::filter("Standardized name" != "-")

    if(nrow(refMetResults) != 0 & nrow(pubChemCID) != 0){
      toRef <- as.data.frame(cbind(pubChemCID[,metab_col], pubChemCID[,CID_col]))
      colnames(toRef) <- c("Metabolite Name", "CID")
      toAdd <- toRef[(toRef$`CID` %in% refMetResults$`Input name`),]

      `Original Input` <- refMetResults$`Input name`
      `Synonym Direct Use` <- rep("NA", nrow(refMetResults))
      `Class or Species Based` <- rep("NA", nrow(refMetResults))
      refMetResults <- cbind(`Original Input`, refMetResults, `Synonym Direct Use`, `Class or Species Based`)

      filteredMetMappings <- rbind(filteredMetMappings, refMetResults)


      for (mets in toAdd$`Metabolite Name`){
        if (mets %in% fullStillNeed){
          fullStillNeed <- fullStillNeed[-which(fullStillNeed == mets)]
        }else{

          int <- pubChemCID[pubChemCID[,CID_col] == toAdd$CID[toAdd$`Metabolite Name` == mets],]
          hmdbToSee <- stats::na.omit(int[,hmdbCols])[1]
          if (!is.na(hmdbToSee)){
            metToSee <- hmdbToSee
          }else{
            int <- pubChemCID[pubChemCID[,CID_col] == toAdd$CID[toAdd$`Metabolite Name` == mets],]
            metToSee <- stats::na.omit(int[,metab_col])[1]
          }
          metToEx <- refToOrig[refToOrig$origin %in% metToSee,]
          fullStillNeed <- fullStillNeed[-which(fullStillNeed == metToEx$stillNeed)]
        }

      }

    }
  }

  # print(paste0("Achieved ", round(100*(nrow(filteredMetMappings)-multimappings)/nrow(input_df),1),
  #              "% mapping with RaMP-DB"))

  print(paste0("Achieved ", round(100*(1-(length(fullStillNeed))/nrow(input_df)),1),
               "% mapping with RaMP-DB"))
  print(paste0("Found ", multimappings, " multi-mapped metabolites"))

  wNotFound <- list(filteredMetMappings, fullStillNeed)
  saveRDS(wNotFound, file = paste0(filename, ".rds"))
  return(wNotFound)
}



##################################################
# function to harmonize
# Idea: take refmet_name column of each file and output to a new dataframe (1st column)
# if it does not already exist, add new column
# if exist then skip

# once all unique refmet_names have a row (different data frames)
# run thru each file to retrieve input name for that particular refmet output to detail in final harmonized file
##################################################

# pseudocode
# paste first file directly over to final harmonized file as first column
# loop thru second file's refmet names & compare with first column of final harmonized file
# if the same, paste what was queried from second file (the HMDB, metabolite name or CID) AND type of identifier
# column bind with exisitng final harmonized file

# similar process for rest of the files


utils::globalVariables(".")
#' Harmonizing metabolites across multiple files
#'
#' @param fileList list of dataframes of metabolites with RefMet query results
#' @param filterOnlyOne boolean for whether to include or exclude metabolites mapped in only one file
#'
#' @return dataframe of metabolites harmonized across all files

harmonizefiles <- function(fileList, filterOnlyOne) {

  nothingdf <- data.frame(matrix(ncol = 1, nrow = 0))
  colnames(nothingdf) <- c("HarmonizedName") #naming empty data frame
  #fileList <- list(filetoharmonize1, filetoharmonize2, filetoharmonize3, filetoharmonize4, filetoharmonize5)

  ## Fill nothingdf with all unique standardized names across inputs
  for (file in fileList) {
    unique_vals <- file$`Standardized name`[!(stringr::fixed(file$`Standardized name`) %in% nothingdf$HarmonizedName)]
    unique_vals <- unique_vals[!duplicated(unique_vals)]

    if (length(unique_vals) > 0){
      nothingdf <- rbind(nothingdf, data.frame("HarmonizedName" = unique_vals))
    }
  }
  rm()

  results_list <- vector("list", length(fileList))

  for (i in seq_along(fileList)) {
    file <- fileList[[i]]
    matching_rows <- file[!is.na(file$`Standardized name`) & file$`Standardized name` %in% nothingdf$HarmonizedName, ]
    norampuse <- matching_rows[matching_rows$`Original Input` == matching_rows$`Input name`,]
    rampuse <- matching_rows[matching_rows$`Original Input` != matching_rows$`Input name`,]

    inputs <- data.frame(matching_rows$`Original Input`, matching_rows$`Standardized name`)
    names(inputs) <- c("OrigInput", "HarmonizedName")

    mapped_inputs <- sapply(nothingdf$"HarmonizedName", function(x){
      inputFiltered <- inputs %>% dplyr::filter("HarmonizedName"==x)
      uniqueValues <- unique(inputFiltered$OrigInput)
      toReturn <- paste0(uniqueValues,collapse = ";")
      return (toReturn)
    })

    intermediatedf<-cbind(nothingdf,mapped_inputs)
    intermediatedf <- replace(intermediatedf, is.na(intermediatedf), "-")

    inputsToEx <- intermediatedf$mapped_inputs


    colOfOrigin <- ifelse(grepl(";", inputsToEx), "toExamine",
                          ifelse(!is.na(inputsToEx) & inputsToEx %in% norampuse$`Original Input`,
                                 ifelse(grepl("HMDB", inputsToEx), "fromHMDB", "fromMetaboliteName"),
                                 ifelse(!is.na(inputsToEx) & inputsToEx %in% rampuse$`Original Input`,
                                        ifelse(grepl("HMDB", inputsToEx), "fromrampusingHMDB", "fromrampusingMetaboliteName"), "-")))


    toAdd <- as.data.frame(cbind(intermediatedf$mapped_inputs, colOfOrigin))
    toExamine <- toAdd %>% dplyr::filter(toAdd$colOfOrigin == "toExamine")
    toExamineIndices <- which(toAdd$colOfOrigin == "toExamine")

    newCol <- apply(toExamine, 1, function (x){
      inputsToEx <- unlist(strsplit(x[1], ";"))
      origToLookAt <- ifelse(!is.na(inputsToEx) & inputsToEx %in% norampuse$`Original Input`,
                             ifelse(grepl("HMDB", inputsToEx), "fromHMDB", "fromMetaboliteName"),
                             ifelse(!is.na(inputsToEx) & inputsToEx %in% rampuse$`Original Input`,
                                    ifelse(grepl("HMDB", inputsToEx), "fromrampusingHMDB", "fromrampusingMetaboliteName"), "-"))
      toAdd <- paste(unique(origToLookAt), collapse = ";")
      return(toAdd)

    })

    toAdd$colOfOrigin[toExamineIndices] <- newCol

    names(toAdd) <- c(paste0("file", i),paste0("origin_file", i))

    nothingdf <- cbind(intermediatedf, toAdd)
    nothingdf <- nothingdf %>% dplyr::select(-"mapped_inputs")
  }
  nothingdf <- nothingdf[!duplicated(nothingdf), ]
  nothingdf <- replace(nothingdf,nothingdf=="","-")

  if (filterOnlyOne){
    nothingdf <- nothingdf %>% dplyr::filter(rowSums(. == "-") != 2*(length(fileList) - 1))
    return(nothingdf)
  }else{
    return(nothingdf)
  }

}


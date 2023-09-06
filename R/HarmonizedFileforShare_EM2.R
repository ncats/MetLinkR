# Jan 2023
# Fred Lee
# Script to perform metabolite naming harmonization on multiple files at once
# Workflow:
# map to Refmet (priority hierarchy: HMDB, metaboltie name, CID *if included in file)
# if not RefMet hit, use HMDB to link to a metabolite name using RaMP; query RefMet with said metabolite name
# otherwise defined as "unharmonizable"


##################################################
# Set up environment (Load in required packages, remove existing data, load in relevant functions)
##################################################
# library(magrittr)
# library(jsonlite)
# library(rlang)
# library(tidyverse)
# library(httr)
# rm(list=ls())
# ## source("HarmonizeFunctions_EM2.R")
# source("HarmonizeFunctions_ACP.R")

##################################################
# Read in Input Files
# Each input file is required to have the following:
# - Columns with metabolite names, and/or HMDB, and/or PubChem CID
# - Metabolites are in rows, descriptions are in columns
# - HMDB and PubChem CIDs should NOT be prepended with anything (e.g. "HMDB000001" is allowed, "hmdb:HMDB0000001" is not)
# - All files should be CSVs
##################################################

#' Run total harmonization of metabolites
#'
#' @param inputcsv character of pathway to csv containing names of files with metabolite names
#' @param outputFileName character of output file
#' @param writecsv boolean for whether or not to write a csv file of the dataframe with harmonized files
#'
#' @return a dataframe containing harmonized files across all input files listed in inputcsv
#' @examples
#' \dontrun{
#' finalHarmonized <- runHarmonization(inputcsv = "HarmInputFiles.csv",
#' outputFileName = "exampleharmonizedfile")
#' }
#' @export

runHarmonization <- function(inputcsv, outputFileName, writecsv){
  myinputfiles <- utils::read.csv(inputcsv, header=T)

  list_input_files <- list()
  for (i in 1:nrow(myinputfiles)) {
    temp <- utils::read.csv(myinputfiles$FileNames[i], header=T)
    list_input_files[[i]] <- temp
  }
  names(list_input_files) <- myinputfiles$ShortFileName

  mapped_list_input_files <- list()

  for (i in 1:length(list_input_files)) {
    mapped_list_input_files[[i]] <- mapMetabolitesToRefMet(input_df=list_input_files[[i]],
                                                           filename=names(list_input_files)[i],
                                                           HMDB_col=myinputfiles$HMDB[i],
                                                           metab_col=myinputfiles$Metabolite_Name[i],
                                                           CID_col=myinputfiles$PubChem_CID[i]
    )
  }

  mapped_list_input_filesmets <- lapply(mapped_list_input_files, function(x) x[1])
  names(mapped_list_input_filesmets)<-myinputfiles$ShortFileName

  outputFileName <- harmonizefiles(fileList = mapped_list_input_filesmets,
                                   filterOnlyOne = FALSE)

  row.names(outputFileName) <- NULL
  outputFileName <- outputFileName[!duplicated(outputFileName),]
  browser()
  if (writecsv){
    utils::write.csv(outputFileName, paste0(outputFileName, ".csv"))
  }
  return(outputFileName)

}

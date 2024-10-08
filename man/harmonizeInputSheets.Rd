% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workflow.R
\name{harmonizeInputSheets}
\alias{harmonizeInputSheets}
\title{Run total harmonization of metabolites}
\usage{
harmonizeInputSheets(inputcsv, mapping_library_format = "both", n_cores = 1)
}
\arguments{
\item{inputcsv}{character of pathway to csv containing names of files with metabolite names}

\item{mapping_library_format}{string, print library in "wide" format, "long" format, or "both"}
}
\value{
a dataframe containing harmonized files across all input files listed in inputcsv
}
\description{
Run total harmonization of metabolites
}
\examples{
\dontrun{
finalHarmonized <- runHarmonization(
  inputcsv = "HarmInputFiles.csv"
)
}

}

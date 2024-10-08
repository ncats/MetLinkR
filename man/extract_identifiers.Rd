% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{extract_identifiers}
\alias{extract_identifiers}
\title{Extract identifiers}
\usage{
extract_identifiers(
  input_df,
  HMDB_col,
  CID_col,
  KEGG_col = NA,
  LM_col = NA,
  CHEBI_col = NA,
  metab_col,
  ramp_prefixes = FALSE
)
}
\arguments{
\item{input_df}{a dataframe with metabolite names, ids (HMDB, CID), information}

\item{HMDB_col}{name of column containing HMDB IDs}

\item{CID_col}{name of column containing PubChem IDs}

\item{KEGG_col}{name of column containing KEGG IDs}

\item{LM_col}{name of column containing LipidMaps IDs}

\item{CHEBI_col}{name of column containing ChEBI IDs}

\item{metab_col}{name of column containing metabolite names}

\item{ramp_prefixes}{boolean for whether to add RaMP prefixes to IDs}
}
\value{
one id per metabolite based on preferred ID types
}
\description{
Extract identifiers
}
\author{
Andrew Patt
}

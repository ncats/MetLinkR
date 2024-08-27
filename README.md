# MetLinkR

MetLinkR is an R package for automated cross-linking of human metabolite identities across datasets to facilitate metaanalysis. It relies on [RefMet](https://www.metabolomicsworkbench.org/databases/refmet/index.php) and [RaMP-DB](https://rampdb.nih.gov/) to standardize names and find cross-database links. MetLinkR is capable of linking common names, HMDB IDs, KEGG IDs, LIPID MAPS IDs, PubChem IDs, and ChEBI IDs for human metabolites. Our vignette (see below) contains instructions on how to format your identifiers so that metLinkR can cross-link them for you.

# Installation Instructions

You'll need to install the RaMP package with the following code:

```
# Locally install RaMP
install.packages("devtools")
library(devtools)
install_github("ncats/RAMP-DB")

# Load the package
library(RaMP)

# initializes the RaMP database object, downloading and caching the latest SQLite database
# if no version already exists in local cache.
rampDB <- RaMP()

# note that you can use the following method to check database versions hosted in your
# computer's local cache and databases that are available to download in our remote repository.
RaMP::listAvailableRaMPDbVersions()

# using that list of available RaMP DB versions, one can specify the database version to use
# if the selected version is not available on your computer, but is in our remote repository at GitHub,
# the SQLite DB file will be automatically downloaded into local file cache.
# RaMP is using the BiocFileCache package to manage a local file cache.
rampDB <- RaMP(version = "2.5.4")
```

You can install metLinkR by running:

```
devtools::install_github("ncats/MetLinkR")
```

# Citation

The manuscript for metLinkR is currently in preparation. 

# Vignette

See our vignette [here](https://github.com/ncats/MetLinkR/blob/gh-pages/Vignette.Rmd). MetLinkR is a relatively simple package to use with only one exported function. The vignette explains how to format your inputs and provides an overview of the outputs expected from a metLinkR analysis.

# Contact

For questions or support, please contact the author at <andrew.patt@nih.gov>. Please report any bugs to the Github issues for this repository.

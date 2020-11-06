[![R-CMD-check](https://github.com/timoast/signac/workflows/R-CMD-check/badge.svg)](https://github.com/timoast/signac/actions)
[![CRAN Version](https://www.r-pkg.org/badges/version/Signac)](https://cran.r-project.org/package=Signac)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Signac)](https://cran.r-project.org/package=Signac)

# Signac 1.1.0

Signac is an extension of [Seurat](https://satijalab.org/seurat) for the analysis of single-cell chromatin data.

Documentation can be found at https://satijalab.org/signac/

## Install

```r
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
setRepositories(ind=1:2)

# To install the current release
install.packages("Signac")

# To install the development version
install.packages("devtools")
devtools::install_github("timoast/signac", ref = "develop")
```

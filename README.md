[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/timoast/signac?svg=true)](https://ci.appveyor.com/project/timoast/signac)
[![CRAN Version](https://www.r-pkg.org/badges/version/Signac)](https://cran.r-project.org/package=Signac)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Signac)](https://cran.r-project.org/package=Signac)

# Signac v1.0.0

Signac is an extension of [Seurat](https://satijalab.org/seurat) for the analysis of single-cell chromatin data.

Documentation can be found at https://satijalab.org/signac/

## Install

```{r}
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

install.packages("Signac")

# To install the development version
install.packages("devtools")
devtools::install_github("timoast/signac", ref = "develop")
```

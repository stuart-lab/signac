[![Build Status](https://travis-ci.com/timoast/signac.svg?branch=master)](https://travis-ci.com/timoast/signac)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/timoast/signac?svg=true)](https://ci.appveyor.com/project/timoast/signac)

# Signac v0.1.3

Signac is an extension of [Seurat](https://satijalab.org/seurat) for the analysis of single-cell chromatin data.

Documentation can be found at https://satijalab.org/signac/

## Install

```{r}
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# Tell R to also check bioconductor when installing dependencies
setRepositories(ind=1:2)

# Install Signac
install.packages("devtools")
devtools::install_github("timoast/signac")
```

[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/timoast/signac?svg=true)](https://ci.appveyor.com/project/timoast/signac)
# Signac v0.2.0

Signac is an extension of [Seurat](https://satijalab.org/seurat) for the analysis of single-cell chromatin data.

Documentation can be found at https://satijalab.org/signac/

## Install

```{r}
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# To install Signac from GitHub
install.packages("devtools")
devtools::install_github("timoast/signac")

# To install the development version
devtools::install_github("timoast/signac", ref = "develop")
```

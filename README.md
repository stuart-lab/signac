[![Build Status](https://dev.azure.com/timstuart/timstuart/_apis/build/status/timoast.signac?branchName=master&jobName=macOS_Mojave)](https://dev.azure.com/timstuart/timstuart/_build/latest?definitionId=1&branchName=master)
[![Build Status](https://dev.azure.com/timstuart/timstuart/_apis/build/status/timoast.signac?branchName=master&jobName=Ubuntu_1804)](https://dev.azure.com/timstuart/timstuart/_build/latest?definitionId=1&branchName=master)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/timoast/signac?svg=true)](https://ci.appveyor.com/project/timoast/signac)
# Signac v0.1.6

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

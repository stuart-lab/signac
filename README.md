# Signac

[![R-CMD-check](https://github.com/timoast/signac/workflows/R-CMD-check/badge.svg)](https://github.com/timoast/signac/actions)
[![CRAN
Version](https://www.r-pkg.org/badges/version/Signac)](https://cran.r-project.org/package=Signac)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/Signac)](https://cran.r-project.org/package=Signac)

## Overview

Signac is a comprehensive R package for the analysis of single-cell
chromatin data. Signac includes functions for quality control,
normalization, dimension reduction, clustering, differential activity,
and more.

Documentation and tutorials can be found at
<https://satijalab.org/signac/>

## Installation

Signac requires that [Bioconductor](https://www.bioconductor.org/) is
installed:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
setRepositories(ind=1:2)
```

To install the latest release of Signac from CRAN:

``` r
install.packages("Signac")
```

To release the latest develop version from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("timoast/signac", ref = "develop")
```

## Release notes

For a changelog please see the [NEWS
file](https://github.com/timoast/signac/blob/master/NEWS.md), also
available on the [Signac
website](https://satijalab.org/signac/news/index.html).

## Getting help

If you encounter a bug or have a feature request, please open an
[issue](https://github.com/timoast/signac/issues).

If you would like to discuss questions related to single-cell analysis,
you can open a
[discussion](https://github.com/timoast/signac/discussions).

## Citing Signac

If you use the Signac package in your work please cite [Stuart et
al. 2020](https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1)

    @UNPUBLISHED{signac,
      title    = "Multimodal single-cell chromatin analysis with Signac",
      author   = "Stuart, Tim and Srivastava, Avi and Lareau, Caleb and Satija,
                  Rahul",
      journal  = "bioRxiv",
      pages    = "2020.11.09.373613",
      month    =  nov,
      year     =  2020,
      url      = "https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1",
      language = "en"
    }

## Related packages

-   [Seurat](https://github.com/satijalab/seurat)
-   [SeuratObject](https://github.com/mojaveazure/seurat-object)
-   [SeuratDisk](https://github.com/mojaveazure/seurat-disk)
-   [SeuratData](https://github.com/satijalab/seurat-data)
-   [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
-   [Azimuth](https://github.com/satijalab/azimuth)

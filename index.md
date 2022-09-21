[![R-CMD-check](https://github.com/stuart-lab/signac/workflows/R-CMD-check/badge.svg)](https://github.com/stuart-lab/signac/actions)
[![CRAN Version](https://www.r-pkg.org/badges/version/Signac)](https://cran.r-project.org/package=Signac)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Signac)](https://cran.r-project.org/package=Signac)

# Signac

Signac is a comprehensive R package for the analysis of single-cell
chromatin data.

## Features

Signac is designed for the analysis of single-cell chromatin data, including scATAC-seq,
single-cell targeted tagmentation methods such as [scCUT&Tag](https://www.nature.com/articles/s41467-019-09982-5)
and [scNTT-seq](https://www.biorxiv.org/content/10.1101/2022.03.08.483436v1),
and multimodal datasets that jointly measure chromatin state alongside other
modalities.

Signac currently supports the following features:

* Calling peaks
* Quantifying per-cell counts in different genomic regions
* Calculating single-cell QC metrics
* Dimensional reduction, visualization, and clustering
* Identifying cell-type-specific peaks
* Visualizing 'pseudo-bulk' coverage tracks
* Integration of multiple single-cell datasets
* Integration with single-cell RNA-seq datasets
* Sequence motif enrichment analysis
* Transcription factor footprinting analysis
* Linking peaks to correlated genes
* Parallelization through the [future](https://cran.r-project.org/package=future) package
* Seamless interface with [Seurat](https://satijalab.org/seurat), [SeuratWrappers](https://github.com/satijalab/seurat-wrappers), [SeuratDisk](https://github.com/mojaveazure/seurat-disk), and [SeuratData](https://github.com/satijalab/seurat-data) functionality
* Interoperability with [Bioconductor](https://bioconductor.org/) tools

Please see the Signac [vignettes](articles/overview.html) page for examples.

For installation instructions see the [install](articles/install.html) page.



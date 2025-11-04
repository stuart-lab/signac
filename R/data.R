#' A small example scATAC-seq dataset
#'
#' A subsetted version of 10x Genomics 10k human (hg38) PBMC scATAC-seq dataset
#'
#' @format A Seurat object with the following assays
#' \describe{
#'   \item{peaks}{A peak x cell dataset}
#'   \item{RNA}{A gene x cell dataset}
#' }
#'
#' @concept data
#' @source \url{https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard}
"atac_small"

#' A small example scATAC-seq dataset (old)
#'
#' A subsetted version of 10x Genomics 10k human (hg19) PBMC scATAC-seq dataset.
#' This object was included in Signac v1 as the small test data object.
#'
#' @format A Seurat object with the following assays
#' \describe{
#'   \item{peaks}{A peak x cell dataset}
#'   \item{bins}{A 5 kb genome bin x cell dataset}
#'   \item{RNA}{A gene x cell dataset}
#' }
#'
#' @concept data
#' @source \url{https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k}
"atac_small_old"

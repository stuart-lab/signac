#' A small example scATAC-seq dataset
#'
#' A subsetted version of 10x Genomics 10k human (hg19) PBMC scATAC-seq dataset
#'
#' @format A Seurat object with the following assays
#'
#' \describe{
#'   \item{peaks}{A peak x cell dataset}
#'   \item{bins}{A 5 kb genome bin x cell dataset}
#'   \item{RNA}{A gene x cell dataset}
#' }
#'
#' @source \url{https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k}
"atac_small"

#' Genomic blacklist regions for Human GRCh38
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_hg38"

#' Genomic blacklist regions for Human hg19
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_hg19"

#' Genomic blacklist regions for Mouse mm10
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_mm10"

#' Genomic blacklist regions for Drosophila dm3
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_dm3"

#' Genomic blacklist regions for Drosophila dm6
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_dm6"

#' Genomic blacklist regions for C. elegans ce10
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_ce10"

#' Genomic blacklist regions for C. elegans ce11
#'
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \url{https://doi.org/10.1038/s41598-019-45839-z}
"blacklist_ce11"

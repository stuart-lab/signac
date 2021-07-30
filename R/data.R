#' A small example scATAC-seq dataset
#'
#' A subsetted version of 10x Genomics 10k human (hg19) PBMC scATAC-seq dataset
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
"atac_small"

#' Genomic blacklist regions for Human GRCh38
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_hg38"

#' Unified genomic blacklist regions for Human GRCh38
#'
#' Manually curated genomic blacklist regions for the hg38 genome by Anshul
#' Kundaje and Anna Shcherbina. See
#' \url{https://www.encodeproject.org/files/ENCFF356LFX/} for a description of
#' how this blacklist was curated.
#'
#' @concept data
#' @author Anshul Kundaje
#' @author Anna Shcherbina
#' @format A GRanges object
#' @source \url{https://www.encodeproject.org/files/ENCFF356LFX/}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_hg38_unified"

#' Genomic blacklist regions for Human hg19 (0-based)
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_hg19"

#' Genomic blacklist regions for Mouse mm10 (0-based)
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_mm10"

#' Genomic blacklist regions for Drosophila dm3 (0-based)
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_dm3"

#' Genomic blacklist regions for Drosophila dm6 (0-based)
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_dm6"

#' Genomic blacklist regions for C. elegans ce10 (0-based)
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_ce10"

#' Genomic blacklist regions for C. elegans ce11 (0-based)
#' @concept data
#' @format A GRanges object
#' @source \url{https://github.com/Boyle-Lab/Blacklist}
#' @source \doi{10.1038/s41598-019-45839-z}
"blacklist_ce11"

#' Add DNA sequence motif information
#'
#' Construct a \code{\link{Motif}} object containing DNA sequence motif
#' information and add it to an existing Seurat object or ChromatinAssay.
#' If running on a Seurat object, \code{AddMotifs} will also run
#' \code{\link{RegionStats}} to compute the GC content of each peak and store
#' the results in the feature metadata.
#'
#' @param object A Seurat object or ChromatinAssay object
#' @param ... Additional arguments passed to other methods
#' @export AddMotifs
#' @rdname AddMotifs
#' @return When running on a \code{ChromatinAssay} or \code{Seurat} object,
#' returns a modified version of the input object. When running on a matrix,
#' returns a \code{Motif} object.
AddMotifs <- function(object, ...) {
  UseMethod(generic = "AddMotifs", object = object)
}

#' Quantify aggregated genome tiles
#'
#' Quantifies fragment counts per cell in fixed-size genome bins across the
#' whole genome, then removes bins with less than a desired minimum number of
#' counts in the bin, then merges adjacent tiles into a single region.
#'
#' @param object A Seurat object or ChromatinAssay object
#' @param ... Additional arguments passed to other methods
#' @export AggregateTiles
#' @rdname AggregateTiles
AggregateTiles <- function(object, ...) {
  UseMethod(generic = "AggregateTiles", object = object)
}

#' Convert objects to a ChromatinAssay
#' @param x An object to convert to class \code{\link{ChromatinAssay}}
#' @param ... Arguments passed to other methods
#' @rdname as.ChromatinAssay
#' @export as.ChromatinAssay
as.ChromatinAssay <- function(x, ...) {
  UseMethod(generic = "as.ChromatinAssay", object = x)
}

#' Compute allele frequencies per cell
#'
#' Collapses allele counts for each strand and normalize by the total number of
#' counts at each nucleotide position.
#'
#' @param object A Seurat object, Assay, or matrix
#' @param variants A character vector of informative variants to keep. For
#' example, \code{c("627G>A","709G>A","1045G>A","1793G>A")}.
#' @param ... Arguments passed to other methods
#'
#' @export
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object with a new assay
#' containing the allele frequencies for the informative variants.
AlleleFreq <- function(object, ...) {
  UseMethod(generic = "AlleleFreq", object = object)
}

#' Annotation
#'
#' Get the annotation from a ChromatinAssay
#'
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' if the annotation data is present, otherwise returns NULL
#' @rdname Annotation
#' @export Annotation
Annotation <- function(object, ...) {
  UseMethod(generic = "Annotation", object = object)
}

#' @param value A value to set. Can be NULL, to remove the current annotation
#' information, or a \code{\link[GenomicRanges]{GRanges}} object. If a
#' \code{GRanges} object is supplied and the genome information is stored in the
#' assay, the genome of the new annotations must match the genome of the assay.
#'
#' @rdname Annotation
#' @export Annotation<-
#'
"Annotation<-" <- function(object, ..., value) {
  UseMethod(generic = 'Annotation<-', object = object)
}

#' Binarize counts
#'
#' Set counts >1 to 1 in a count matrix
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname BinarizeCounts
#' @export BinarizeCounts
BinarizeCounts <- function(object, ...) {
  UseMethod(generic = "BinarizeCounts", object = object)
}

#' Call peaks
#'
#' Call peaks using MACS. Fragment files linked to the specified assay will be
#' used to call peaks. If multiple fragment files are present, all will be used
#' in a single MACS invocation. Returns the \code{.narrowPeak} MACS output as a
#' \code{GRanges} object.
#'
#' See \url{https://macs3-project.github.io/MACS/} for MACS documentation.
#'
#' If you call peaks using MACS2 please cite:
#' \doi{10.1186/gb-2008-9-9-r137}
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' @rdname CallPeaks
#' @export CallPeaks
CallPeaks <- function(object, ...) {
  UseMethod(generic = "CallPeaks", object = object)
}

#' Set and get cell barcode information for a Fragment object
#'
#' @param x A Seurat object
#' @param value A character vector of cell barcodes
#' @param ... Arguments passed to other methods
#' @export Cells<-
"Cells<-" <- function(x, ..., value) {
  UseMethod(generic = "Cells<-", object = x)
}

#' Convert between motif name and motif ID
#'
#' Converts from motif name to motif ID or vice versa. To convert common names
#' to IDs, use the \code{name} parameter. To convert IDs to common names, use
#' the \code{id} parameter.
#'
#' @param object A Seurat, ChromatinAssay, or Motif object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a character vector with the same length and order as the
#' input. Any names or IDs that were not found will be stored as \code{NA}.
#'
#' @rdname ConvertMotifID
#' @export ConvertMotifID
#'
ConvertMotifID <- function(object, ...) {
  UseMethod(generic = "ConvertMotifID", object = object)
}

#' Find most frequently observed features
#'
#' Find top features for a given assay based on total number of counts for the
#' feature. Can specify a minimum cell count, or a lower percentile
#' bound to determine the set of variable features. Running this function will
#' store the total counts and percentile rank for each feature in the feature
#' metadata for the assay. To only compute the feature metadata, without
#' changing the variable features for the assay, set \code{min.cutoff=NA}.
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname FindTopFeatures
#' @export FindTopFeatures
FindTopFeatures <- function(object, ...) {
  UseMethod(generic = "FindTopFeatures", object = object)
}

#' Transcription factor footprinting analysis
#'
#' Compute the normalized observed/expected Tn5 insertion frequency
#' for each position surrounding a set of motif instances.
#'
#' @param object A Seurat or ChromatinAssay object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname Footprint
#' @export Footprint
Footprint <- function(object, ...) {
  UseMethod(generic = "Footprint", object = object)
}

#' Get the Fragment objects
#'
#' @param ... Arguments passed to other methods
#' @return Returns a list of \code{\link{Fragment}} objects. If there are
#' no Fragment objects present, returns an empty list.
#' @rdname Fragments
#' @export Fragments
Fragments <- function(object, ...) {
  UseMethod(generic = "Fragments", object = object)
}

#' @param value A \code{\link{Fragment}} object or list of Fragment objects
#'
#' @rdname Fragments
#' @export Fragments<-
#'
"Fragments<-" <- function(object, ..., value) {
  UseMethod(generic = 'Fragments<-', object = object)
}

#' Compute Tn5 insertion bias
#'
#' Counts the Tn5 insertion frequency for each DNA hexamer.
#' @param object A Seurat or ChromatinAssay object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname InsertionBias
#' @export InsertionBias
InsertionBias <- function(object, ...) {
  UseMethod(generic = "InsertionBias", object = object)
}

#' Retrieve a motif matrix
#'
#' Get motif matrix for given assay
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname GetMotifData
#' @export GetMotifData
GetMotifData <- function(object, ...) {
  UseMethod(generic = "GetMotifData", object = object)
}

#' Get or set a motif information
#'
#' Get or set the Motif object for a Seurat object or ChromatinAssay.
#'
#' @param ... Arguments passed to other methods
#' @rdname Motifs
#' @export Motifs
Motifs <- function(object, ...) {
  UseMethod(generic = "Motifs", object = object)
}

#' @param value A \code{\link{Motif}} object
#' @rdname Motifs
#' @export Motifs<-
"Motifs<-" <- function(object, ..., value) {
  UseMethod(generic = 'Motifs<-', object = object)
}

#' Get or set links information
#'
#' Get or set the genomic link information for a Seurat object or ChromatinAssay
#'
#' @param ... Arguments passed to other methods
#' @rdname Links
#' @export Links
Links <- function(object, ...) {
  UseMethod(generic = "Links", object = object)
}

#' @param value A \code{\link[GenomicRanges]{GRanges}} object
#' @rdname Links
#' @export Links<-
"Links<-" <- function(object, ..., value) {
  UseMethod(generic = "Links<-", object = object)
}

#' Identify mitochondrial variants
#'
#' Identify mitochondrial variants present in single cells.
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @rdname IdentifyVariants
#' @export IdentifyVariants
IdentifyVariants <- function(object, ...) {
  UseMethod(generic = "IdentifyVariants", object = object)
}

#' Compute base composition information for genomic ranges
#'
#' Compute the GC content, region lengths, and dinucleotide base frequencies
#' for regions in the assay and add to the feature metadata.
#'
#' @param object A Seurat object, Assay object, or set of genomic ranges
#' @param ... Arguments passed to other methods
#' @return Returns a dataframe
#' @rdname RegionStats
#' @export RegionStats
RegionStats <- function(object, ...) {
  UseMethod(generic = "RegionStats", object = object)
}

#' Run chromVAR
#'
#' Wrapper to run \code{\link[chromVAR]{chromVAR}} on an assay with a motif
#' object present. Will return a new Seurat assay with the motif activities
#' (the deviations in chromatin accessibility across the set of regions) as
#' a new assay.
#'
#' See the chromVAR documentation for more information:
#' \url{https://greenleaflab.github.io/chromVAR/index.html}
#'
#' See the chromVAR paper: \url{https://www.nature.com/articles/nmeth.4401}
#'
#' @param object A Seurat object
#' @param genome A \code{BSgenome}, \code{DNAStringSet}, \code{FaFile}, or
#' string stating the genome build recognized by \code{getBSgenome}.
#' @param motif.matrix A peak x motif matrix. If NULL, pull the peak x motif
#' matrix from a Motif object stored in the assay.
#' @param verbose Display messages
#' @param ... Additional arguments passed to
#' \code{\link[chromVAR]{getBackgroundPeaks}}
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object with a new assay
#' @rdname RunChromVAR
#' @export RunChromVAR
RunChromVAR <- function(object, ...) {
  UseMethod(generic = "RunChromVAR", object = object)
}

#' Run singular value decomposition
#'
#' Run partial singular value decomposition using \code{\link[irlba]{irlba}}
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname RunSVD
#' @export RunSVD
RunSVD <- function(object, ...) {
  UseMethod(generic = "RunSVD", object = object)
}

#' Compute the term-frequency inverse-document-frequency
#'
#' Run term frequency inverse document frequency (TF-IDF) normalization on a
#' matrix.
#'
#' Four different TF-IDF methods are implemented. We recommend using method 1
#' (the default).
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname RunTFIDF
#' @export RunTFIDF
#' @references \url{https://en.wikipedia.org/wiki/Latent_semantic_analysis#Latent_semantic_indexing}
RunTFIDF <- function(object, ...) {
  UseMethod(generic = "RunTFIDF", object = object)
}

#' Set motif data
#'
#' Set motif matrix for given assay
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @rdname SetMotifData
#' @export SetMotifData
SetMotifData <- function(object, ...) {
  UseMethod(generic = "SetMotifData", object = object)
}

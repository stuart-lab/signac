#' Convert objects to a ChromatinAssay
#' @param x An object to convert to class \code{\link{ChromatinAssay}}
#' @param ... Arguments passed to other methods
#' @rdname as.ChromatinAssay
#' @export as.ChromatinAssay
as.ChromatinAssay <- function(x, ...) {
  UseMethod(generic = "as.ChromatinAssay", object = x)
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
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @rdname BinarizeCounts
#' @export BinarizeCounts
BinarizeCounts <- function(object, ...) {
  UseMethod(generic = "BinarizeCounts", object = object)
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
#' Find top binary features for a given assay based on total number of cells
#' containing feature. Can specify a minumum cell count, or a lower percentile
#' bound.
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[Seurat]{Seurat}} object
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
#' @return Returns a \code{\link[Seurat]{Seurat}} object
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
#' @return Returns a \code{\link[Seurat]{Seurat}} object
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

#' Run singular value decomposition
#'
#' Run partial singular value decomposition using \code{\link[irlba]{irlba}}
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[Seurat]{Seurat}} object
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
#' @return Returns a \code{\link[Seurat]{Seurat}} object
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
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @rdname SetMotifData
#' @export SetMotifData
SetMotifData <- function(object, ...) {
  UseMethod(generic = "SetMotifData", object = object)
}

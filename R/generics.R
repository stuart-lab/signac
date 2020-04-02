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

#' Motifs
#' Retrieve a Motif object
#'
#' Get the Motif object
#'
#' @param ... Arguments passed to other methods
#' @rdname Motifs
#' @export Motifs
Motifs <- function(object, ...) {
  UseMethod(generic = "Motifs", object = object)
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
#' Run term frequency inverse document frequency (TF-IDF) normalization
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @rdname RunTFIDF
#' @export RunTFIDF
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

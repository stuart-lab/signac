#' BinarizeCounts
#'
#' Set counts >1 to 1 in a count matrix
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname BinarizeCounts
#' @export BinarizeCounts
BinarizeCounts <- function(object, ...) {
  UseMethod(generic = 'BinarizeCounts', object = object)
}

#' FindTopFeatures
#'
#' Find top binary features for a given assay based on total number of cells containing feature.
#' Can specify a minumum cell count, or a lower percentile bound.
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname FindTopFeatures
#' @export FindTopFeatures
FindTopFeatures <- function(object, ...) {
  UseMethod(generic = 'FindTopFeatures', object = object)
}

#' GetMotifData
#'
#' Get motif matrix for given assay
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname GetMotifData
#' @export GetMotifData
GetMotifData <- function(object, ...) {
  UseMethod(generic = 'GetMotifData', object = object)
}

#' RunSVD
#'
#' Run singular value decomposition
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname RunSVD
#' @export RunSVD
RunSVD <- function(object, ...) {
  UseMethod(generic = 'RunSVD', object = object)
}


#' RunTFIDF
#'
#' Run term frequency inverse document frequency normalization
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname RunTFIDF
#' @export RunTFIDF
RunTFIDF <- function(object, ...) {
  UseMethod(generic = 'RunTFIDF', object = object)
}

#' SetMotifData
#'
#' Set motif matrix for given assay
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname SetMotifData
#' @export SetMotifData
SetMotifData <- function(object, ...) {
  UseMethod(generic = 'SetMotifData', object = object)
}

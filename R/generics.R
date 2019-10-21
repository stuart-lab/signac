#' Convert objects to a ChromatinAssay
#' @param x An object to convert to class \code{\link{ChromatinAssay}}
#' @param ... Arguments passed to other methods
#' @rdname as.ChromatinAssay
#' @export as.ChromatinAssay
as.ChromatinAssay <- function(x, ...) {
  UseMethod(generic = 'as.ChromatinAssay', object = x)
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
  UseMethod(generic = 'Annotation', object = object)
}

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

#' ClusterMotifs
#'
#' Cluster motifs by co-occurrence in genomic regions.
#' Computes the Jaccard similarity between motifs based on the proportion of
#' regions that they co-occur in. This is normalized for the overall frequency of the motifs.
#' Motifs are then clustered, and cluster identities stored in the Motif object meta data.
#'
#' @param object A Seurat object
#' @return Returns a Seurat object
#' @rdname ClusterMotifs
#' @export ClusterMotifs
ClusterMotifs <- function(object, ...) {
  UseMethod(generic = 'ClusterMotifs', object = object)
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
#' @return Returns a sparse matrix
#' @rdname GetMotifData
#' @export GetMotifData
GetMotifData <- function(object, ...) {
  UseMethod(generic = 'GetMotifData', object = object)
}

#' Motifs
#'
#' Get the Motif object
#'
#' @param ... Arguments passed to other methods
#' @rdname Motifs
#' @export Motifs
Motifs <- function(object, ...) {
  UseMethod(generic = 'Motifs', object = object)
}

#' Compute base composition information for genomic ranges
#'
#' Compute the GC content, region lengths, and dinucleotide base frequencies
#' for regions in the assay and add to the feature metadata.
#'
#' @param object A Seurat object, Assay object, or set of genomic ranges
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname RegionStats
#' @export RegionStats
RegionStats <- function(object, ...) {
  UseMethod(generic = 'RegionStats', object = object)
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

#' RunMotifTSNE
#'
#' Run tSNE on a Mofif object. This will project the motifs into two tSNE dimensions based on
#' a neighbor graph stored in the Motif object (commonly computed from the Jaccard similarity between motifs).
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname RunMotifTSNE
#' @export RunMotifTSNE
RunMotifTSNE <- function(object, ...) {
  UseMethod(generic = 'RunMotifTSNE', object = object)
}

#' RunMotifUMAP
#'
#' Run UMAP on a Mofif object. This will project the motifs into two UMAP dimensions based on
#' a neighbor graph stored in the Motif object (commonly computed from the Jaccard similarity between motifs).
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname RunMotifUMAP
#' @export RunMotifUMAP
RunMotifUMAP <- function(object, ...) {
  UseMethod(generic = 'RunMotifUMAP', object = object)
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

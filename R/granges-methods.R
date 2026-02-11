#' @include generics.R
#' @importFrom SeuratObject DefaultAssay
#' @importFrom GenomicRanges granges
NULL

#' Access genomic ranges for GRangesAssay objects
#'
#' Methods for accessing [GenomicRanges::GRanges()] object
#' information stored in a [GRangesAssay-class] object.
#'
#' @name granges-methods
#' @param x A [GRangesAssay-class] object
#' @param use.names Whether the names on the genomic ranges should be
#' propagated to the returned object.
#' @param use.mcols Not supported for [GRangesAssay-class] objects
#' @param ... Additional arguments
#'
#' @return Returns a [GenomicRanges::GRanges()] object
#'
#' @aliases granges granges,GRangesAssay-method
#' @seealso
#' \itemize{
#'   \item{[granges][GenomicRanges::granges] in the \pkg{GenomicRanges} package.}
#'   \item{[GRangesAssay-class]}
#'  }
#' @exportMethod granges
#' @concept granges
#' @examples
#' granges(atac_small)
setMethod(
  f = "granges",
  signature = "GRangesAssay",
  definition = function(x, use.names = TRUE, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for GRangesAssay objects ",
           "does not support the 'use.mcols' argument")
    }
    slot(object = x, name = "ranges")
  }
)

#' @describeIn granges-methods method for Seurat objects
#' @concept granges
setMethod(
  f = "granges",
  signature = "Seurat",
  definition = function(x, use.names = TRUE, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for Seurat objects ",
           "does not support the 'use.mcols' argument")
    }
    assay <- DefaultAssay(object = x)
    granges(x = x[[assay]])
  }
)

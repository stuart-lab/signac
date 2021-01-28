#' @include generics.R
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomicRanges granges
NULL

#' Access genomic ranges for ChromatinAssay objects
#'
#' Methods for accessing \code{\link[GenomicRanges]{GRanges}} object
#' information stored in a \code{\link{ChromatinAssay}} object.
#'
#' @name granges-methods
#' @param x A \code{\link{ChromatinAssay}} object
#' @param use.names Whether the names on the genomic ranges should be
#' propagated to the returned object.
#' @param use.mcols Not supported for \code{\link{ChromatinAssay}} objects
#' @param ... Additional arguments
#'
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#'
#' @aliases granges granges,ChromatinAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[GenomicRanges]{granges} in the \pkg{GenomicRanges} package.}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod granges
#' @concept granges
#' @examples
#' granges(atac_small)
setMethod(
  f = "granges",
  signature = "ChromatinAssay",
  definition = function(x, use.names = TRUE, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for ChromatinAssay objects ",
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

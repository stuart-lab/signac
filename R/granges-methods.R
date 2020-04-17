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
#' @aliases granges
#' @seealso
#' \itemize{
#'   \item{\link[GenomicRanges]{granges} in the \pkg{GenomicRanges} package.}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod granges
setMethod(
  f = "granges",
  signature = "ChromatinAssay",
  definition = function(x, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for ChromatinAssay objects ",
           "does not support the 'use.mcols' argument")
    }
    slot(object = x, name = "ranges")
  }
)

#' @describeIn granges-methods method for Seurat objects
setMethod(
  f = "granges",
  signature = "Seurat",
  definition = function(x, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for Seurat objects ",
           "does not support the 'use.mcols' argument")
    }
    assay <- DefaultAssay(object = x)
    granges(x = x[[assay]])
  }
)

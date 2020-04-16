#' @importFrom GenomicRanges granges
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

#' @importFrom GenomicRanges granges
#' @importFrom Seurat DefaultAssay
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

#' @docType package
#' @name Signac-package
#' @rdname Signac-package
#' @concept assay
#'
"_PACKAGE"

.onLoad <- function(...) {
  vctrs::s3_register(generic = "Seurat::FoldChange", class = "ChromatinAssay")
}
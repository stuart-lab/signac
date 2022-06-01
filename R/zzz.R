#' @docType package
#' @name Signac-package
#' @rdname Signac-package
#'
"_PACKAGE"

.onLoad <- function(...) {
  vctrs::s3_register(generic = "Seurat::FoldChange", class = "ChromatinAssay")
}
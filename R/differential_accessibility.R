# Re-export Seurat generic
#' @importFrom Seurat FoldChange
#' @export
Seurat::FoldChange

#' @rdname FoldChange
#' @export
#' @method FoldChange ChromatinAssay
#' @importFrom Seurat FoldChange GetAssayData
#' @importFrom Matrix rowMeans
FoldChange.ChromatinAssay <- function(
  object,
  cells.1,
  cells.2,
  features = NULL,
  slot = "data",
  pseudocount.use = 1,
  fc.name = NULL,
  mean.fxn = NULL,
  base = 2,
  ...
) {
  mean.fxn <-  function(x) {
    return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
  }
  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- SetIfNull(x = fc.name, y = paste0("avg_log", base.text, "FC"))
  data <- GetAssayData(object = object, slot = slot)
  FoldChange(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )
}

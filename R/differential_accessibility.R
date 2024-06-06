# dynamically exported, see zzz.R
#' @method FoldChange ChromatinAssay
#' @importFrom SeuratObject GetAssayData
#' @importFrom Matrix rowMeans
FoldChange.ChromatinAssay <- function(
  object,
  cells.1,
  cells.2,
  features = NULL,
  slot = "data",
  fc.name = NULL,
  mean.fxn = NULL,
  base = 2,
  ...
) {
  if (!requireNamespace(package = "Seurat", quietly = TRUE)) {
    stop("Please install Seurat: install.packages('Seurat')")
  }
  if (is.null(x = mean.fxn)) {
    mean.fxn <-  function(x) {
      return(log(x = rowMeans(x = x) + 1/10000, base = base))
    }
  }
  # Omit the decimal value of e from the column name if base == exp(1)
  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- SetIfNull(x = fc.name, y = paste0("avg_log", base.text, "FC"))
  data <- GetAssayData(object = object, layer = slot)
  Seurat::FoldChange(
    object = data,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    mean.fxn = mean.fxn,
    fc.name = fc.name
  )
}

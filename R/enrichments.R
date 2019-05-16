#' FindMotifs
#'
#' Find motifs enriched in a given set of peaks
#'
#' @param object A Seurat object
#' @param cells A vector of cells to test
#' @param assay Which assay to use. Default is the active assay
#'
#' @importFrom Matrix rowSums
#' @return Returns a data.frame
#' @export
FindMotifs <- function(
  object,
  cells = NULL,
  features = NULL,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  cells <- cells %||% colnames(x = object)
  features <- features %||% rownames(x = object)
  feature.matrix <- GetAssayData(object = object, assay = assay, slot = 'data')[features, cells]
  motif.matrix <- GetMotifData(object = object, assay = assay, slot = 'data')[features, ]
  enrichment <- t(motif.matrix) %*% feature.matrix
  # TODO Check the Cusanovich papers to see how they did this.
  top.motifs <- sort(x = rowSums(x = enrichment), decreasing = TRUE)
  results <- data.frame(motif = names(x = top.motifs), score = top.motifs, row.names = names(top.motifs))
  return(results)
}

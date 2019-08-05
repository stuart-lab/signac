#' Integrate multiple datasets using Harmony
#'
#' Runs Harmony on a set of cell embeddings for the integration of multiple datasets.
#' Harmony was developed by Ilya Korsunsky et al., if you use Harmony please consider citing:
#' \url{https://doi.org/10.1101/461954}
#'
#' @param object A Seurat object
#' @param assay Name of the assay to use. If NULL, use the default assay.
#' @param group.by Grouping variable.
#' @param reduction Name of dimension reduction to harmonize
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[harmony]{HarmonyMatrix}}
#' @importFrom Seurat Embeddings CreateDimReducObject DefaultAssay
#' @return Returns a Seurat object
#' @export
Harmonize <- function(
  object,
  group.by,
  assay = NULL,
  reduction = 'lsi',
  verbose = TRUE,
  ...
) {
  if (!requireNamespace('harmony', quietly = TRUE)) {
    stop("Please install harmony. https://github.com/immunogenomics/harmony")
  }
  assay <- assay %||% DefaultAssay(object = object)
  embed <- Embeddings(
    object = object,
    reduction = reduction
  )
  metadata <- object[[]]
  harmonyEmbed <- harmony::HarmonyMatrix(
    data_mat = embed,
    meta_data = metadata,
    vars_use = group.by,
    do_pca = FALSE,
    npcs = 0,
    verbose = verbose,
    ...
  )
  rownames(harmonyEmbed) <- rownames(embed)
  colnames(harmonyEmbed) <- paste0("harmony_", seq_len(ncol(harmonyEmbed)))
  harmonydata <- CreateDimReducObject(
    embeddings = harmonyEmbed,
    assay = assay,
    key = "harmony"
  )
  object[['harmony']] <- harmonydata
  return(object)
}

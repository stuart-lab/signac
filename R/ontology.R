#' Find enriched ontology terms
#' 
#' Run ontology term enrichment testing on differential features for each of the
#' identity classes. Runs \code{\link[fgsea]{fgsea}} on differential testing
#' results from each identity class.
#' 
#' @param object A Seurat object.
#' @param terms Ontology term list. The name of each element in the list should
#' be an ontology term, and the list elements a vector of feature names present
#' in the assay.
#' @param group.by Name of grouping variable to use. If NULL, use the active cell
#' identities.
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param var.features Subset to only variable features for ontology term
#' enrichment.
#' @param verbose Display messages.
#' @param ... Additional arguments passed to \code{\link[Seurat]{FindMarkers}}
#' 
#' @importFrom SeuratObject VariableFeatures Idents RenameIdents DefaultAssay
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats setNames
#' 
#' @return Returns a named list of dataframes. Each element of the list contains
#' a dataframe with the term enrichment results for an identity class.
#' 
#' @export
EnrichedTerms <- function(
    object,
    terms,
    group.by = NULL,
    assay = NULL,
    var.features = TRUE,
    verbose = TRUE,
    ...
) {
  
  if (!requireNamespace(package = "Seurat", quietly = TRUE)) {
    stop("Please install Seurat: install.packages('Seurat')")
  }
  if (!requireNamespace(package = "fgsea", quietly = TRUE)) {
    stop("Please install fgsea: BiocManager::install('fgsea')")
  }
  
  assay <- assay %||% DefaultAssay(object = object)
  
  if (var.features) {
    terms <- lapply(
      X = terms,
      FUN = function(x) {
        x[x %in% VariableFeatures(object = object[[assay]])]
      }
    )
  }
  
  # for each identity class in object, find the top terms
  if (is.null(x = group.by)) {
    cellgroups <- unique(x = Idents(object = object))
  }
  else {
    meta.data <- object[[]]
    cellgroups <- unique(x = meta.data[[group.by]])
  }
  pred <- list()
  if (verbose) pb <- txtProgressBar(
    min = 1,
    max = length(x = cellgroups),
    style = 3
  )
  for (i in seq_along(along.with = cellgroups)) {
    mk <- Seurat::FindMarkers(
      object = object,
      assay = assay,
      ident.1 = cellgroups[[i]],
      group.by = group.by,
      ...
    )
    mk$rank_score <- sign(mk$avg_log2FC) * -log10(mk$p_val)
    ranked_list <- setNames(object = mk$rank_score, nm = rownames(x = mk))
    max_non_infinite <- max(ranked_list[!is.infinite(x = ranked_list)])
    ranked_list[is.infinite(x = ranked_list)] <- max_non_infinite
    
    # run fgsea
    fgsea_results <- fgsea::fgsea(
      pathways = terms,
      stats = ranked_list,
      scoreType = "pos"
    )
    fgsea_results <- fgsea_results[fgsea_results$NES > 0, ]
    fgsea_results <- fgsea_results[fgsea_results$padj < 0.05, ]
    fgsea_results <- fgsea_results[order(fgsea_results$NES, fgsea_results$padj, decreasing = c(TRUE, FALSE)), ]
    pred[[as.character(cellgroups[i])]] <- fgsea_results
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  return(pred)
}

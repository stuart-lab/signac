#' Find enriched ontology terms
#'
#' Run ontology term enrichment testing on differential features for each of the
#' identity classes. Runs [fgsea::fgsea()] on differential testing
#' results from each identity class.
#'
#' @param object A Seurat object.
#' @param terms Ontology term list. The name of each element in the list should
#' be an ontology term, and the list elements a vector of feature names present
#' in the assay.
#' @param group.by Name of grouping variable to use. If `NULL`, use the active
#' cell identities.
#' @param assay Name of assay to use. If `NULL`, use the default assay.
#' @param var.features Restrict the analysis to variable features. When `TRUE`
#' (the default), differential testing is performed only on the variable
#' features of the assay, so the ranked list scored by [fgsea::fgsea()] (the
#' enrichment universe) and the term gene sets are evaluated within the same
#' variable-feature space. Requires variable features to be set for the assay
#' (e.g. with [FindTopFeatures()]).
#' @param scoreType `scoreType` parameter for [fgsea::fgseaSimple()].
#' Options are "std", "pos", "neg" (two-tailed or one-tailed tests).
#' @param direction Which direction of enrichment to retain. `"up"` (default)
#' keeps terms with positive NES (enriched in the identity class), `"down"`
#' keeps terms with negative NES (depleted), and `"both"` retains both and
#' orders terms by `abs(NES)`. When `scoreType` is `"pos"` or `"neg"` only
#' one direction is testable, so set `direction` accordingly.
#' @param top.n Number of top enriched terms to retain for each set of cells. If
#' NULL, retain all terms.
#' @param verbose Display messages.
#' @param ... Additional arguments passed to [Seurat::FindMarkers()]
#'
#' @importFrom SeuratObject VariableFeatures Idents RenameIdents DefaultAssay
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats setNames
#'
#' @return Returns a named list of dataframes. Each element of the list contains
#' a dataframe with the term enrichment results for an identity class.
#' @concept utilities
#' @export
EnrichedTerms <- function(
  object,
  terms,
  group.by = NULL,
  assay = NULL,
  var.features = TRUE,
  scoreType = "std",
  direction = c("up", "down", "both"),
  top.n = NULL,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace(package = "Seurat", quietly = TRUE)) {
    stop("Please install Seurat: install.packages('Seurat')")
  }
  if (!requireNamespace(package = "fgsea", quietly = TRUE)) {
    stop("Please install fgsea: BiocManager::install('fgsea')")
  }
  direction <- match.arg(arg = direction)

  assay <- assay %||% DefaultAssay(object = object)

  marker.args <- list(...)
  if (var.features) {
    var.feat <- VariableFeatures(object = object[[assay]])
    if (length(x = var.feat) == 0) {
      stop(
        "var.features = TRUE but no variable features are set for assay '",
        assay, "'. Compute variable features first (for example with ",
        "FindTopFeatures()) or set var.features = FALSE."
      )
    }
    # Restrict the analysis to variable features
    marker.args$features <- if (is.null(x = marker.args$features)) {
      var.feat
    } else {
      intersect(x = marker.args$features, y = var.feat)
    }
    terms <- lapply(
      X = terms,
      FUN = function(x) x[x %in% var.feat]
    )
  }

  # for each identity class in object, find the top terms
  if (is.null(x = group.by)) {
    cellgroups <- unique(x = Idents(object = object))
  } else {
    meta.data <- object[[]]
    cellgroups <- unique(x = meta.data[[group.by]])
  }
  pred <- list()
  if (verbose) {
    pb <- txtProgressBar(
      min = 1,
      max = length(x = cellgroups),
      style = 3
    )
  }
  for (i in seq_along(along.with = cellgroups)) {
    mk <- do.call(
      what = Seurat::FindMarkers,
      args = c(
        list(
          object = object,
          assay = assay,
          ident.1 = cellgroups[[i]],
          group.by = group.by
        ),
        marker.args
      )
    )
    ranked_list <- RankFeatures(markers = mk)

    # run fgsea
    fgsea_results <- fgsea::fgsea(
      pathways = terms,
      stats = ranked_list,
      scoreType = scoreType
    )
    fgsea_results <- switch(
      EXPR = direction,
      up = fgsea_results[fgsea_results$NES > 0, ],
      down = fgsea_results[fgsea_results$NES < 0, ],
      both = fgsea_results
    )
    fgsea_results <- fgsea_results[fgsea_results$padj < 0.05, ]
    sort_score <- switch(
      EXPR = direction,
      up = fgsea_results$NES,
      down = -fgsea_results$NES,
      both = abs(x = fgsea_results$NES)
    )
    fgsea_results <- fgsea_results[order(
      sort_score,
      fgsea_results$padj,
      decreasing = c(TRUE, FALSE)
    ), ]
    if (!is.null(x = top.n)) {
      n.use <- min(nrow(x = fgsea_results), top.n)
      fgsea_results <- fgsea_results[1:n.use, ]
    }
    pred[[as.character(cellgroups[i])]] <- fgsea_results
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }
  return(pred)
}

# Build a named ranked statistic vector for fgsea from a FindMarkers result.
# The score for each feature is sign(avg_log2FC) * -log10(p_val), so that
# strongly up-regulated features get large positive scores and strongly
# down-regulated features get large negative scores. A p-value of exactly zero
# produces an infinite score (+Inf when up-regulated, -Inf when down-regulated);
# each infinite score is capped at the most extreme finite score of the same
# sign. Capping the two signs separately keeps up-regulated features at the top
# of the ranking and down-regulated features at the bottom
# @param markers A data frame of differential test results with columns
#   `avg_log2FC` and `p_val` and feature names as row names (as returned by
#   [Seurat::FindMarkers()]).
# @return A named numeric vector of ranking statistics.
RankFeatures <- function(markers) {
  rank_score <- sign(x = markers$avg_log2FC) * -log10(x = markers$p_val)
  ranked_list <- setNames(object = rank_score, nm = rownames(x = markers))
  finite_scores <- ranked_list[is.finite(x = ranked_list)]
  pos_cap <- if (length(x = finite_scores) > 0) max(finite_scores) else 1
  neg_cap <- if (length(x = finite_scores) > 0) min(finite_scores) else -1
  ranked_list[ranked_list == Inf] <- pos_cap
  ranked_list[ranked_list == -Inf] <- neg_cap
  return(ranked_list)
}

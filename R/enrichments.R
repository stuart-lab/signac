#' MotifCellEnrichment
#'
#' Find motifs enriched in a given set of peaks for a set of cells.
#'
#' Performs matrix multiplication between a motif x feature and feature x cell matrix to produce motif
#' counts per cell. This is repeated using a set of background features to calculate the relative enrichment
#' of each motif in tested features. A permutation test can also be performed by setting the \code{permute}
#' argument, and the enrichment testing will be repeated \emph{n} times using features sampled at random from
#' the background set. A p-value is then computed for each motif, which is equal to the number of times
#' a permuted enrichment score greater or equal to the tested enrichment score was observed.
#'
#' @param object A Seurat object
#' @param features A vector of features to test
#' @param cells A vector of cells to include in enrichment tests
#' @param assay Which assay to use. Default is the active assay
#' @param permute Number of times to permute the set of features. If NULL, do not perform any permutation
#' @param background Set of background features to use when calculating relative enrichment. If NULL
#' use all features present in the object
#' @param verbose Display messages
#'
#' @importFrom Matrix rowSums
#' @importFrom future plan
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#'
#' @return Returns a data.frame
#' @export
MotifCellEnrichment <- function(
  object,
  cells = NULL,
  features = NULL,
  assay = NULL,
  permute = NULL,
  background = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  cells <- cells %||% colnames(x = object)
  features <- features %||% rownames(x = object)
  background <- background %||% rownames(x = object)
  if (verbose) {
    message('Testing motif enrichment in ', length(features), ' regions')
  }
  data.use <- GetAssayData(object = object, assay = assay, slot = 'counts')[, cells]
  motif.all <- GetMotifData(object = object, assay = assay, slot = 'data')
  pwm <- GetMotifData(object = object, assay = assay, slot = 'pwm')
  if (class(x = pwm) == 'PFMatrixList') {
    motif.names <- name(x = pwm)
  } else {
    motif.names <- NULL
  }
  top.motifs <- TestEnrichment(
    motif.matrix = motif.all,
    feature.matrix = data.use[features, ]
  )
  background.scores <- TestEnrichment(
    motif.matrix = motif.all,
    feature.matrix = data.use[background, ]
  )
  enrichment <- top.motifs / background.scores
  results <- data.frame(
    motif = names(x = enrichment),
    score = enrichment,
    row.names = names(enrichment)
  )
  if (!is.null(x = motif.names)) {
    results$motif.name <- motif.names
  }
  if (is.null(x = permute)) {
    return(results[order(results$enrichment, decreasing = TRUE), ])
  } else {
    message("Permuting feature sets ", permute, " times")
    n.features <- length(features)
    if (PlanThreads() > 1) {
      mysapply <- future_sapply
    } else {
      mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
    }
    permuted.scores <- mysapply(
      X = 1:permute,
      FUN = function(x) {
        rand.sample <- sample(x = background, size = n.features, replace = FALSE)
        permuted.enrichment <- TestEnrichment(
          motif.matrix = motif.all,
          feature.matrix = data.use[rand.sample, ]
        )
        return(permuted.enrichment / background.scores)
      })
    test.results <- permuted.scores >= enrichment
    p.vals <- rowSums(test.results) / permute
    results$pvalue <- p.vals
    return(results[with(data = results, expr = order(pvalue, -score)), ])
  }
}

#' FindMotifs
#'
#' Find motifs overrepresented in a given set of genomic features. Computes the number of features
#' containing the motif (observed) and compares this to the total number of features containing the
#' motif (background) using the hypergeometric test.
#'
#' @param object A Seurat object
#' @param features A vector of features to test for enrichments over background
#' @param assay Which assay to use. Default is the active assay
#' @param background Vector of features to use as the background set. If NULL, use all features in the assay.
#' @param verbose Display messages
#'
#' @return Returns a data frame
#'
#' @importFrom Matrix colSums
#' @importFrom stats phyper
#' @importFrom future plan
#'
#' @export
FindMotifs <- function(
  object,
  features,
  assay = NULL,
  background = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  background <- background %||% rownames(x = object)
  if (verbose) {
    message('Testing motif enrichment in ', length(features), ' regions')
  }
  motif.all <- GetMotifData(object = object, assay = assay, slot = 'data')
  pwm <- GetMotifData(object = object, assay = assay, slot = 'pwm')
  if (class(x = pwm) == 'PFMatrixList') {
    motif.names <- name(x = pwm)
  } else {
    motif.names <- NULL
  }
  subs.motifs <- motif.all[features, ]
  subs.counts <- colSums(x = subs.motifs)
  all.counts <- colSums(x = motif.all)
  obs.expect <- subs.counts / all.counts
  p.list <- c()
  for (i in seq_along(along.with = subs.counts)) {
    p.list[[i]] <- phyper(
      q = subs.counts[[i]]-1,
      m = all.counts[[i]],
      n = nrow(x = motif.all) - all.counts[[i]],
      k = length(x = features),
      lower.tail = FALSE
    )
  }
  results <- data.frame(
    motif = names(x = subs.counts),
    observed = subs.counts,
    background = all.counts,
    enrichment = obs.expect,
    pvalue = p.list
  )
  if (!is.null(x = motif.names)) {
    results$motif.name <- motif.names
  }
  return(results[with(data = results, expr = order(pvalue, -enrichment)), ])
}

#' TestEnrichment
#'
#' Compute motif counts per cell for a given set of features
#'
#' @param motif.matrix A feature x motif matrix
#' @param feature.matrix A feature x cell matrix
#' @return Returns a vector of motif counts
TestEnrichment <- function(motif.matrix, feature.matrix) {
  features = rownames(x = feature.matrix)
  motif.matrix <- motif.matrix[features, ]
  enrichment <- crossprod(x = motif.matrix, y = feature.matrix)
  motif.counts <- rowSums(x = enrichment)
  return(motif.counts)
}

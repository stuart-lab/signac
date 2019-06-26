#' @include generics.R
#'
NULL

#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[Seurat]{FindNeighbors}}
#' and \code{\link[Seurat]{FindClusters}}
#' @rdname ClusterMotifs
#' @method ClusterMotifs Motif
#' @importFrom Matrix crossprod colSums
#' @importFrom Seurat FindNeighbors FindClusters
#' @export
ClusterMotifs.Motif <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  data.use <- t(x = GetMotifData(object = object, slot = 'data'))
  motif.jaccard <- Jaccard(x = data.use, y = data.use)
  object <- SetMotifData(
    object = object,
    slot = 'neighbors',
    new.data = FindNeighbors(
      object = 1/motif.jaccard,
      distance.matrix = TRUE,
      verbose = verbose,
      ...
    )
  )
  clusters <- FindClusters(
    object = GetMotifData(object = object, slot = 'neighbors')$nn,
    verbose = verbose,
    ...
  )
  meta.data <- GetMotifData(object = object, slot = 'meta.data')
  if (nrow(x = meta.data) == 0) {
    meta.data <- clusters
  } else {
    meta.data[[colnames(x = clusters)]] <- clusters[, 1]
  }
  object <- SetMotifData(object = object, slot = 'meta.data', new.data = meta.data)
  return(object)
}

#' @rdname ClusterMotifs
#' @method ClusterMotifs Assay
#' @export
ClusterMotifs.Assay <- function(
  object,
  verbose = TRUE,
  ...
) {
  motif.obj <- GetMotifObject(object = object)
  motif.obj <- ClusterMotifs(object = motif.obj, verbose = verbose, ...)
  object <- AddMotifObject(object = object, motif.object = motif.obj, verbose = FALSE)
  return(object)
}

#' @param assay Which assay to use
#' @rdname ClusterMotifs
#' @method ClusterMotifs Seurat
#' @export
ClusterMotifs.Seurat <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- ClusterMotifs(
    object = assay.data,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' CreateMotifActivityMatrix
#'
#' Create a matrix of the normalized motif accessibility per cell.
#' Fits a loess model per motif between the total cell accessibility and the sequencing depth for the cell,
#' and stores the model residual as the normalized motif activity in the cell.
#'
#' @param object A Seurat object
#' @param assay Which assay to use. This must contain a Motif object.
#' @param span Degree of loess smoothing, passed to \code{\link[stats]{loess}}. Default is 0.5.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[stats]{loess}}
#'
#' @importFrom Matrix crossprod colSums
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_sapply
#' @importFrom pbapply pbsapply
#' @importFrom stats loess
#'
#' @return Returns a matrix
#' @export
CreateMotifActivityMatrix <- function(
  object,
  assay = NULL,
  span = 0.5,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  motifs <- GetMotifData(object = object, assay = assay, slot = 'data')
  accessibility <- GetAssayData(object = object, assay = assay, slot = 'counts')
  if (verbose) {
    message("Computing motif accessibility per cell")
  }
  motif.accessibility <- as.matrix(x = crossprod(x = motifs, y = accessibility))
  seq.depth <- colSums(x = accessibility)
  if (nbrOfWorkers() > 1) {
    mysapply <- future_sapply
  } else {
    mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  }
  if (verbose) {
    message("Fitting models")
  }
  residuals <- mysapply(
    X = seq_along(along.with = rownames(x = motif.accessibility)),
    FUN = function(x) {
      loess(formula = motif.accessibility[x, ] ~ seq.depth, span = span, ...)$residuals
    }
  )
  colnames(residuals) <- rownames(x = motif.accessibility)
  return(t(residuals))
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
    message('Testing motif enrichment in ', length(x = features), ' regions')
  }
  motif.all <- GetMotifData(object = object, assay = assay, slot = 'data')
  pwm <- GetMotifData(object = object, assay = assay, slot = 'pwm')
  if (is(object = pwm, class2 = 'PFMatrixList')) {
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
    pvalue = p.list,
    stringsAsFactors = FALSE
  )
  if (!is.null(x = motif.names)) {
    results$motif.name <- motif.names
  }
  return(results[with(data = results, expr = order(pvalue, -enrichment)), ])
}

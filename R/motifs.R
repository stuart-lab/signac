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

#' Run chromVAR
#'
#' Wrapper to run \code{\link[chromVAR]{chromVAR}} on an assay with a motif object present.
#' Will return a new Seurat assay with the motif activities stored.
#'
#' @param object A Seurat object
#' @param genome A BSgenome object
#' @param assay Name of assay to use
#' @param new.assay.name Name of new assay used to store the chromVAR results. Default is "chromvar".
#' @param motif.matrix A peak x motif matrix. If NULL, pull the peak x motif matrix from a Motif object stored in the assay.
#' @param sep A length-2 character vector containing the separators passed to \code{\link{StringToGRanges}}.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[chromVAR]{getBackgroundPeaks}}
#'
#' @importFrom Seurat GetAssayData DefaultAssay CreateAssayObject
#' @importFrom Matrix rowSums
#'
#' @return Returns a \code{\link[Seurat]{Seurat}} object with a new assay
#'
#' @export
#'
RunChromVAR <- function(
  object,
  genome,
  new.assay.name = 'chromvar',
  motif.matrix = NULL,
  assay = NULL,
  sep = c(":", "-"),
  verbose = TRUE,
  ...
) {
  if (!requireNamespace('chromVAR', quietly = TRUE)) {
    stop("Please install chromVAR. https://greenleaflab.github.io/chromVAR/")
  }
  if (!requireNamespace('SummarizedExperiment', quietly = TRUE)) {
    stop("Please install SummarizedExperiment")
  }
  assay <- assay %||% DefaultAssay(object = object)
  motif.matrix <- motif.matrix %||% GetMotifData(object = object, assay = assay, slot = 'data')
  peak.matrix <- GetAssayData(object = object, assay = assay, slot = 'counts')
  peak.matrix <- peak.matrix[rowSums(x = peak.matrix) > 0, ]
  motif.matrix <- motif.matrix[rownames(x = peak.matrix), ]
  peak.ranges <- StringToGRanges(regions = rownames(peak.matrix), sep = sep)

  chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak.matrix),
    rowRanges = peak.ranges
  )
  if (verbose) {
    message("Computing GC bias per region")
  }
  chromvar.obj <- chromVAR::addGCBias(
    object = chromvar.obj,
    genome = genome
  )
  if (verbose) {
    message("Selecting background regions")
  }
  bg <- chromVAR::getBackgroundPeaks(
    object = chromvar.obj,
    ...
  )
  if (verbose) {
    message("Computing motif deviations from background")
  }
  dev <- chromVAR::computeDeviations(
    object = chromvar.obj,
    annotations = motif.matrix,
    background_peaks = bg
  )
  chromvar.z <- SummarizedExperiment::assays(dev)[[2]]
  if (verbose) {
    message("Constructing chromVAR assay")
  }
  object[['chromvar']] <- CreateAssayObject(data = chromvar.z)
  return(object)
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
#' @param background Either a vector of features to use as the background set,
#' or a number specify the number of features to randomly select as a background set.
#' If a number is provided, regions will be selected to match the sequence characteristics
#' of the query features. To match the sequence characteristics, these characteristics
#' must be stored in the feature metadata for the assay. This can be added using the
#'  \code{\link{RegionStats}} function. If NULL, use all features in the assay.
#' @param verbose Display messages
#' @param ... Arguments passed to \code{\link{MatchRegionStats}}.
#'
#' @return Returns a data frame
#'
#' @importFrom Matrix colSums
#' @importFrom stats phyper
#' @importFrom methods is
#'
#' @export
FindMotifs <- function(
  object,
  features,
  background = 40000,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  background <- background %||% rownames(x = object)
  if (is(object = background, class2 = 'numeric')) {
    if (verbose) {
      message("Selecting background regions to match input sequence characteristics")
    }
    background <- MatchRegionStats(
      meta.feature = GetAssayData(object = object, assay = assay, slot = 'meta.features'),
      regions = features,
      n = background,
      verbose = verbose,
      ...
    )
  }
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
  query.motifs <- motif.all[features, ]
  background.motifs <- motif.all[background, ]
  query.counts <- colSums(x = query.motifs)
  background.counts <- colSums(x = background.motifs)
  percent.observed <- query.counts / length(x = features) * 100
  percent.background <- background.counts / length(x = background) * 100
  fold.enrichment <- percent.observed / percent.background
  p.list <- c()
  for (i in seq_along(along.with = query.counts)) {
    p.list[[i]] <- phyper(
      q = query.counts[[i]]-1,
      m = background.counts[[i]],
      n = nrow(x = background.motifs) - background.counts[[i]],
      k = length(x = features),
      lower.tail = FALSE
    )
  }
  results <- data.frame(
    motif = names(x = query.counts),
    observed = query.counts,
    background = background.counts,
    percent.observed = percent.observed,
    percent.background = percent.background,
    fold.enrichment = fold.enrichment,
    pvalue = p.list,
    stringsAsFactors = FALSE
  )
  if (!is.null(x = motif.names)) {
    results$motif.name <- motif.names
  }
  return(results[with(data = results, expr = order(pvalue, -fold.enrichment)), ])
}

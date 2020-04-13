#' @include generics.R
#'
NULL

#' Run chromVAR
#'
#' Wrapper to run \code{\link[chromVAR]{chromVAR}} on an assay with a motif
#' object present. Will return a new Seurat assay with the motif activities
#' (the deviations in chromatin accessibility across the set of regions) as
#' a new assay.
#'
#' See the chromVAR documentation for more information:
#' \url{https://greenleaflab.github.io/chromVAR/index.html}
#'
#' See the chromVAR paper: \url{https://www.nature.com/articles/nmeth.4401}
#'
#' @param object A Seurat object
#' @param genome A BSgenome object
#' @param assay Name of assay to use
#' @param new.assay.name Name of new assay used to store the chromVAR results.
#' Default is "chromvar".
#' @param motif.matrix A peak x motif matrix. If NULL, pull the peak x motif
#' matrix from a Motif object stored in the assay.
#' @param sep A length-2 character vector containing the separators passed to
#' \code{\link{StringToGRanges}}.
#' @param verbose Display messages
#' @param ... Additional arguments passed to
#' \code{\link[chromVAR]{getBackgroundPeaks}}
#'
#' @importFrom Seurat GetAssayData DefaultAssay CreateAssayObject
#' @importFrom Matrix rowSums
#'
#' @return Returns a \code{\link[Seurat]{Seurat}} object with a new assay
#'
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RunChromVAR(object = atac_small, genome = BSgenome.Hsapiens.UCSC.hg19)
#' }
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
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  motif.matrix <- SetIfNull(
    x = motif.matrix,
    y = GetMotifData(object = object, assay = assay, slot = 'data')
  )
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
  rownames(x = chromvar.z) <- colnames(x = motif.matrix)
  if (verbose) {
    message("Constructing chromVAR assay")
  }
  object[['chromvar']] <- CreateAssayObject(data = chromvar.z)
  return(object)
}

globalVariables(names = 'pvalue', package = 'Signac')
#' FindMotifs
#'
#' Find motifs overrepresented in a given set of genomic features.
#' Computes the number of features containing the motif (observed) and
#' compares this to the total number of features containing the
#' motif (background) using the hypergeometric test.
#'
#' @param object A Seurat object
#' @param features A vector of features to test for enrichments over background
#' @param assay Which assay to use. Default is the active assay
#' @param background Either a vector of features to use as the background set,
#' or a number specify the number of features to randomly select as a background
#' set. If a number is provided, regions will be selected to match the sequence
#' characteristics of the query features. To match the sequence characteristics,
#' these characteristics must be stored in the feature metadata for the assay.
#' This can be added using the
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
#' @examples
#' de.motif <- head(rownames(atac_small))
#' bg.peaks <- tail(rownames(atac_small))
#' FindMotifs(
#' object = atac_small,
#' features = de.motif,
#' background = bg.peaks
#' )
FindMotifs <- function(
  object,
  features,
  background = 40000,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  background <- SetIfNull(x = background, y = rownames(x = object))
  if (is(object = background, class2 = 'numeric')) {
    if (verbose) {
      message("Selecting background regions to match input
              sequence characteristics")
    }
    background <- MatchRegionStats(
      meta.feature = GetAssayData(
        object = object,
        assay = assay,
        slot = 'meta.features'
      ),
      regions = features,
      n = background,
      verbose = verbose,
      ...
    )
  }
  if (verbose) {
    message('Testing motif enrichment in ', length(x = features), ' regions')
  }
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = 'data'
  )
  motif.names <- GetMotifData(
    object = object, assay = assay, slot = 'motif.names'
  )
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
      q = query.counts[[i]] - 1,
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
    motif.name = as.vector(
      x = unlist(x = motif.names[names(x = query.counts)])
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(x = results) == 0) {
    return(results)
  } else {
    return(results[order(-results[, 7], -results[, 6]), ])
  }
}

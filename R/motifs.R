#' @include generics.R
#'
NULL

#' @importFrom Seurat GetAssayData CreateAssayObject
#' @importFrom Matrix rowSums
#'
#' @concept motifs
#' @method RunChromVAR ChromatinAssay
#' @rdname RunChromVAR
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RunChromVAR(object = atac_small[["peaks"]], genome = BSgenome.Hsapiens.UCSC.hg19)
#' }
RunChromVAR.ChromatinAssay <- function(
  object,
  genome,
  motif.matrix = NULL,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("chromVAR", quietly = TRUE)) {
    stop("Please install chromVAR. https://greenleaflab.github.io/chromVAR/")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment")
  }
  motif.matrix <- SetIfNull(
    x = motif.matrix,
    y = GetMotifData(object = object, slot = "data")
  )
  peak.matrix <- GetAssayData(object = object, slot = "counts")
  if (!(all(peak.matrix@x == floor(peak.matrix@x)))) {
    warning("Count matrix contains non-integer values.
            ChromVAR should only be run on integer counts.")
  }
  idx.keep <- rowSums(x = peak.matrix) > 0
  peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
  motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
  peak.ranges <- granges(x = object)
  peak.ranges <- peak.ranges[idx.keep]
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
    message("Computing deviations from background")
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
  obj <- CreateAssayObject(data = chromvar.z)
  return(obj)
}

#' @param assay Name of assay to use
#' @param new.assay.name Name of new assay used to store the chromVAR results.
#' Default is "chromvar".
#' @method RunChromVAR Seurat
#' @rdname RunChromVAR
#' @export
#' @importFrom Seurat DefaultAssay
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RunChromVAR(object = atac_small, genome = BSgenome.Hsapiens.UCSC.hg19)
#' }
RunChromVAR.Seurat <- function(
  object,
  genome,
  motif.matrix = NULL,
  assay = NULL,
  new.assay.name = "chromvar",
  ...
) {
  if (!requireNamespace("chromVAR", quietly = TRUE)) {
    stop("Please install chromVAR. https://greenleaflab.github.io/chromVAR/")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment")
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  chromvar.assay <- RunChromVAR(
    object = object[[assay]],
    genome = genome,
    motif.matrix = motif.matrix,
    ...
  )
  object[[new.assay.name]] <- chromvar.assay
  return(object)
}

globalVariables(names = "pvalue", package = "Signac")
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
#' @importFrom qvalue qvalue
#'
#' @export
#' @concept motifs
#' @examples
#' de.motif <- head(rownames(atac_small))
#' bg.peaks <- tail(rownames(atac_small))
#' FindMotifs(
#'   object = atac_small,
#'   features = de.motif,
#'   background = bg.peaks
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
  if (is(object = background, class2 = "numeric")) {
    if (verbose) {
      message("Selecting background regions to match input
              sequence characteristics")
    }
    background <- MatchRegionStats(
      meta.feature = GetAssayData(
        object = object,
        assay = assay,
        slot = "meta.features"
      ),
      regions = features,
      n = background,
      verbose = verbose,
      ...
    )
  }
  if (verbose) {
    message("Testing motif enrichment in ", length(x = features), " regions")
  }
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  motif.names <- GetMotifData(
    object = object, assay = assay, slot = "motif.names"
  )
  query.motifs <- motif.all[features, ]
  background.motifs <- motif.all[background, ]
  query.counts <- colSums(x = query.motifs)
  background.counts <- colSums(x = background.motifs)
  percent.observed <- query.counts / length(x = features) * 100
  percent.background <- background.counts / length(x = background) * 100
  fold.enrichment <- percent.observed / percent.background
  p.list <- vector(mode = "numeric")
  for (i in seq_along(along.with = query.counts)) {
    p.list[[i]] <- phyper(
      q = query.counts[[i]] - 1,
      m = background.counts[[i]],
      n = nrow(x = background.motifs) - background.counts[[i]],
      k = length(x = features),
      lower.tail = FALSE
    )
  }
  qv <- qvalue(p = p.list)
  results <- data.frame(
    motif = names(x = query.counts),
    observed = query.counts,
    background = background.counts,
    percent.observed = percent.observed,
    percent.background = percent.background,
    fold.enrichment = fold.enrichment,
    pvalue = p.list,
    qvalue = qv$qvalues,
    motif.name = as.vector(
      x = unlist(x = motif.names[names(x = query.counts)])
    ),
    stringsAsFactors = FALSE
  )
  if (nrow(x = results) == 0) {
    return(results)
  } else {
    return(results[order(results[, 7], -results[, 6]), ])
  }
}

#' @param name A vector of motif names
#' @param id A vector of motif IDs. Only one of \code{name} and \code{id} should
#' be supplied
#' @rdname ConvertMotifID
#' @concept motifs
#' @importFrom methods hasArg
#' @export
ConvertMotifID.default <- function(object, name, id, ...) {
  if (hasArg(name = name) & hasArg(name = id)) {
    stop("Supply either name or ID, not both")
  } else if (!hasArg(name = name) & !(hasArg(name = id))) {
    stop("Supply vector of names or IDs to convert")
  } else {
    if (hasArg(name = name)) {
      # convert name to ID
      # construct a new vector for conversion
      name.to.id <- names(x = object)
      names(x = name.to.id) <- object
      converted.names <- as.vector(x = name.to.id[name])
    } else {
      # convert ID to name
      tmp <- object[id]
      # for missing motif, change from NULL to NA
      tmp[is.na(x = names(x = tmp))] <- NA
      converted.names <- unlist(x = tmp, use.names = FALSE)
    }
    return(converted.names)
  }
}

#' @method ConvertMotifID Motif
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.Motif <- function(object, ...) {
  motif.names <- GetMotifData(object = object, slot = "motif.names")
  return(ConvertMotifID(object = motif.names, ...))
}

#' @method ConvertMotifID ChromatinAssay
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.ChromatinAssay <- function(object, ...) {
  motifs <- Motifs(object = object)
  return(ConvertMotifID(object = motifs, ...))
}

#' @param assay For \code{Seurat} objectd. Name of assay to use.
#' If NULL, use the default assay
#'
#' @importFrom Seurat DefaultAssay
#'
#' @method ConvertMotifID Seurat
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.Seurat <- function(object, assay = NULL, ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  return(ConvertMotifID(object = object[[assay]], ...))
}

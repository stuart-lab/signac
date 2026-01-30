#' @include generics.R
#'
NULL

#' @rdname AddMotifs
#' @method AddMotifs default
#' @concept motifs
#' @importFrom methods slot
#' @importFrom Seqinfo seqlevels seqnames
#' @export
AddMotifs.default <- function(
  object,
  genome,
  pfm,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop(
      "Please install motifmatchr.\n",
      "https://www.bioconductor.org/packages/motifmatchr/"
    )
  }
  if (is.null(x = names(x = pfm))) {
    warning("No 'names' attribute found in PFMatrixList. ",
      "Extracting names from individual entries.",
      immediate. = TRUE
    )
    names(x = pfm) <- vapply(
      X = pfm, FUN = slot, FUN.VALUE = "character", "name"
    )
  }
  if (verbose) {
    message("Building motif matrix")
  }
  # genome can be string
  if (is.character(x = genome)) {
    if (!requireNamespace("BSgenome", quietly = TRUE)) {
      stop("Please install BSgenome.
             https://www.bioconductor.org/packages/BSgenome/")
    }
    genome <- BSgenome::getBSgenome(genome = genome)
  }
  motif.matrix <- CreateMotifMatrix(
    features = object,
    pwm = pfm,
    genome = genome,
    use.counts = FALSE
  )
  if (verbose) {
    message("Finding motif positions")
  }

  # for positions, a list of granges is returned
  # each element of list is a PFM name
  # each entry in granges is the position within a feature that matches motif
  obj_keep <- as.character(seqnames(x = object)) %in% seqlevels(x = genome)
  motif.positions <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = object[obj_keep],
    out = "positions",
    genome = genome
  )
  if (verbose) {
    message("Creating Motif object")
  }
  motif <- CreateMotifObject(
    data = motif.matrix,
    positions = motif.positions,
    pwm = pfm
  )
  return(motif)
}

#' @rdname AddMotifs
#' @method AddMotifs GRangesAssay
#' @concept motifs
#' @export
AddMotifs.GRangesAssay <- function(
  object,
  genome,
  pfm,
  verbose = TRUE,
  ...
) {
  motif <- AddMotifs(
    object = granges(x = object),
    genome = genome,
    pfm = pfm,
    verbose = verbose
  )
  object <- SetAssayData(
    object = object,
    layer = "motifs",
    new.data = motif
  )
  return(object)
}

#' @rdname AddMotifs
#' @method AddMotifs Assay
#' @concept motifs
#' @export
AddMotifs.Assay <- function(
  object,
  genome,
  pfm,
  verbose = TRUE,
  ...
) {
  stop(
    "Attempting to run AddMotifs on a standard Assay.\n",
    "Please supply a GRangesAssay instead."
  )
}

#' @rdname AddMotifs
#' @method AddMotifs StdAssay
#' @concept motifs
#' @export
AddMotifs.StdAssay <- function(
  object,
  genome,
  pfm,
  verbose = TRUE,
  ...
) {
  stop(
    "Attempting to run AddMotifs on an Assay5 assay.\n",
    "Please supply a GRangesAssay instead."
  )
}

#' @param assay Name of assay to use. If NULL, use the default assay
#' @param genome A `BSgenome`, `DNAStringSet`, `FaFile`, or
#' string stating the genome build recognized by `getBSgenome`.
#' @param pfm A `PFMatrixList` or `PWMatrixList` object containing
#' position weight/frequency matrices to use
#' @param verbose Display messages
#' @importFrom SeuratObject DefaultAssay
#' @rdname AddMotifs
#' @method AddMotifs Seurat
#' @concept motifs
#' @export
AddMotifs.Seurat <- function(
  object,
  genome,
  pfm,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- AddMotifs(
    object = object[[assay]],
    genome = genome,
    pfm = pfm,
    verbose = verbose
  )
  object <- RegionStats(
    object = object,
    assay = assay,
    genome = genome,
    verbose = verbose
  )
  return(object)
}

#' Create motif matrix
#'
#' Create a motif x feature matrix from a set of genomic ranges,
#' the genome, and a set of position weight matrices.
#'
#' Requires that motifmatchr is installed
#' <https://www.bioconductor.org/packages/motifmatchr/>.
#'
#' @param features A GRanges object containing a set of genomic features
#' @param pwm A [TFBSTools::PFMatrixList()] or
#' [TFBSTools::PWMatrixList()]
#' object containing position weight/frequency matrices to use
#' @param genome Any object compatible with the `genome` argument
#' in [motifmatchr::matchMotifs()]
#' @param score Record the motif match score, rather than presence/absence
#' (default FALSE)
#' @param use.counts Record motif counts per region. If FALSE (default),
#' record presence/absence of motif. Only applicable if `score=FALSE`.
#' @param ... Additional arguments passed to
#' [motifmatchr::matchMotifs()]
#'
#' @return Returns a sparse matrix
#' @export
#' @concept motifs
#' @concept preprocessing
#' @examples
#' \dontrun{
#' library(JASPAR2018)
#' library(TFBSTools)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' pwm <- getMatrixSet(
#'   x = JASPAR2018,
#'   opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
#' )
#' motif.matrix <- CreateMotifMatrix(
#'   features = granges(atac_small),
#'   pwm = pwm[1:10],
#'   genome = BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
CreateMotifMatrix <- function(
  features,
  pwm,
  genome,
  score = FALSE,
  use.counts = FALSE,
  ...
) {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop("Please install motifmatchr.
         https://www.bioconductor.org/packages/motifmatchr/")
  }

  # genome can be string
  if (is.character(x = genome)) {
    if (!requireNamespace("BSgenome", quietly = TRUE)) {
      stop("Please install BSgenome.
             https://www.bioconductor.org/packages/BSgenome/")
    }
    genome <- BSgenome::getBSgenome(genome = genome)
  }

  # check that all seqnames in features are in genome
  # remove missing, replace later with zeros and show warning
  miss_sn <- !(as.character(seqnames(x = features)) %in% seqlevels(x = genome))
  if (sum(miss_sn) > 0) {
    warning("Not all seqlevels present in supplied genome",
      immediate. = TRUE
    )
    # remove from features and remember original order
    feature_order <- features
    features <- features[!miss_sn]
  }
  motif_ix <- motifmatchr::matchMotifs(
    pwms = pwm,
    subject = features,
    genome = genome,
    out = "scores",
    ...
  )
  if (score) {
    motif.matrix <- motifmatchr::motifScores(object = motif_ix)
  } else {
    if (use.counts) {
      motif.matrix <- motifmatchr::motifCounts(object = motif_ix)
    } else {
      motif.matrix <- motifmatchr::motifMatches(object = motif_ix)
      motif.matrix <- as(Class = "CsparseMatrix", object = motif.matrix)
    }
  }
  rownames(motif.matrix) <- as.character(x = features)
  if (is.null(x = names(x = pwm))) {
    warning(
      "No 'names' attribute found in PFMatrixList. ",
      "Extracting names from individual entries."
    )
    colnames(x = motif.matrix) <- vapply(
      X = pwm, FUN = slot, FUN.VALUE = "character", "name"
    )
  }
  # add missing features
  if (sum(miss_sn) > 0) {
    replacement_matrix <- sparseMatrix(
      i = sum(miss_sn),
      j = ncol(x = motif.matrix)
    )
    rownames(x = replacement_matrix) <- as.character(x = feature_order[miss_sn])
    colnames(x = replacement_matrix) <- colnames(x = motif.matrix)
    motif.matrix <- rbind(motif.matrix, replacement_matrix)
    motif.matrix <- motif.matrix[as.character(x = feature_order), ]
  }
  return(motif.matrix)
}

#' @importFrom SeuratObject LayerData CreateAssayObject DefaultLayer
#' @importFrom Matrix rowSums
#'
#' @concept motifs
#' @method RunChromVAR GRangesAssay
#' @rdname RunChromVAR
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RunChromVAR(object = atac_small[["peaks"]], genome = BSgenome.Hsapiens.UCSC.hg19)
#' }
RunChromVAR.GRangesAssay <- function(
  object,
  genome,
  layer = NULL,
  motif.matrix = NULL,
  verbose = TRUE,
  ...
) {
  layer <- layer %||% DefaultLayer(object = object)
  if (!requireNamespace("chromVAR", quietly = TRUE)) {
    stop("Please install chromVAR. https://greenleaflab.github.io/chromVAR/")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment")
  }
  motif.matrix <- motif.matrix %||% GetMotifData(object = object, slot = "data")
  peak.matrix <- LayerData(object = object, layer = layer)
  peak.matrix <- as(object = peak.matrix, Class = "sparseMatrix")
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
  # Remove NA values https://github.com/GreenleafLab/chromVAR/issues/26
  row.data <- data.frame(SummarizedExperiment::rowData(x = chromvar.obj))
  row.data[is.na(x = row.data)] <- 0
  SummarizedExperiment::rowData(x = chromvar.obj) <- row.data
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
  obj <- CreateAssay5Object(data = chromvar.z)
  return(obj)
}

#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param layer Name of layer to use. If NULL, use the default layer.
#' @param new.assay.name Name of new assay used to store the chromVAR results.
#' Default is "chromvar".
#' @method RunChromVAR Seurat
#' @rdname RunChromVAR
#' @export
#' @importFrom SeuratObject DefaultAssay
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
  layer = NULL,
  new.assay.name = "chromvar",
  ...
) {
  if (!requireNamespace("chromVAR", quietly = TRUE)) {
    stop("Please install chromVAR. https://greenleaflab.github.io/chromVAR/")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Please install SummarizedExperiment")
  }
  assay <- assay %||% DefaultAssay(object = object)
  chromvar.assay <- RunChromVAR(
    object = object[[assay]],
    layer = layer,
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
#' Find motifs over-represented in a given set of genomic features.
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
#'  [RegionStats()] function. If NULL, use all features in the assay.
#' @param verbose Display messages
#' @param p.adjust.method Multiple testing correction method to be applied.
#' Passed to [stats::p.adjust()].
#' @param ... Arguments passed to [MatchRegionStats()].
#'
#' @return Returns a data frame
#'
#' @importFrom Matrix colSums
#' @importFrom stats phyper p.adjust
#' @importFrom methods is
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
  p.adjust.method = "BH",
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  background <- background %||% rownames(x = object)
  # TODO allow FindMotifs on a ChromatinAssay5 object
  if (!inherits(x = object[[assay]], what = "GRangesAssay")) {
    stop("Cannot run FindMotifs on ", class(x = object[[assay]]))
  }
  if (is(object = background, class2 = "numeric")) {
    if (verbose) {
      message(
        "Selecting background regions to match input ",
        "sequence characteristics"
      )
    }
    meta.feature <- object[[assay]][[]]
    mf.choose <- meta.feature[
      setdiff(x = rownames(x = meta.feature), y = features), ,
      drop = FALSE
    ]
    missing.features <- setdiff(x = features, y = rownames(x = meta.feature))
    if (length(x = missing.features) > 0) {
      warning(
        "The following features were not found in the assay: ",
        missing.features,
        "\nRemoving missing features",
        immediate. = TRUE
      )
      features <- intersect(x = features, y = rownames(x = meta.feature))
    }
    mf.query <- meta.feature[features, , drop = FALSE]
    
    # we can run FindMotifs on some object that does not have genomic ranges as the features
    # but the user will have to supply the meta.features containing GC.percent
    # to identify the matched background set of features
    background <- MatchRegionStats(
      meta.feature = mf.choose,
      query.feature = mf.query,
      regions = features,
      n = background,
      verbose = verbose,
      ...
    )
  }
  if (verbose) {
    msg <- ifelse(
      test = length(x = features) > 1,
      yes = " regions",
      no = " region"
    )
    message("Testing motif enrichment in ", length(x = features), msg)
  }
  if (length(x = features) < 10) {
    warning(
      "Testing motif enrichment using a small number of regions is ",
      "not recommended"
    )
  }
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  motif.names <- GetMotifData(
    object = object, assay = assay, slot = "motif.names"
  )
  
  # TODO update this for new ChromatinAssay5 class definition
  # no longer require that features in motif matrix match features in the data layer
  # will need to check that all features are in the motif object
  
  query.motifs <- motif.all[features, , drop = FALSE]
  background.motifs <- motif.all[background, , drop = FALSE]
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
    p.adjust = p.adjust(p = p.list, method = p.adjust.method),
    stringsAsFactors = FALSE
  )
  if (nrow(x = results) == 0) {
    return(results)
  } else {
    return(results[order(results[, 7], -results[, 6]), ])
  }
}

#' @param name A vector of motif names
#' @param id A vector of motif IDs. Only one of `name` and `id` should
#' be supplied
#' @rdname ConvertMotifID
#' @concept motifs
#' @importFrom methods hasArg
#' @export
ConvertMotifID.default <- function(object, name, id, ...) {
  if (hasArg(name = name) && hasArg(name = id)) {
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

#' @method ConvertMotifID GRangesAssay
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.GRangesAssay <- function(object, ...) {
  motifs <- Motifs(object = object)
  if (is.null(x = motifs)) {
    stop("No motif information present in assay")
  }
  return(ConvertMotifID(object = motifs, ...))
}

#' @method ConvertMotifID Assay
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.Assay <- function(object, ...) {
  stop("Cannot run ConvertMotifID on a standard Assay object")
}

#' @method ConvertMotifID StdAssay
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.StdAssay <- function(object, ...) {
  stop("Cannot run ConvertMotifID on an Assay5 object")
}

#' @param assay For `Seurat` object. Name of assay to use.
#' If NULL, use the default assay
#'
#' @importFrom SeuratObject DefaultAssay
#'
#' @method ConvertMotifID Seurat
#' @rdname ConvertMotifID
#' @concept motifs
#' @export
ConvertMotifID.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(ConvertMotifID(object = object[[assay]], ...))
}

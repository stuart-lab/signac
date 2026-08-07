#' @include generics.R
#'
NULL

#' @rdname AddMotifs
#' @method AddMotifs default
#' @concept motifs
#' @importFrom methods slot
#' @importFrom Seqinfo seqlevels seqnames seqinfo seqinfo<-
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
  # Since motifmatchr::matchMotifs returns a GenomicRanges without seqinfo
  seqinfo(motif.positions) <- seqinfo(genome)[seqlevels(motif.positions)]
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
#' @importFrom SeuratObject SetAssayData
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
#'   opts = list(
#'     collection = "CORE",
#'     tax_group = "vertebrates",
#'     all_versions = FALSE
#'   )
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
  # features on seqlevels absent from the genome were dropped above; add them
  # back as all-zero rows and restore the original feature order
  if (sum(miss_sn) > 0) {
    motif.matrix <- PadMissingFeatures(
      motif.matrix = motif.matrix,
      feature_order = as.character(x = feature_order),
      missing.features = as.character(x = feature_order[miss_sn])
    )
  }
  return(motif.matrix)
}

# Re-insert dropped features into a motif match matrix as all-zero rows,
# restoring the original feature order. Features on seqlevels not present in the
# genome cannot be scored by motifmatchr and are removed before matching (see
# CreateMotifMatrix); they carry no motif hits, so they are added back as zeros.
# @param motif.matrix Sparse motif match matrix for the scored features.
# @param feature_order Character vector of all feature names in original order.
# @param missing.features Character vector of the dropped feature names.
# @return The motif matrix with all features present, ordered by feature_order.
PadMissingFeatures <- function(motif.matrix, feature_order, missing.features) {
  # an empty sparse matrix (no nonzero entries) of the correct dimensions
  # one all-zero row per missing feature
  replacement_matrix <- sparseMatrix(
    i = integer(length = 0L),
    j = integer(length = 0L),
    dims = c(length(x = missing.features), ncol(x = motif.matrix))
  )
  rownames(x = replacement_matrix) <- missing.features
  colnames(x = replacement_matrix) <- colnames(x = motif.matrix)
  motif.matrix <- rbind(motif.matrix, replacement_matrix)
  return(motif.matrix[feature_order, ])
}

#' @importFrom SeuratObject LayerData CreateAssayObject DefaultLayer as.sparse
#' @importFrom Matrix rowSums
#'
#' @concept motifs
#' @method RunChromVAR GRangesAssay
#' @rdname RunChromVAR
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RunChromVAR(
#'   object = atac_small[["peaks"]],
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
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
  if (inherits(x = peak.matrix, what = "IterableMatrix")) {
    peak.matrix <- as.sparse(x = peak.matrix)
  }
  idx.keep <- rowSums(x = peak.matrix) > 0
  peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
  peak.ranges <- granges(x = object)
  peak.ranges <- peak.ranges[idx.keep]
  # a ChromatinAssay5 does not require the motif matrix to mirror the data
  # layer, so align on feature name rather than assuming the rows correspond
  missing.features <- setdiff(
    x = rownames(x = peak.matrix), y = rownames(x = motif.matrix)
  )
  if (length(x = missing.features) > 0) {
    stop(
      length(x = missing.features),
      " features are not present in the motif matrix, for example: ",
      paste(head(x = missing.features, n = 3), collapse = ", ")
    )
  }
  motif.matrix <- motif.matrix[rownames(x = peak.matrix), , drop = FALSE]
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

#' Resolve motifs sharing the same name
#'
#' Given a list of motifs and a quality score for each motif (lower is better),
#' either make all names unique by appending a suffix, or retain only the
#' best-scoring motif for each name.
#'
#' @param motifs A list of motif objects
#' @param names Motif names (may contain duplicates)
#' @param quality Numeric vector giving the quality rank of each motif, where
#'   lower values indicate higher-quality motifs. Ties are broken by the order
#'   the motifs appear in `motifs`.
#' @param keep Either "all" or "best"
#' @return The named list of motifs, in the original order
#' @keywords internal
#' @noRd
ResolveMotifNames <- function(motifs, names, quality, keep = c("all", "best")) {
  keep <- match.arg(arg = keep)
  if (keep == "all") {
    names(x = motifs) <- make.unique(names = names)
    return(motifs)
  }
  # order by name, then quality, then original position: the first entry
  # within each name group is the one to retain
  idx <- order(names, quality, seq_along(along.with = names))
  best <- idx[!duplicated(x = names[idx])]
  best <- sort(x = best)
  motifs <- motifs[best]
  names(x = motifs) <- names[best]
  return(motifs)
}

#' Parse motif IDs
#'
#' Extract the transcription factor name and relative motif quality from motif
#' IDs. Both are encoded in the ID, but the encoding differs between motif
#' collections:
#'
#' * `mode = "hocomoco"`: IDs are structured as
#'   `<TF>.<collection>.<model index>.[<source>.]<quality>`, for example
#'   `ALX3.H14CORE.0.SM.B` or `AHR_HUMAN.H11MO.0.B` (v11 IDs omit the source
#'   field, so the number of fields is not fixed). The quality letter is the
#'   final field and the model index is the third field. Per
#'   <https://hocomoco14.autosome.org/help>, models carry "a quality rating
#'   from A to D where A represents motifs with the highest confidence", and
#'   since v11 alternative motifs for a TF are "ranked from 0 (the primary
#'   model) to 1,2,.. (the alternative motifs)", where rank 0 models "are the
#'   most 'general' variants with the best performance across available data".
#'   The source field, when present, is any combination of the experiment-type
#'   abbreviations `P`, `S`, `M`, `G`, `I` and `B`.
#' * `mode = "jaspar"`: IDs are structured as `<base ID>.<version>`, for
#'   example `MA1104.2`. Higher versions are more recent. JASPAR IDs do not
#'   contain the TF name, which is instead given as a separate field on the
#'   header line.
#'
#' @param ids Character vector of motif IDs
#' @param mode Which convention to assume, either "jaspar" or "hocomoco". IDs
#'   that cannot be parsed under the requested convention carry no quality
#'   information and are ranked last.
#' @return A list with elements:
#'
#'   * `name`: the first `.`-separated field of each ID
#'   * `ranked`: whether any quality information could be extracted
#'   * `letter`: HOCOMOCO quality letter as an integer, or 0
#'   * `index`: HOCOMOCO model index or negated JASPAR version, so that lower
#'     values are always better
#' @keywords internal
#' @noRd
ParseMotifID <- function(ids, mode = c("jaspar", "hocomoco")) {
  mode <- match.arg(arg = mode)
  parts <- strsplit(x = ids, split = ".", fixed = TRUE)
  nfield <- lengths(x = parts)
  field <- function(i) {
    vapply(
      X = seq_along(along.with = parts), FUN.VALUE = character(1L),
      FUN = function(j) {
        pos <- if (i < 0) nfield[[j]] + 1 + i else i
        if (pos < 1 || pos > nfield[[j]]) NA_character_ else parts[[j]][[pos]]
      }
    )
  }
  last <- field(i = -1L)
  # quality is documented as running from A to D only, so a final field
  # outside that range is a source abbreviation or some other convention
  # rather than a quality rating
  letter <- match(x = toupper(x = last), table = LETTERS[1:4])
  index <- suppressWarnings(expr = as.numeric(x = field(i = 3L)))
  version <- suppressWarnings(expr = as.numeric(x = last))
  none <- rep_len(x = FALSE, length.out = length(x = ids))
  hocomoco <- if (mode == "hocomoco") !is.na(x = letter) else none
  jaspar <- if (mode == "jaspar") {
    nfield >= 2 & !is.na(x = version)
  } else {
    none
  }
  list(
    name = field(i = 1L),
    ranked = hocomoco | jaspar,
    letter = ifelse(test = hocomoco, yes = letter, no = 0),
    index = ifelse(
      test = hocomoco,
      yes = ifelse(test = is.na(x = index), yes = Inf, no = index),
      # higher JASPAR version is better, so negate
      no = ifelse(test = jaspar, yes = -version, no = 0)
    )
  )
}

#' Rank motif quality from the motif ID
#'
#' @param ids Character vector of motif IDs
#' @param names Motif names, used to detect whether any choice between motifs
#'   is actually needed
#' @param mode Which convention to assume, either "jaspar" or "hocomoco"
#' @return Integer vector of quality ranks, where lower values are better.
#'   IDs carrying no quality information are ranked last, and ties are broken
#'   by the position of the ID in `ids`.
#' @keywords internal
#' @noRd
MotifQuality <- function(ids, names, mode = c("jaspar", "hocomoco")) {
  mode <- match.arg(arg = mode)
  parsed <- ParseMotifID(ids = ids, mode = mode)
  # a motif has to be chosen for each duplicated name, so warn if none of the
  # IDs involved carry quality information: this usually means the wrong
  # convention was requested
  duplicated_names <- names %in% names[duplicated(x = names)]
  if (any(duplicated_names) && !any(parsed$ranked[duplicated_names])) {
    warning(
      "No motif quality information could be read from the motif IDs using ",
      "mode = \"", mode, "\". The first motif for each name will be kept. ",
      "Check that the correct value of mode is set for this motif collection.",
      call. = FALSE
    )
  }
  order(order(!parsed$ranked, parsed$letter, parsed$index))
}

#' Determine a motif name from its ID
#'
#' The first `.`-separated field of the ID holds the TF name for HOCOMOCO IDs
#' and the base matrix ID for JASPAR IDs, so this does not depend on which
#' convention the ID follows.
#'
#' @param ids Character vector of motif IDs
#' @param short_names Use the first `.`-separated field of the ID as the name
#' @return Character vector of motif names
#' @keywords internal
#' @noRd
MotifNameFromID <- function(ids, short_names) {
  if (!short_names) {
    return(ids)
  }
  vapply(
    X = strsplit(x = ids, split = ".", fixed = TRUE),
    FUN = function(x) x[[1L]],
    FUN.VALUE = character(1L)
  )
}

#' Read PWM files into a PWMatrixList
#'
#' Read position weight matrices from a directory of `.pwm` files
#' and return as a [TFBSTools::PWMatrixList]. Each `.pwm` file
#' should contain a header line starting with `>` followed by the motif
#' ID, and subsequent lines containing a position x 4 nucleotide matrix
#' (columns: A, C, G, T).
#'
#' @param pwm_dir Path to directory containing `.pwm` files
#' @param short_names Use the first section of the motif ID (before the first
#'   `.`) as the motif name. For example, `AHR.H14CORE.0.P.B` becomes `AHR`.
#'   Default is `TRUE`.
#' @param keep How to handle motifs that share the same name. This commonly
#'   happens when `short_names = TRUE`, as many collections contain several
#'   motifs per transcription factor. Options are:
#'
#'   * `"all"`: retain all motifs, appending a numeric suffix to duplicated
#'     names (`AHR`, `AHR.1`, `AHR.2`). This is the default.
#'   * `"best"`: retain only the highest-quality motif for each name. Quality
#'     is taken from the HOCOMOCO motif ID, using the quality letter (A is
#'     best, D is worst) and then the model index (lower is better). For
#'     example, `AHR.H14CORE.0.P.A` is preferred over both
#'     `AHR.H14CORE.1.P.A` and `AHR.H14CORE.0.P.B`. Both the four-field v11 ID
#'     format (`AHR_HUMAN.H11MO.0.B`) and the five-field format used from v12
#'     onwards are recognised. Motifs with IDs that do not follow this
#'     convention are ranked last.
#' @param mode Which motif ID convention to use when ranking motif quality for
#'   `keep = "best"`. Has no effect when `keep = "all"`.
#' @return A [TFBSTools::PWMatrixList]
#' @export
#' @concept motifs
ReadPWM <- function(
  pwm_dir,
  short_names = TRUE,
  keep = c("all", "best"),
  mode = c("hocomoco", "jaspar")
) {
  keep <- match.arg(arg = keep)
  mode <- match.arg(arg = mode)
  if (!requireNamespace("TFBSTools", quietly = TRUE)) {
    stop("Please install TFBSTools.
         https://www.bioconductor.org/packages/TFBSTools/")
  }
  pwm_files <- list.files(
    path = pwm_dir, pattern = "\\.pwm$", full.names = TRUE
  )
  pwm_list <- lapply(X = pwm_files, FUN = function(f) {
    lines <- readLines(con = f)
    motif_id <- sub(pattern = "^>", replacement = "", x = lines[1])
    motif_name <- MotifNameFromID(ids = motif_id, short_names = short_names)
    mat <- do.call(what = rbind, args = lapply(
      X = lines[-1],
      FUN = function(l) {
        as.numeric(x = strsplit(x = trimws(x = l), split = "\\s+")[[1]])
      }
    ))
    colnames(x = mat) <- c("A", "C", "G", "T")
    mat <- t(x = mat)
    TFBSTools::PWMatrix(
      ID = motif_id,
      name = motif_name,
      profileMatrix = mat
    )
  })
  raw_names <- vapply(
    X = pwm_list, FUN = TFBSTools::name, FUN.VALUE = character(1L)
  )
  raw_ids <- vapply(
    X = pwm_list, FUN = TFBSTools::ID, FUN.VALUE = character(1L)
  )
  pwm_list <- ResolveMotifNames(
    motifs = pwm_list,
    names = raw_names,
    quality = MotifQuality(ids = raw_ids, names = raw_names, mode = mode),
    keep = keep
  )
  do.call(what = TFBSTools::PWMatrixList, args = pwm_list)
}

#' Read JASPAR-format PFMs and convert to PWMs
#'
#' Read position frequency matrices (PFMs) from a JASPAR-format file and
#' convert to position weight matrices (PWMs). Each motif entry should have
#' a header line starting with `>` followed by 4 rows (A, C, G, T). Rows
#' may optionally include nucleotide labels and brackets
#' (e.g. `A  [ 4 19 0 0 ]`).
#'
#' The header line should contain the matrix ID, optionally followed by
#' whitespace and the transcription factor name (for example
#' `>MA1104.2 GATA6`). The full matrix ID is stored as the motif ID, and the
#' transcription factor name is used to name the returned motifs. Collections
#' distributed in JASPAR format do not always include the transcription factor
#' name as a separate field: HOCOMOCO, for example, encodes it in the matrix
#' ID (`>ALX3.H14CORE.0.SM.B`). In that case the name is taken from the ID
#' according to `short_names`.
#'
#' @param file Path to JASPAR-format PFM file
#' @param pseudocount Pseudocount added during PFM to PWM conversion
#' @param short_names Use the first section of the motif ID (before the first
#'   `.`) as the motif name when the header line does not include a separate
#'   transcription factor name. For example, `ALX3.H14CORE.0.SM.B` becomes
#'   `ALX3` and `MA0004.1` becomes `MA0004`. If `FALSE`, the full matrix ID is
#'   always used as the name. Default is `TRUE`.
#' @param keep How to handle motifs that share the same name, which happens
#'   when a file contains several matrices for one transcription factor.
#'   Options are:
#'
#'   * `"all"`: retain all motifs, appending a numeric suffix to duplicated
#'     names (`ALX3`, `ALX3.1`). This is the default.
#'   * `"best"`: retain only the highest-quality motif for each name. Quality
#'     is taken from the matrix ID, interpreted according to `mode`. With
#'     `mode = "jaspar"` motifs are ranked on the matrix version (higher is
#'     better), so `MA1104.2` is preferred over `MA1104.1`. With
#'     `mode = "hocomoco"` they are ranked on the quality letter (A is best,
#'     D is worst) and then the model index (lower is better), so
#'     `ALX3.H14CORE.0.SM.B` is preferred over `ALX3.H14CORE.1.S.B`. Motifs
#'     with IDs carrying no quality information are ranked last.
#' @param mode Which motif ID convention to use when ranking motif quality for
#'   `keep = "best"`. The default is `"jaspar"`; set `"hocomoco"` when reading
#'   a HOCOMOCO collection distributed in JASPAR format, since those IDs
#'   encode a quality rating rather than a version. A warning is given if no
#'   quality information can be read from the IDs. Has no effect when
#'   `keep = "all"`.
#' @return A [TFBSTools::PWMatrixList]
#' @importFrom methods is
#' @export
#' @concept motifs
ReadJASPAR <- function(
  file,
  pseudocount = 1,
  short_names = TRUE,
  keep = c("all", "best"),
  mode = c("jaspar", "hocomoco")
) {
  if (!requireNamespace("TFBSTools", quietly = TRUE)) {
    stop("Please install TFBSTools.
         https://www.bioconductor.org/packages/TFBSTools/")
  }
  keep <- match.arg(arg = keep)
  mode <- match.arg(arg = mode)
  lines <- trimws(x = readLines(con = file))
  motif_indices <- grep(pattern = "^>", x = lines)
  motifs <- list()
  motif_names <- character()
  motif_ids <- character()
  for (i in seq_along(along.with = motif_indices)) {
    start <- motif_indices[i]
    end <- if (i < length(x = motif_indices)) {
      motif_indices[i + 1] - 1
    } else {
      length(x = lines)
    }
    header <- strsplit(
      x = sub(pattern = "^>", replacement = "", x = lines[start]),
      split = "\\s+"
    )[[1]]
    id <- header[[1]]
    # a separate TF name field takes precedence over the ID, but only when
    # short names were requested: otherwise the full ID is used as the name
    name <- if (length(x = header) > 1 && short_names) {
      header[[2]]
    } else {
      MotifNameFromID(ids = id, short_names = short_names)
    }
    if (end <= start) {
      stop("Motif ", id, " has no matrix rows")
    }
    mat_lines <- lines[(start + 1):end]
    mat <- do.call(what = rbind, args = lapply(
      X = mat_lines,
      FUN = function(x) {
        # strip nucleotide labels and brackets (e.g. "A  [ 1 2 3 ]")
        x <- gsub(pattern = "[A-Za-z]|\\[|\\]", replacement = "", x = x)
        as.numeric(x = strsplit(x = trimws(x = x), split = "\\s+")[[1]])
      }
    ))
    if (nrow(x = mat) != 4) {
      stop(
        "Motif ", id, " has ", nrow(x = mat),
        " rows; expected 4 (A,C,G,T)."
      )
    }
    rownames(x = mat) <- c("A", "C", "G", "T")
    pfm <- TFBSTools::PFMatrix(
      ID = id, name = name, profileMatrix = mat
    )
    pwm <- suppressWarnings(
      expr = TFBSTools::toPWM(x = pfm, pseudocounts = pseudocount)
    )
    if (is(object = pwm, class2 = "PWMatrix")) {
      motifs[[length(x = motifs) + 1]] <- pwm
      motif_names <- c(motif_names, name)
      motif_ids <- c(motif_ids, id)
    }
  }
  motifs <- ResolveMotifNames(
    motifs = motifs,
    names = motif_names,
    quality = MotifQuality(ids = motif_ids, names = motif_names, mode = mode),
    keep = keep
  )
  do.call(what = TFBSTools::PWMatrixList, args = motifs)
}

globalVariables(names = "pvalue", package = "Signac")
#' Find over-represented motifs in genomic regions
#'
#' Identify DNA sequence motifs that are over-represented in a set of genomic
#' features (for example, a set of differentially accessible peaks) relative to
#' a background set of features. Enrichment is quantified with a one-sided
#' hypergeometric test.
#'
#' @details
#' For each motif, `FindMotifs` compares how often the motif occurs in the
#' query features to how often it occurs in a background set of features. The
#' motif occurrences themselves are taken from the `motifs` data of the assay
#' (a binary feature-by-motif matrix created with [AddMotifs()]). The procedure is:
#'
#' 1. **Restrict to scored features.** Query and background features that are
#'    not present in the motif matrix are dropped, since motif occurrences are
#'    only known for the features that were scored.
#' 2. **Select the background set.** If `background` is a single number, that
#'    many features are selected to match the sequence characteristics (by
#'    default GC content) of the query using [MatchRegionStats()], drawn from
#'    features other than the query. If `background` is a vector of feature
#'    names it is used directly; if it is `NULL`, all features in the assay are
#'    used. Matching the background to the query's sequence composition is
#'    important, since a motif can appear enriched simply because the query
#'    features differ in base composition from the genomic average (for example,
#'    a GC-rich query will appear enriched for GC-rich motifs) rather than
#'    because of genuine biological enrichment.
#' 3. **Count motif occurrences.** For every motif the function counts the
#'    number of query features containing it (`observed`) and the number of
#'    background features containing it (`background`).
#' 4. **Compute the enrichment p-value.** The query and background features
#'    together define a population of features (their union). For each motif,
#'    let `q` be the number of query features containing the motif, `k` the
#'    number of query features, `m` the number of features in the population
#'    containing the motif, and `N` the total number of features in the
#'    population. The one-sided hypergeometric p-value is
#'    `phyper(q - 1, m, N - m, k, lower.tail = FALSE)`, the probability of
#'    observing at least `q` motif-containing features when `k` features are
#'    drawn at random from the population. Defining the population as the union
#'    of the query and background guarantees that `m >= q` and `N >= k`, so the
#'    test is always well defined whether or not the background overlaps the
#'    query.
#'    
#' 5. **Summarize and correct for multiple testing.** A fold enrichment is
#'    reported as the percentage of query features containing the motif divided
#'    by the percentage of background features containing the motif, and the
#'    p-values are adjusted across motifs with [stats::p.adjust()] using
#'    `p.adjust.method`. Results are ordered by increasing p-value, breaking
#'    ties by decreasing fold enrichment.
#'
#' @param object A Seurat object
#' @param features A vector of features to test for enrichment over the
#' background set. These should be present in the `motifs` data of the assay.
#' @param assay Which assay to use. Default is the active assay
#' @param background Either a vector of features to use as the background set,
#' or a number specifying the number of features to select as a background set.
#' If a number is provided, regions will be selected to match the sequence
#' characteristics of the query features. To match the sequence characteristics,
#' these characteristics must be stored in the feature metadata for the assay.
#' This can be added using the [RegionStats()] function. If `NULL`, use all
#' features in the assay.
#' @param verbose Display messages
#' @param p.adjust.method Multiple testing correction method to be applied.
#' Passed to [stats::p.adjust()].
#' @param ... Arguments passed to [MatchRegionStats()].
#'
#' @return Returns a data frame with one row per motif and the following
#' columns:
#'   - `motif`: the motif ID
#'   - `observed`: number of query features containing the motif
#'   - `background`: number of background features containing the motif
#'   - `percent.observed`: percentage of query features containing the motif
#'   - `percent.background`: percentage of background features containing the
#'   motif
#'   - `fold.enrichment`: `percent.observed` divided by `percent.background`
#'   - `pvalue`: the hypergeometric enrichment p-value
#'   - `motif.name`: the motif name
#'   - `p.adjust`: the p-value adjusted for multiple testing
#'
#'  Rows are ordered by increasing `pvalue`, breaking ties by decreasing
#'  `fold.enrichment`.
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
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("Cannot run FindMotifs on ", class(x = object[[assay]]))
  }
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  motif.names <- GetMotifData(
    object = object, assay = assay, slot = "motif.names"
  )
  motif.features <- rownames(x = motif.all)

  # features and background must correspond to rows of the motif matrix.
  # ChromatinAssay5 does not require the motif matrix to mirror the data layer,
  # so drop anything not present in the motif matrix before subsetting.
  missing.query <- setdiff(x = features, y = motif.features)
  if (length(x = missing.query) > 0) {
    warning(
      "The following features are not in the motif matrix ",
      "and will be ignored: ",
      paste(missing.query, collapse = ", "),
      immediate. = TRUE
    )
    features <- intersect(x = features, y = motif.features)
  }
  if (length(x = features) == 0) {
    stop("No query features are present in the motif matrix")
  }

  if (is(object = background, class2 = "numeric")) {
    if (verbose) {
      message(
        "Selecting background regions to match input ",
        "sequence characteristics"
      )
    }
    meta.feature <- object[[assay]][[]]
    # only sample candidates that are present in the motif matrix
    meta.feature <- meta.feature[
      intersect(x = rownames(x = meta.feature), y = motif.features), ,
      drop = FALSE
    ]
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

    # we can run FindMotifs on some object that does not have genomic ranges as
    # the features but the user will have to supply the meta.features containing
    # GC.percent to identify the matched background set of features
    background <- MatchRegionStats(
      meta.feature = mf.choose,
      query.feature = mf.query,
      regions = features,
      n = background,
      verbose = verbose,
      ...
    )
  } else {
    missing.bg <- setdiff(x = background, y = motif.features)
    if (length(x = missing.bg) > 0) {
      warning(
        "The following background features are not in the motif matrix ",
        "and will be ignored: ",
        paste(missing.bg, collapse = ", "),
        immediate. = TRUE
      )
      background <- intersect(x = background, y = motif.features)
    }
    if (length(x = background) == 0) {
      stop("No background features are present in the motif matrix")
    }
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

  query.motifs <- motif.all[features, , drop = FALSE]
  background.motifs <- motif.all[background, , drop = FALSE]
  query.counts <- colSums(x = query.motifs)
  background.counts <- colSums(x = background.motifs)
  percent.observed <- query.counts / length(x = features) * 100
  percent.background <- background.counts / length(x = background) * 100
  fold.enrichment <- percent.observed / percent.background
  test.features <- union(x = features, y = background)
  test.counts <- colSums(x = motif.all[test.features, , drop = FALSE])
  p.list <- phyper(
    q = query.counts - 1,
    m = test.counts,
    n = length(x = test.features) - test.counts,
    k = length(x = features),
    lower.tail = FALSE
  )
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
    return(results[
      order(results$pvalue, -results$fold.enrichment),
    ])
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
  } else if (!hasArg(name = name) && !(hasArg(name = id))) {
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

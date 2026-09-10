################################################################################
### Exported functions
################################################################################


#' Clear cached AlphaGenome Atlas results
#'
#' Discard the AlphaGenome Atlas API results cached in memory by
#' [LoadAVIFromAtlas()] during the current R session.
#'
#' @return Invisibly returns `NULL`
#' @seealso [LoadAVIFromAtlas()]
#' @export
#' @concept alphagenome
#' @examples
#' ClearAtlasCache()
ClearAtlasCache <- function() {
  rm(
    list = ls(envir = atlas.cache, all.names = TRUE),
    envir = atlas.cache
  )
  invisible(x = NULL)
}

#' List AlphaGenome Atlas variant scorers
#'
#' List the variant scorers available through the AlphaGenome Atlas API,
#' along with whether each scorer produces signed scores and the number of
#' scores (tracks or features) it returns per variant. See
#' [LoadAVIFromAtlas()] for the software and API key requirements.
#'
#' @param api.key AlphaGenome API key. If `NULL`, the
#' `ALPHA_GENOME_API_KEY` environment variable is used.
#'
#' @return Returns a `data.frame` with columns `name`,
#' `is_signed`, and `n_scores`.
#' @seealso [LoadAVIFromAtlas()]
#' @export
#' @concept alphagenome
#' @examples
#' \dontrun{
#' ListAtlasScorers()
#' }
ListAtlasScorers <- function(api.key = NULL) {
  client <- AtlasClient(api.key = api.key)
  return(AtlasScorerTable(client = client))
}

#' Query AVI scores from the AlphaGenome Atlas API
#'
#' Retrieve AlphaGenome Variant Impact (AVI) scores or AVI SHAP feature
#' attributions for all SNVs in a genomic region directly from the AlphaGenome
#' Atlas API, without downloading the Atlas files.
#'
#' The Atlas API is a gRPC service that is accessed through the `alphagenome`
#' Python package. This function uses reticulate to call the Python client, so
#' both reticulate and the `alphagenome` Python package must be installed.
#'
#' An AlphaGenome API key is required (see
#' <https://deepmind.google.com/science/alphagenome>). The key can be
#' passed using the `api.key` argument or stored in the
#' `ALPHA_GENOME_API_KEY` environment variable (for example in your
#' `~/.Renviron` file or by running `Sys.setenv(ALPHA_GENOME_API_KEY = ...)`.
#'
#' Each position in the region carries three SNVs. The Atlas API serves these
#' in 32 bp chunks (about 0.5 s per request), so querying large regions can
#' take some time.
#'
#' Results are cached in memory for the duration of the R session
#' (`cache = TRUE`), so plotting the same or an overlapping region again
#' only queries the API for the parts of the region that have not been
#' fetched before. Use [ClearAtlasCache()] to discard the cached results.
#'
#' The Atlas API responds with "The service is currently unavailable" when a
#' key exceeds its request rate. Requests that fail this way are retried with
#' exponential backoff (shared across the concurrent workers) up to
#' `max.retries` times before the query fails. Lowering
#' `max.workers` reduces the request rate.
#'
#' @param region Genomic region ([GenomicRanges::GRanges] or a string that can
#' be converted to `GRanges`, such as "chr11:5225000-5230000").
#' @param api.key AlphaGenome API key. If `NULL`, the
#' `ALPHA_GENOME_API_KEY` environment variable is used.
#' @param scorer Name of the Atlas variant scorer to retrieve. The default,
#' `"AVI_SCORE_FEATURE_IMPORTANCE"`, is the AVI SHAP feature attributions.
#' Other useful scorers are `"AVI_SCORE"` (the summary AVI score) and
#' `"AVI_SCORE_MODEL_FEATURES"` (the 18 raw AVI input features). Use
#' [ListAtlasScorers()] to see all available scorer names.
#' @param max.workers Number of requests sent to the Atlas API concurrently.
#' Lower this if the API responds with "The service is currently unavailable",
#' which can indicate rate limiting.
#' @param max.retries Maximum number of times a throttled request is retried
#' before the query fails (see Details).
#' @param cache Reuse results fetched earlier in the session and only query
#' the API for parts of the region that have not been fetched yet.
#' @param verbose Display progress messages
#'
#' @return Returns a `data.frame` with one row per variant, containing
#' the columns `chromosome`, `position`, `ref`, `alt`,
#' and one numeric column per score returned by the scorer (for AVI SHAP,
#' one column per feature). The data frame can be passed directly to
#' [AVITrack()] or the `avi` argument of [CoveragePlot()].
#'
#' @seealso [ListAtlasScorers()], [ClearAtlasCache()], [AVITrack()]
#' @references Cheng et al. (2026). AlphaGenome Atlas: in silico mutagenesis of
#' the entire human genome improves prioritization and interpretation of
#' non-coding variants
#' <https://deepmind.google.com/science/alphagenome/atlas>
#'
#' @importFrom GenomicRanges GRanges
#' @export
#' @concept alphagenome
#' @examples
#' \dontrun{
#' # requires an AlphaGenome API key in ALPHA_GENOME_API_KEY
#' avi <- LoadAVIFromAtlas(region = "chr11:5225000-5230000")
#' AVITrack(avi = avi, region = "chr11:5225000-5230000")
#' }
LoadAVIFromAtlas <- function(
  region,
  api.key = NULL,
  scorer = "AVI_SCORE_FEATURE_IMPORTANCE",
  max.workers = 4,
  max.retries = 5,
  cache = TRUE,
  verbose = TRUE
) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- GRanges(region)
  }
  if (length(x = region) != 1) {
    stop("region must describe a single genomic range")
  }
  client <- AtlasClient(api.key = api.key)
  avi <- QueryAtlasCached(
    client = client,
    region = region,
    scorer = scorer,
    cache = cache,
    max.workers = max.workers,
    max.retries = max.retries,
    verbose = verbose
  )
  return(avi)
}

################################################################################
### Internal functions
################################################################################

# Session-level cache of Atlas API results, one entry per scorer holding the
# fetched variants and the genomic ranges they cover
atlas.cache <- new.env(parent = emptyenv())


# Resolve the AlphaGenome API key
#
# @param api.key API key supplied by the user, or NULL
# @return A non-empty API key string
AlphaGenomeAPIKey <- function(api.key = NULL) {
  api.key <- api.key %||% Sys.getenv("ALPHA_GENOME_API_KEY")
  if (!nzchar(x = api.key)) {
    api.key <- Sys.getenv("ALPHAGENOME_API_KEY")
  }
  if (!is.character(x = api.key) || length(x = api.key) != 1 || !nzchar(x = api.key)) {
    stop(
      "An AlphaGenome API key is required. Supply it using the api.key ",
      "argument or set the ALPHA_GENOME_API_KEY environment variable. ",
      "See https://deepmind.google.com/science/alphagenome to obtain a key."
    )
  }
  return(api.key)
}

# Convert an Atlas AnnData result to a data.frame
#
# The Python Atlas client returns one AnnData object per scorer, with one row
# per variant (the Variant objects are stored in obs$variant) and one column
# per score, named in var$name.
#
# @param ad A Python anndata.AnnData object
# @param scorer Scorer name, used to name the score column when the scorer
# returns a single unnamed score per variant
# @return A data.frame with columns chromosome, position, ref, alt, and one
# numeric column per score
AtlasAnnDataToDataFrame <- function(ad, scorer = "score") {
  scores <- as.matrix(x = reticulate::py_to_r(x = ad$X))
  var <- as.data.frame(x = reticulate::py_to_r(x = ad$var))
  # inspect the obs column names without converting the (large) obs frame
  obs <- reticulate::py_get_attr(x = ad, name = "obs")
  obs.columns <- as.character(x = reticulate::py_to_r(
    x = reticulate::py_get_attr(x = obs, name = "columns")$tolist()
  ))
  score.names <- if ("name" %in% colnames(x = var)) {
    as.character(x = var[["name"]])
  } else {
    rownames(x = var)
  }
  if (length(x = score.names) != ncol(x = scores)) {
    score.names <- paste0(scorer, "_", seq_len(length.out = ncol(x = scores)))
  }
  if (length(x = score.names) == 1 && score.names %in% c("score", "0", "")) {
    score.names <- scorer
  }
  colnames(x = scores) <- score.names
  if (!("variant" %in% obs.columns)) {
    stop("Atlas result does not contain variant information")
  }
  # extract the variant fields in Python in one pass rather than one
  # reticulate call per field per variant
  avi <- reticulate::py_to_r(x = AtlasHelper()$variants_to_frame(ad))
  avi <- as.data.frame(x = avi, stringsAsFactors = FALSE)
  avi[["position"]] <- as.integer(x = avi[["position"]])
  avi[["chromosome"]] <- as.character(x = avi[["chromosome"]])
  avi[["ref"]] <- as.character(x = avi[["ref"]])
  avi[["alt"]] <- as.character(x = avi[["alt"]])
  avi <- cbind(avi, as.data.frame(x = scores, stringsAsFactors = FALSE))
  rownames(x = avi) <- NULL
  avi <- avi[order(avi[["position"]], avi[["alt"]]), , drop = FALSE]
  rownames(x = avi) <- NULL
  return(avi)
}

# Create an AlphaGenome Atlas API client
#
# @param api.key API key, or NULL to read from the environment
# @return A Python AtlasClient object
AtlasClient <- function(api.key = NULL) {
  api.key <- AlphaGenomeAPIKey(api.key = api.key)
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Please install reticulate to query the AlphaGenome Atlas API: ",
      "install.packages('reticulate')"
    )
  }
  if (!reticulate::py_module_available(module = "alphagenome")) {
    stop(
      "The alphagenome Python package is required to query the AlphaGenome ",
      "Atlas API. Install it with reticulate::py_install('alphagenome', ",
      "pip = TRUE), or configure reticulate to use a Python environment ",
      "where it is installed."
    )
  }
  atlas <- reticulate::import(module = "alphagenome.atlas.atlas")
  client <- tryCatch(
    expr = atlas$create(api_key = api.key),
    error = function(e) {
      stop("Could not connect to the AlphaGenome Atlas API: ", conditionMessage(e))
    }
  )
  return(client)
}

# Import the Python helper module shipped with Signac
AtlasHelper <- function() {
  # Keep package libraries read-only and restore the user's Python setting.
  sys <- reticulate::import(module = "sys")
  previous <- sys$dont_write_bytecode
  sys$dont_write_bytecode <- TRUE
  on.exit(sys$dont_write_bytecode <- previous, add = TRUE)
  reticulate::import_from_path(
    module = "signac_atlas",
    path = system.file("python", package = "Signac")
  )
}

# Summarize the scorers available from an Atlas client
#
# @param client A Python AtlasClient object
# @return A data.frame with columns name, is_signed, n_scores
AtlasScorerTable <- function(client) {
  metadata <- client$scorer_metadata()
  scorer.names <- names(x = metadata)
  is.signed <- vapply(
    X = scorer.names,
    FUN = function(x) isTRUE(x = metadata[[x]]$is_signed),
    FUN.VALUE = logical(length = 1)
  )
  n.scores <- vapply(
    X = scorer.names,
    FUN = function(x) {
      # read the shape without converting the pandas DataFrame to R
      tm <- reticulate::py_get_attr(x = metadata[[x]], name = "track_metadata")
      shape <- reticulate::py_to_r(x = reticulate::py_get_attr(x = tm, name = "shape"))
      as.integer(x = shape[[1]])
    },
    FUN.VALUE = integer(length = 1)
  )
  scorers <- data.frame(
    name = scorer.names,
    is_signed = unname(obj = is.signed),
    n_scores = unname(obj = n.scores),
    stringsAsFactors = FALSE
  )
  rownames(x = scorers) <- NULL
  return(scorers)
}

# Identify the AVI SHAP feature columns in an AVI table
#
# Feature attribution columns are the numeric columns whose names map to one
# of the known AVI modality groups (see AVIFeatureGroups()). Other numeric
# columns (aggregate scores such as AVI_SCORE, quality metrics, identifiers)
# are never treated as features, but can still be requested explicitly by
# name in AVITrack().
#
# @param avi A data.frame with chromosome, position, ref and alt columns
# @param allow.empty Return an empty vector rather than erroring when no
# feature columns are found
# @return A character vector of column names
AVIFeatureColumns <- function(avi, allow.empty = FALSE) {
  reserved <- c("chromosome", "position", "ref", "alt")
  candidates <- setdiff(x = colnames(x = avi), y = reserved)
  is.numeric.col <- vapply(
    X = avi[, candidates, drop = FALSE],
    FUN = is.numeric,
    FUN.VALUE = logical(length = 1)
  )
  candidates <- candidates[is.numeric.col]
  groups <- AVIFeatureGroups(features = candidates)
  features <- candidates[groups %in% AVIFeatureTable()[["group"]]]
  if (length(x = features) == 0 && !allow.empty) {
    stop(
      "No AVI SHAP feature columns found in the AVI data (numeric columns ",
      "matching a known modality group, see AVIFeatureGroups()). Found ",
      "numeric columns: ", paste(candidates, collapse = ", "),
      ". Score columns can be requested explicitly by name with `features`."
    )
  }
  return(features)
}

# Group AVI SHAP features by modality
#
# Assign each of the 18 AlphaGenome Atlas AVI input features to a modality
# group for display (see AVIFeatureTable()). Matching is exact apart from
# case; names that are not one of the 18 features are returned unchanged.
#
# @param features Character vector of feature names
# @return A named character vector of group labels, with names equal to the
# input feature names
AVIFeatureGroups <- function(features) {
  table <- AVIFeatureTable()
  idx <- match(x = toupper(x = features), table = table[["feature"]])
  groups <- ifelse(
    test = is.na(x = idx), yes = features, no = table[["group"]][idx]
  )
  names(x = groups) <- features
  return(groups)
}

# The 18 AVI input features, their modality groups and display colors
#
# Feature names are those of the AlphaGenome Atlas publication (Cheng et al.
# 2026, fig. S3 and Methods) and of the Atlas API. Groups follow the input
# feature categories of its Fig. 1, with the transcription readouts (CAGE,
# PRO-cap, RNA-seq) combined as "Expression" and AlphaMissense combined with
# the protein-termination consequences as "Protein".
#
# @return A data.frame with columns feature, group and color
AVIFeatureTable <- function() {
  groups <- c(
    "Accessibility" = "#EFE804",
    "TF binding" = "#17BECF",
    "Histone" = "#F4B9F5",
    "Expression" = "#555AC6",
    "Polyadenylation" = "#BCBD22",
    "Splicing" = "#F3A27F",
    "Contact maps" = "#8C564B",
    "Protein" = "#51606D",
    "Conservation" = "#81DAB6",
    "Indel" = "#E377C2"
  )
  features <- c(
    "MAX_ABS_ATAC" = "Accessibility",
    "MAX_ABS_DNASE" = "Accessibility",
    "MAX_ABS_CHIP_TF" = "TF binding",
    "MAX_ABS_CHIP_HISTONE" = "Histone",
    "MAX_ABS_CAGE" = "Expression",
    "MAX_ABS_PROCAP" = "Expression",
    "MAX_ABS_RNA_SEQ" = "Expression",
    "MAX_ABS_POLYADENYLATION" = "Polyadenylation",
    "MERGED_SPLICING" = "Splicing",
    "MAX_ABS_CONTACT_MAPS" = "Contact maps",
    "ALPHAMISSENSE" = "Protein",
    "PROTEIN_TERMINATION" = "Protein",
    "START_LOST" = "Protein",
    "STOP_LOST" = "Protein",
    "PHASTCONS_470_WAY" = "Conservation",
    "CACTUS_241_WAY" = "Conservation",
    "IS_INSERTION" = "Indel",
    "IS_DELETION" = "Indel"
  )
  data.frame(
    feature = names(x = features),
    group = unname(obj = features),
    color = unname(obj = groups[features]),
    stringsAsFactors = FALSE
  )
}

# Bin per-position AVI values across a region
#
# Positions without a variant are assigned a value of zero. The region is
# divided into `bins` equal-width bins and each bin takes the value of largest
# magnitude (keeping its sign) among its positions, as in the AlphaGenome Atlas
# publication. This preserves sharp SHAP peaks while keeping the number of
# plotted values fixed. With rank.by supplied, all tracks take their value from
# the same position per bin so that stacked contributions remain additive.
#
# @param values A data.frame with a position column and one column per track
# @param region A GRanges object (a single range)
# @param tracks Names of the track columns in `values`
# @param bins Number of bins across the region (at most one per position)
# @param rank.by Optional column used to select one common position per bin.
# @return A data.frame in long format with columns position (bin midpoint),
# group and value, with the bin width in bp stored in the "bin.size" attribute
#' @importFrom GenomicRanges start end width
BinAVI <- function(values, region, tracks, bins = 200, rank.by = NULL) {
  positions <- start(x = region):end(x = region)
  n <- length(x = positions)
  idx <- match(x = positions, table = values[["position"]])
  bins <- if (is.finite(x = bins)) min(as.integer(x = bins), n) else n
  bins <- max(1L, bins)
  bin.size <- n / bins
  bin <- pmin(floor(x = (positions - start(x = region)) / bin.size), bins - 1)
  bin.mid <- start(x = region) + (bin + 0.5) * bin.size - 0.5
  prepare <- function(track) {
    v <- values[[track]][idx]
    v[is.na(x = v)] <- 0
    return(v)
  }
  selected <- NULL
  if (!is.null(x = rank.by)) {
    o <- order(bin, -abs(x = prepare(rank.by)))
    selected <- o[!duplicated(x = bin[o])]
  }
  long <- list()
  for (track in tracks) {
    v <- prepare(track)
    keep <- selected
    if (is.null(x = keep)) {
      o <- order(bin, -abs(x = v))
      keep <- o[!duplicated(x = bin[o])]
    }
    long[[track]] <- data.frame(
      position = bin.mid[keep],
      group = track,
      value = v[keep],
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(what = rbind, args = long)
  rownames(x = long) <- NULL
  attr(x = long, which = "bin.size") <- bin.size
  return(long)
}

# Check that an AVI table has the required variant columns
#
# @param avi A data.frame
# @return The data.frame with position coerced to integer
CheckAVIColumns <- function(avi) {
  required <- c("chromosome", "position", "ref", "alt")
  missing <- setdiff(x = required, y = colnames(x = avi))
  if (length(x = missing) > 0) {
    stop(
      "AVI data must contain the columns ", paste(required, collapse = ", "),
      ". Missing: ", paste(missing, collapse = ", ")
    )
  }
  avi[["position"]] <- as.integer(x = avi[["position"]])
  avi[["chromosome"]] <- as.character(x = avi[["chromosome"]])
  return(avi)
}

# Collapse the alternate alleles at each position
#
# The Atlas carries three SNVs at each position (one per alternate allele).
# Following the AlphaGenome Atlas publication, the variant with the largest
# absolute total AVI score is kept at each position. The total used for
# ranking is the "Total" column when present (the sum of all SHAP features), so
# that the selected allele does not depend on which features are displayed;
# otherwise the sum of the requested features is used and returned as "Total".
#
# @param avi A data.frame with chromosome, position, ref and alt columns
# @param features Character vector of feature columns
# @return A data.frame with one row per position containing the position, ref
# and alt alleles of the selected variant, the Total column, and the feature
# columns
CollapseAVIAlleles <- function(avi, features) {
  if (!("Total" %in% colnames(x = avi))) {
    avi[["Total"]] <- rowSums(
      x = as.matrix(x = avi[, features, drop = FALSE]), na.rm = TRUE
    )
  }
  o <- order(avi[["position"]], -abs(x = avi[["Total"]]))
  avi <- avi[o, , drop = FALSE]
  avi <- avi[!duplicated(x = avi[["position"]]), , drop = FALSE]
  keep <- intersect(
    x = c("position", "ref", "alt", unique(x = c(features, "Total"))),
    y = colnames(x = avi)
  )
  collapsed <- avi[, keep, drop = FALSE]
  rownames(x = collapsed) <- NULL
  return(collapsed)
}

# Query an Atlas client for the variants in a region
#
# @param client A Python AtlasClient object (or a compatible test double)
# @param region A GRanges object (a single range)
# @param scorer Scorer name
# @param max.workers Number of concurrent requests
# @param max.retries Maximum retries per throttled request
# @param verbose Display progress messages
# @return A data.frame with chromosome, position, ref, alt and score columns
#' @importFrom Seqinfo seqnames
#' @importFrom GenomicRanges start end width
QueryAtlas <- function(
  client,
  region,
  scorer = "AVI_SCORE_FEATURE_IMPORTANCE",
  max.workers = 4,
  max.retries = 5,
  verbose = TRUE
) {
  ValidateAtlasQuery(max.workers, max.retries)
  if (length(x = region) != 1) {
    stop("region must describe a single genomic range")
  }
  chromosome <- as.character(x = seqnames(x = region))
  if (verbose) {
    message(
      "Querying ", width(x = region) * 3, " variants for scorer '", scorer,
      "' from the AlphaGenome Atlas"
    )
  }
  helper <- AtlasHelper()
  ad <- tryCatch(
    expr = helper$query_interval_chunked(
      client = client,
      chromosome = chromosome,
      start = as.integer(x = start(x = region)),
      end = as.integer(x = end(x = region)),
      scorer = scorer,
      max_workers = as.integer(x = max.workers),
      max_retries = as.integer(x = max.retries)
    ),
    error = function(e) {
      hint <- ""
      if (grepl(pattern = "UNAVAILABLE", x = conditionMessage(e))) {
        hint <- paste0(
          "\nThe Atlas API refused the request after ", max.retries,
          " retries. Try again later, try a smaller region, or lower ",
          "max.workers (currently ", max.workers, ")."
        )
      }
      stop("AlphaGenome Atlas query failed: ", conditionMessage(e), hint)
    }
  )
  if (is.null(x = ad)) {
    stop(
      "The AlphaGenome Atlas returned no scores for scorer '", scorer,
      "' in ", chromosome, ":", start(x = region), "-", end(x = region),
      ". Use ListAtlasScorers() to see the available scorer names."
    )
  }
  return(AtlasAnnDataToDataFrame(ad = ad, scorer = scorer))
}

# Query the Atlas, reusing cached results for parts of the region already
# fetched for the same scorer
#
# @inheritParams QueryAtlas
# @param cache Use the session cache
# @return A data.frame of variants in the region
#' @importFrom GenomicRanges GRanges reduce setdiff start end
#' @importFrom IRanges IRanges
#' @importFrom data.table rbindlist
QueryAtlasCached <- function(
  client,
  region,
  scorer = "AVI_SCORE_FEATURE_IMPORTANCE",
  cache = TRUE,
  ...
) {
  if (!cache) {
    return(QueryAtlas(client = client, region = region, scorer = scorer, ...))
  }
  if (length(x = region) != 1) {
    stop("region must describe a single genomic range")
  }
  chromosome <- as.character(x = seqnames(x = region))
  region <- GRanges(
    seqnames = chromosome,
    ranges = IRanges(start = start(x = region), end = end(x = region))
  )
  entry <- atlas.cache[[scorer]]
  if (is.null(x = entry)) {
    entry <- list(covered = GRanges(), data = NULL)
  }
  missing <- setdiff(x = region, y = entry$covered, ignore.strand = TRUE)
  if (length(x = missing) > 0) {
    fetched <- lapply(
      X = seq_len(length.out = length(x = missing)),
      FUN = function(i) {
        QueryAtlas(client = client, region = missing[i], scorer = scorer, ...)
      }
    )
    pieces <- c(list(entry$data), fetched)
    pieces <- pieces[!vapply(X = pieces, FUN = is.null, FUN.VALUE = logical(1))]
    entry$data <- as.data.frame(
      x = rbindlist(l = pieces, use.names = TRUE, fill = TRUE),
      stringsAsFactors = FALSE
    )
    entry$covered <- reduce(x = c(entry$covered, missing))
    atlas.cache[[scorer]] <- entry
  }
  keep <- entry$data[["chromosome"]] == chromosome &
    entry$data[["position"]] >= start(x = region) &
    entry$data[["position"]] <= end(x = region)
  avi <- entry$data[keep, , drop = FALSE]
  avi <- avi[order(avi[["position"]], avi[["alt"]]), , drop = FALSE]
  rownames(x = avi) <- NULL
  return(avi)
}

# Resolve the AVI SHAP features to plot
#
# Features can be requested by column name or by modality group label (see
# AVIFeatureGroups()), matched case-insensitively. Explicit column names may
# also refer to aggregate score columns such as AVI_SCORE. The keyword "Total"
# refers to the total AVI score, the sum of all feature attributions, which
# AVITrack() adds to the table.
#
# @param avi A data.frame with chromosome, position, ref and alt columns
# @param features Character vector of feature columns and/or group labels. If
# NULL, all feature attribution columns are used.
# @return A character vector of column names
ResolveAVIFeatures <- function(avi, features = NULL) {
  available <- AVIFeatureColumns(avi = avi, allow.empty = !is.null(x = features))
  if (is.null(x = features)) {
    return(available)
  }
  is.total <- tolower(x = features) %in% c("total", "total impact")
  features[is.total] <- "Total"
  groups <- AVIFeatureGroups(features = available)
  numeric.cols <- colnames(x = avi)[vapply(
    X = avi, FUN = is.numeric, FUN.VALUE = logical(length = 1)
  )]
  known.groups <- c(unique(x = AVIFeatureTable()[["group"]]), "Total")
  selected <- character()
  unknown <- character()
  absent <- character()
  for (f in features) {
    if (f %in% numeric.cols) {
      selected <- c(selected, f)
      next
    }
    hit <- available[tolower(x = groups) == tolower(x = f)]
    if (length(x = hit) == 0) {
      hit <- numeric.cols[tolower(x = numeric.cols) == tolower(x = f)]
    }
    if (length(x = hit) > 0) {
      selected <- c(selected, hit)
    } else if (tolower(x = f) %in% tolower(x = known.groups)) {
      # a valid group with no features in this data: skip it
      absent <- c(absent, f)
    } else {
      unknown <- c(unknown, f)
    }
  }
  if (length(x = unknown) > 0) {
    stop(
      "Requested features not found in AVI data: ",
      paste(unknown, collapse = ", "),
      ". Available features: ", paste(available, collapse = ", "),
      ". Available groups: ", paste(unique(x = groups), collapse = ", ")
    )
  }
  if (length(x = selected) == 0) {
    stop(
      "None of the requested feature groups (",
      paste(features, collapse = ", "), ") are present in the AVI data. ",
      "Available groups: ", paste(unique(x = groups), collapse = ", ")
    )
  }
  if (length(x = absent) > 0) {
    message(
      "No features found for group(s): ", paste(absent, collapse = ", ")
    )
  }
  return(unique(x = selected))
}

# Validate query limits before integer conversion or network access
ValidateAtlasQuery <- function(max.workers, max.retries) {
  values <- list(max.workers = max.workers, max.retries = max.retries)
  minimum <- c(max.workers = 1, max.retries = 0)
  for (name in names(x = values)) {
    value <- values[[name]]
    if (!is.numeric(x = value) || length(x = value) != 1 ||
        !is.finite(x = value) || value != floor(x = value) ||
        value < minimum[[name]]) {
      stop(name, " must be a finite integer of at least ", minimum[[name]])
    }
  }
}

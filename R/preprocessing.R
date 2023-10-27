#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

#' @param verbose Display messages
#' @rdname BinarizeCounts
#' @importFrom methods is slot "slot<-"
#' @export
#' @concept preprocessing
#' @examples
#' x <- matrix(data = sample(0:3, size = 25, replace = TRUE), ncol = 5)
#' BinarizeCounts(x)
BinarizeCounts.default <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  if (inherits(x = object, what = "CsparseMatrix")) {
    slot(object = object, name = "x") <- rep.int(
      x = 1,
      times = length(
        x = slot(object = object, name = "x")
        )
      )
  } else {
    object[object > 1] <- 1
  }
  return(object)
}

#' @rdname BinarizeCounts
#' @method BinarizeCounts Assay
#' @importFrom SeuratObject GetAssayData SetAssayData
#' @export
#' @concept preprocessing
#' @examples
#' BinarizeCounts(atac_small[['peaks']])
BinarizeCounts.Assay <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  data.matrix <- GetAssayData(object = object, slot = "counts")
  object <- SetAssayData(
    object = object,
    slot = "counts",
    new.data = BinarizeCounts(
      object = data.matrix, assay = assay, verbose = verbose
    )
  )
  return(object)
}

#' @param assay Name of assay to use. Can be a list of assays,
#' and binarization will be applied to each.
#' @rdname BinarizeCounts
#' @method BinarizeCounts Seurat
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept preprocessing
#' @examples
#' BinarizeCounts(atac_small)
BinarizeCounts.Seurat <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  for (i in seq_along(along.with = assay)) {
    assay.data <- object[[assay[[i]]]]
    assay.data <- BinarizeCounts(
      object = assay.data,
      assay = assay[[i]],
      verbose = verbose
    )
    object[[assay[[i]]]] <- assay.data
  }
  return(object)
}

#' Create motif matrix
#'
#' Create a motif x feature matrix from a set of genomic ranges,
#' the genome, and a set of position weight matrices.
#'
#' Requires that motifmatchr is installed
#' \url{https://www.bioconductor.org/packages/motifmatchr/}.
#'
#' @param features A GRanges object containing a set of genomic features
#' @param pwm A \code{\link[TFBSTools]{PFMatrixList}} or
#' \code{\link[TFBSTools]{PWMatrixList}}
#' object containing position weight/frequency matrices to use
#' @param genome Any object compatible with the \code{genome} argument
#' in \code{\link[motifmatchr]{matchMotifs}}
#' @param score Record the motif match score, rather than presence/absence
#' (default FALSE)
#' @param use.counts Record motif counts per region. If FALSE (default),
#' record presence/absence of motif. Only applicable if \code{score=FALSE}.
#' @param sep A length-2 character vector containing the separators to be used
#' when constructing matrix rownames from the GRanges
#' @param ... Additional arguments passed to
#' \code{\link[motifmatchr]{matchMotifs}}
#'
#' @return Returns a sparse matrix
#' @export
#' @concept motifs
#' @concept preprocessing
#' @examples
#' \dontrun{
#' library(JASPAR2018)
#' library(TFBSTools)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' pwm <- getMatrixSet(
#'   x = JASPAR2018,
#'   opts = list(species = 9606, all_versions = FALSE)
#' )
#' motif.matrix <- CreateMotifMatrix(
#'   features = granges(atac_small),
#'   pwm = pwm,
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
#' }
CreateMotifMatrix <- function(
  features,
  pwm,
  genome,
  score = FALSE,
  use.counts = FALSE,
  sep = c("-", "-"),
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
            immediate. = TRUE)
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
  rownames(motif.matrix) <- GRangesToString(grange = features, sep = sep)
  if (is.null(x = names(x = pwm))) {
    warning("No 'names' attribute found in PFMatrixList. ",
            "Extracting names from individual entries.")
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
    rownames(x = replacement_matrix) <- GRangesToString(
      grange = feature_order[miss_sn], sep = sep
    )
    colnames(x = replacement_matrix) <- colnames(x = motif.matrix)
    motif.matrix <- rbind(motif.matrix, replacement_matrix)
    motif.matrix <- motif.matrix[GRangesToString(
      grange = feature_order, sep = sep
      ), ]
  }
  return(motif.matrix)
}

#' Downsample Features
#'
#' Randomly downsample features and assign to VariableFeatures for the object.
#' This will select n features at random.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. Default is the active assay.
#' @param n Number of features to retain (default 20000).
#' @param verbose Display messages
#' @importFrom SeuratObject DefaultAssay GetAssayData "VariableFeatures<-"
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object with
#' \code{\link[SeuratObject]{VariableFeatures}} set to the randomly sampled features.
#' @export
#' @concept preprocessing
#' @examples
#' DownsampleFeatures(atac_small, n = 10)
DownsampleFeatures <- function(
  object,
  assay = NULL,
  n = 20000,
  verbose = TRUE
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (n > nrow(object[[assay]])) {
    stop("Requested more features than present in the assay")
  }
  if (verbose) {
    message("Randomly downsampling features")
  }
  VariableFeatures(object = object) <- sample(
    x = rownames(x = object[[assay]]), size = n, replace = FALSE
  )
  return(object)
}

#' @param assay Name of assay to use
#' @param min.cutoff Cutoff for feature to be included in the VariableFeatures
#' for the object. This can be a percentile specified as 'q' followed by the
#' minimum percentile, for example 'q5' to set the top 95\% most common features
#' as the VariableFeatures for the object. Alternatively, this can be an integer
#' specifying the minimum number of counts for the feature
#' to be included in the set of VariableFeatures. For example, setting to 10
#' will include features with >10 total counts in the set of VariableFeatures. If NULL,
#' include all features in VariableFeatures. If NA, VariableFeatures will not be
#' altered, and only the feature metadata will be updated with the total counts
#' and percentile rank for each feature.
#' @param verbose Display messages
#'
#' @importFrom Matrix rowSums
#' @importFrom stats ecdf
#' @rdname FindTopFeatures
#' @export
#' @concept preprocessing
#' @examples
#' FindTopFeatures(object = atac_small[['peaks']]['data'])
FindTopFeatures.default <- function(
  object,
  assay = NULL,
  min.cutoff = "q5",
  verbose = TRUE,
  ...
) {
  featurecounts <- rowSums(x = object)
  e.dist <- ecdf(x = featurecounts)
  hvf.info <- data.frame(
    row.names = names(x = featurecounts),
    count = featurecounts,
    percentile = e.dist(featurecounts)
  )
  hvf.info <- hvf.info[order(hvf.info$count, decreasing = TRUE), ]
  return(hvf.info)
}

#' @rdname FindTopFeatures
#' @importFrom SeuratObject GetAssayData VariableFeatures
#' @importFrom utils packageVersion
#' @export
#' @method FindTopFeatures Assay
#' @concept preprocessing
#' @examples
#' FindTopFeatures(object = atac_small[['peaks']])
FindTopFeatures.Assay <- function(
  object,
  assay = NULL,
  min.cutoff = "q5",
  verbose = TRUE,
  ...
) {
  data.use <- GetAssayData(object = object, slot = "counts")
  if (IsMatrixEmpty(x = data.use)) {
    if (verbose) {
      message("Count slot empty")
    }
    return(object)
  }
  hvf.info <- FindTopFeatures(
    object = data.use,
    assay = assay,
    min.cutoff = min.cutoff,
    verbose = verbose,
    ...
  )
  object[[names(x = hvf.info)]] <- hvf.info
  if (is.null(x = min.cutoff)) {
    VariableFeatures(object = object) <- rownames(x = hvf.info)
  } else if (is.numeric(x = min.cutoff)) {
    VariableFeatures(object = object) <- rownames(
      x = hvf.info[hvf.info$count > min.cutoff, ]
    )
  } else if (is.na(x = min.cutoff)) {
    # don't change the variable features
    return(object)
  } else {
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = min.cutoff))
    ) / 100
    VariableFeatures(object = object) <- rownames(
      x = hvf.info[hvf.info$percentile > percentile.use, ]
    )
  }
  return(object)
}

#' @rdname FindTopFeatures
#' @importFrom SeuratObject GetAssayData VariableFeatures
#' @importFrom utils packageVersion
#' @export
#' @method FindTopFeatures StdAssay
#' @concept preprocessing
#' @examples
#' FindTopFeatures(object = atac_small[['peaks']])
FindTopFeatures.StdAssay <- function(
    object,
    assay = NULL,
    min.cutoff = "q5",
    verbose = TRUE,
    ...
) {
  FindTopFeatures.Assay(
    object = object,
    assay = assay,
    min.cutoff = min.cutoff,
    verbose = verbose,
    ...
  )
}

#' @rdname FindTopFeatures
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept preprocessing
#' @method FindTopFeatures Seurat
#' @examples
#' FindTopFeatures(atac_small)
FindTopFeatures.Seurat <- function(
  object,
  assay = NULL,
  min.cutoff = "q5",
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object))
  assay.data <- object[[assay]]
  assay.data <- FindTopFeatures(
    object = assay.data,
    assay = assay,
    min.cutoff = min.cutoff,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' Calculate fraction of reads in peaks per cell
#'
#' @param object A Seurat object
#' @param assay Name of the assay containing a peak x cell matrix
#' @param total.fragments Name of a metadata column containing the total number
#' of sequenced fragments for each cell. This can be computed using the
#' \code{\link{CountFragments}} function.
#' @param col.name Name of column in metadata to store the FRiP information.
#' @param verbose Display messages
#'
#' @importFrom Matrix colSums
#' @importFrom SeuratObject GetAssayData AddMetaData
#'
#' @export
#' @concept qc
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @examples
#' FRiP(object = atac_small, assay = 'peaks', total.fragments = "fragments")
FRiP <- function(
  object,
  assay,
  total.fragments,
  col.name = "FRiP",
  verbose = TRUE
) {
  if (verbose) {
    message("Calculating fraction of reads in peaks per cell")
  }
  peak.data <- GetAssayData(object = object, assay = assay, slot = "counts")
  total_fragments_cell <- object[[]][[total.fragments]]
  peak.counts <- colSums(x = peak.data)
  frip <- peak.counts / total_fragments_cell
  object <- AddMetaData(object = object, metadata = frip, col.name = col.name)
  return(object)
}

globalVariables(names = "cell", package = "Signac")
#' NucleosomeSignal
#'
#' Calculate the strength of the nucleosome signal per cell.
#' Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to
#' fragments < 147 bp (nucleosome-free)
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. Only required if a fragment path is not
#' provided. If NULL, use the active assay.
#' @param n Number of lines to read from the fragment file. If NULL, read all
#' lines. Default scales with the number of cells in the object.
#' @param verbose Display messages
#' @param ... Arguments passed to other functions
#'
#' @importFrom dplyr group_by summarize
#' @importFrom stats ecdf
#'
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object with
#' added metadata for the ratio of mononucleosomal to nucleosome-free fragments
#' per cell, and the percentile rank of each ratio.
#' @export
#' @concept qc
#' @importFrom fastmatch fmatch
#' @importFrom SeuratObject AddMetaData
#' @importFrom stats ecdf
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   tolerance = 0.5
#' )
#' NucleosomeSignal(object = atac_small)
NucleosomeSignal <- function(
  object,
  assay = NULL,
  n = ncol(object) * 5e3,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  # first check that fragments are present
  frags <- Fragments(object = object[[assay]])
  if (length(x = frags) == 0) {
    stop("No fragment files present in assay")
  }
  verbose <- as.logical(x = verbose)
  af <- list()
  for (i in seq_along(along.with = frags)) {
    counts <- ExtractFragments(
      fragments = frags[[i]],
      n = n,
      verbose = verbose
    )
    cells.keep <- fmatch(
      x = counts$CB, table = colnames(x = object), nomatch = 0L
    )
    rownames(x = counts) <- counts$CB
    counts <- counts[
      cells.keep > 0, c("mononucleosomal", "nucleosome_free")
    ]
    af[[i]] <- counts
  }
  af <- do.call(what = rbind, args = af)
  af$nucleosome_signal <- af$mononucleosomal / af$nucleosome_free
  e.dist <- ecdf(x = af$nucleosome_signal)
  af$nucleosome_percentile <- round(
    x = e.dist(af$nucleosome_signal),
    digits = 2
  )
  af <- af[, c("nucleosome_signal", "nucleosome_percentile")]
  object <- AddMetaData(object = object, metadata = af)
  return(object)
}


#' @param genome A \code{BSgenome} object or any other object supported by
#' \code{getSeq}. Do \code{showMethods("getSeq")} to get the list of all
#' supported object types.
#' @param verbose Display messages
#'
#' @importMethodsFrom GenomicRanges width
#' @rdname RegionStats
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RegionStats(
#'   object = rownames(atac_small),
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
#' }
RegionStats.default <- function(
  object,
  genome,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace('BSgenome', quietly = TRUE)) {
    stop("Please install BSgenome: BiocManager::install('BSgenome')")
  }
  if (!requireNamespace('Biostrings', quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }
  sequence.length <- width(x = object)
  common.seq <- intersect(x = seqlevels(x = object), y = seqlevels(x = genome))
  if (length(x = common.seq) < length(x = seqlevels(x = object))) {
    warning("Not all seqlevels present in supplied genome", immediate. = TRUE)
  }
  seq.keep <- as.character(x = seqnames(x = object)) %in% common.seq
  enum <- seq_along(along.with = seq.keep)
  object <- object[seq.keep]
  object <- keepSeqlevels(x = object, value = common.seq, pruning.mode = "coarse")
  sequences <- Biostrings::getSeq(x = genome, object)
  gc <- Biostrings::letterFrequency(
    x = sequences, letters = 'CG'
  ) / sequence.length[seq.keep] * 100
  colnames(gc) <- 'GC.percent'
  dinuc <- Biostrings::dinucleotideFrequency(sequences)
  sequence.stats <- cbind(dinuc, gc)
  # fill missing seqnames with NA
  nadf <- as.data.frame(
    x = matrix(nrow = sum(!seq.keep), ncol = ncol(x = sequence.stats))
  )
  colnames(x = nadf) <- colnames(x = sequence.stats)
  rownames(x = nadf) <- enum[!seq.keep]
  rownames(x = sequence.stats) <- enum[seq.keep]
  sequence.stats <- rbind(sequence.stats, nadf)
  sequence.stats <- sequence.stats[enum, ]
  sequence.stats <- cbind(sequence.stats, sequence.length)
  return(sequence.stats)
}

#' @rdname RegionStats
#' @method RegionStats ChromatinAssay
#' @importFrom methods slot
#' @importFrom SeuratObject GetAssayData
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RegionStats(
#'   object = atac_small[['peaks']],
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
#' }
RegionStats.ChromatinAssay <- function(
  object,
  genome,
  verbose = TRUE,
  ...
) {
  regions <- granges(x = object)
  feature.metadata <- RegionStats(
    object = regions,
    genome = genome,
    verbose = verbose,
    ...
  )
  rownames(x = feature.metadata) <- rownames(x = object)
  meta.data <- GetAssayData(object = object, slot = "meta.features")
  feature.metadata <- feature.metadata[rownames(x = meta.data), ]
  meta.data <- cbind(meta.data, feature.metadata)
  slot(object = object, name = "meta.features") <- meta.data
  return(object)
}

#' @param assay Name of assay to use
#' @rdname RegionStats
#' @method RegionStats Seurat
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RegionStats(
#'   object = atac_small,
#'   assay = 'bins',
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
#' }
RegionStats.Seurat <- function(
  object,
  genome,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  assay.data <- object[[assay]]
  assay.data <- RegionStats(
    object = assay.data,
    genome = genome,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' @param method Which TF-IDF implementation to use. Choice of:
#' \itemize{
#'  \item{1}: The TF-IDF implementation used by Stuart & Butler et al. 2019
#'  (\doi{10.1101/460147}). This computes
#'  \eqn{\log(TF \times IDF)}.
#'  \item{2}: The TF-IDF implementation used by Cusanovich & Hill
#'  et al. 2018 (\doi{10.1016/j.cell.2018.06.052}). This
#'  computes \eqn{TF \times (\log(IDF))}.
#'  \item{3}: The log-TF method used by Andrew Hill.
#'  This computes \eqn{\log(TF) \times \log(IDF)}.
#'  \item{4}: The 10x Genomics method (no TF normalization). This computes
#'  \eqn{IDF}.
#' }
#' @param scale.factor Which scale factor to use. Default is 10000.
#' @param idf A precomputed IDF vector to use. If NULL, compute based on the
#' input data matrix.
#' @param verbose Print progress
#' @rdname RunTFIDF
#' @importFrom Matrix colSums rowSums Diagonal tcrossprod
#' @importFrom methods is "slot<-" slot
#' @export
#' @concept preprocessing
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' RunTFIDF(object = mat)
RunTFIDF.default <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  idf = NULL,
  verbose = TRUE,
  ...
) {
  if (inherits(x = object, what = "data.frame")) {
    object <- as.matrix(x = object)
  }
  if (!inherits(x = object, what = "CsparseMatrix")) {
    object <- as(object = object, Class = "CsparseMatrix")
  }
  if (verbose) {
    message("Performing TF-IDF normalization")
  }
  npeaks <- colSums(x = object)
  if (any(npeaks == 0)) {
    warning("Some cells contain 0 total counts")
  }
  if (method == 4) {
    tf <- object
  } else {
    tf <- tcrossprod(x = object, y = Diagonal(x = 1 / npeaks))
  }
  if (!is.null(x = idf)) {
    precomputed_idf <- TRUE
    if (!inherits(x = idf, what = "numeric")) {
      stop("idf parameter must be a numeric vector")
    }
    if (length(x = idf) != nrow(x = object)) {
      stop("Length of supplied IDF vector does not match",
           " number of rows in input matrix")
    }
    if (any(idf == 0)) {
      stop("Supplied IDF values cannot be zero")
    }
    if (verbose) {
      message("Using precomputed IDF vector")
    }
  } else {
    precomputed_idf <- FALSE
    rsums <- rowSums(x = object)
    if (any(rsums == 0)) {
      warning("Some features contain 0 total counts")
    }
    idf <- ncol(x = object) / rsums
  }

  if (method == 2) {
    if (!precomputed_idf) {
      idf <- log(1 + idf)
    }
  } else if (method == 3) {
    slot(object = tf, name = "x") <- log1p(
      x = slot(object = tf, name = "x") * scale.factor
    )
    if (!precomputed_idf) {
      idf <- log(1 + idf)
    }
  }
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  if (method == 1) {
    slot(object = norm.data, name = "x") <- log1p(
      x = slot(object = norm.data, name = "x") * scale.factor
    )
  }
  colnames(x = norm.data) <- colnames(x = object)
  rownames(x = norm.data) <- rownames(x = object)
  # set NA values to 0
  vals <- slot(object = norm.data, name = "x")
  vals[is.na(x = vals)] <- 0
  slot(object = norm.data, name = "x") <- vals
  return(norm.data)
}

#' @rdname RunTFIDF
#' @method RunTFIDF Assay
#' @export
#' @concept preprocessing
#' @examples
#' RunTFIDF(atac_small[['peaks']])
RunTFIDF.Assay <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  idf = NULL,
  verbose = TRUE,
  ...
) {
  new.data <- RunTFIDF(
    object = GetAssayData(object = object, slot = "counts"),
    method = method,
    assay = assay,
    scale.factor = scale.factor,
    idf = idf,
    verbose = verbose,
    ...
  )
  new.data <- as(object = new.data, Class = "CsparseMatrix")
  object <- SetAssayData(
    object = object,
    slot = "data",
    new.data = new.data
  )
  return(object)
}

#' @rdname RunTFIDF
#' @method RunTFIDF StdAssay
#' @export
#' @concept preprocessing
#' @examples
#' RunTFIDF(atac_small[['peaks']])
RunTFIDF.StdAssay <- function(
    object,
    assay = NULL,
    method = 1,
    scale.factor = 1e4,
    idf = NULL,
    verbose = TRUE,
    ...
) {
  RunTFIDF.Assay(
    object = object,
    assay = assay,
    method = method,
    scale.factor = scale.factor,
    idf = idf,
    verbose = verbose,
    ...
  )
}

#' @param assay Name of assay to use
#' @rdname RunTFIDF
#' @method RunTFIDF Seurat
#' @export
#' @concept preprocessing
#' @examples
#' RunTFIDF(object = atac_small)
RunTFIDF.Seurat <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  idf = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object))
  assay.data <- object[[assay]]
  assay.data <- RunTFIDF(
    object = assay.data,
    assay = assay,
    method = method,
    scale.factor = scale.factor,
    idf = idf,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

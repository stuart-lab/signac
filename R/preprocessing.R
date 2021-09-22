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
  if (inherits(x = object, what = "dgCMatrix")) {
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
#' @importFrom Seurat GetAssayData SetAssayData
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
#' @importFrom Seurat GetAssay DefaultAssay
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
    assay.data <- GetAssay(object = object, assay = assay[[i]])
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
      motif.matrix <- as(Class = "dgCMatrix", object = motif.matrix)
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
#' @importFrom Seurat DefaultAssay GetAssayData "VariableFeatures<-"
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
#' specifying the minimum number of cells containing the feature for the feature
#' to be included in the set of VariableFeatures. For example, setting to 10
#' will include features in >10 cells in the set of VariableFeatures. If NULL,
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
#' FindTopFeatures(object = atac_small[['peaks']][])
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
#' @importFrom Seurat GetAssayData VariableFeatures
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
#' @importFrom Seurat DefaultAssay GetAssay
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
  assay.data <- GetAssay(object = object, assay = assay)
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
#' @importFrom Seurat GetAssayData AddMetaData
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
#' @importFrom Seurat AddMetaData
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


#' @param genome A BSgenome object
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
  sequences <- BSgenome::getSeq(x = genome, names = object)
  gc <- Biostrings::letterFrequency(
    x = sequences, letters = 'CG'
  ) / sequence.length * 100
  colnames(gc) <- 'GC.percent'
  dinuc <- Biostrings::dinucleotideFrequency(sequences)
  sequence.stats <- cbind(dinuc, gc, sequence.length)
  return(sequence.stats)
}

#' @rdname RegionStats
#' @method RegionStats ChromatinAssay
#' @importFrom methods slot
#' @importFrom Seurat GetAssayData
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
  assay.data <- GetAssay(object = object, assay = assay)
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
  verbose = TRUE,
  ...
) {
  if (inherits(x = object, what = "data.frame")) {
    object <- as.matrix(x = object)
  }
  if (!inherits(x = object, what = "dgCMatrix")) {
    object <- as(object = object, Class = "dgCMatrix")
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
  rsums <- rowSums(x = object)
  if (any(rsums == 0)) {
    warning("Some features contain 0 total counts")
  }
  idf <- ncol(x = object) / rsums
  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    slot(object = tf, name = "x") <- log1p(
      x = slot(object = tf, name = "x") * scale.factor
    )
    idf <- log(1 + idf)
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
  verbose = TRUE,
  ...
) {
  new.data <- RunTFIDF(
    object = GetAssayData(object = object, slot = "counts"),
    method = method,
    assay = assay,
    scale.factor = scale.factor,
    verbose = verbose,
    ...
  )
  new.data <- as(object = new.data, Class = "dgCMatrix")
  object <- SetAssayData(
    object = object,
    slot = "data",
    new.data = new.data
  )
  return(object)
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
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object))
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- RunTFIDF(
    object = assay.data,
    assay = assay,
    method = method,
    scale.factor = scale.factor,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' Compute TSS enrichment score per cell
#'
#' Compute the transcription start site (TSS) enrichment score for each cell,
#' as defined by ENCODE:
#' \url{https://www.encodeproject.org/data-standards/terms/}.
#'
#' The computed score will be added to the object metadata as "TSS.enrichment".
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param tss.positions A GRanges object containing the TSS positions. If NULL,
#' use the genomic annotations stored in the assay.
#' @param n Number of TSS positions to use. This will select the first _n_
#' TSSs from the set. If NULL, use all TSSs (slower).
#' @param cells A vector of cells to include. If NULL (default), use all cells
#' in the object
#' @param fast Just compute the TSS enrichment score, without storing the
#' base-resolution matrix of integration counts at each site. This reduces the
#' memory required to store the object but does not allow plotting the
#' accessibility profile at the TSS.
#' @param process_n Number of regions to process at a time if using \code{fast}
#' option.
#' @param verbose Display messages
#'
#' @importFrom Matrix rowMeans
#' @importFrom methods slot
#' @importFrom stats ecdf
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges start width strand
#' @importFrom Seurat DefaultAssay
#'
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @export
#' @concept qc
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   tolerance = 0.5
#' )
#' TSSEnrichment(object = atac_small)
#' }
TSSEnrichment <- function(
  object,
  tss.positions = NULL,
  n = NULL,
  fast = TRUE,
  assay = NULL,
  cells = NULL,
  process_n = 2000,
  verbose = TRUE
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
  if (is.null(x = tss.positions)) {
    if (verbose) {
      message("Extracting TSS positions")
    }
    # work out TSS positions from gene annotations
    annotations <- Annotation(object = object[[assay]])
    if (is.null(x = annotations)) {
      stop("No gene annotations present in assay")
    }
    tss.positions <- GetTSSPositions(ranges = annotations)
  }
  if (!is.null(x = n)) {
    if (n > length(x = tss.positions)) {
      n <- length(x = tss.positions)
    }
    tss.positions <- tss.positions[1:n, ]
  }

  # exclude chrM
  sn <- seqnames(x = tss.positions)
  tss.positions <- tss.positions[!as.character(sn) %in% c("chrM", "Mt", "MT")]

  if (fast) {
    # just compute the TSS enrichment score without storing the full matrix
    object <- TSSFast(
      object = object,
      assay = assay,
      tss.positions = tss.positions,
      process_n = process_n,
      verbose = verbose
    )
    return(object)
  }

  tss.positions <- Extend(
    x = tss.positions,
    upstream = 1000,
    downstream = 1000,
    from.midpoint = TRUE
  )
  cutmatrix <- CreateRegionPileupMatrix(
    object = object,
    regions = tss.positions,
    assay = assay,
    cells = cells,
    verbose = verbose
  )

  # compute mean read counts in 100 bp at each flank for each cell
  # (200 bp total averaged)
  if (verbose) {
    message("Computing mean insertion frequency in flanking regions")
  }
  flanking.mean <- rowMeans(x = cutmatrix[, c(1:100, 1902:2001)])

  # if the flanking mean is 0 for any cells, the enrichment score will be zero.
  # instead replace with the mean from the whole population
  flanking.mean[is.na(x = flanking.mean)] <- 0
  flanking.mean[flanking.mean == 0] <- mean(flanking.mean, na.rm = TRUE)

  # compute fold change at each position relative to flanking mean
  # (flanks should start at 1)
  if (verbose) {
    message("Normalizing TSS score")
  }

  norm.matrix <- cutmatrix / flanking.mean

  # Take signal value at center of distribution after normalization as
  # TSS enrichment score, average the 1000 bases at the center
  object$TSS.enrichment <- rowMeans(x = norm.matrix[, 500:1500], na.rm = TRUE)
  e.dist <- ecdf(x = object$TSS.enrichment)
  object$TSS.percentile <- round(
    x = e.dist(object$TSS.enrichment),
    digits = 2
  )

  # store expected as one additional row in the matrix
  expected.insertions <- rep(1, ncol(x = cutmatrix))
  expected.insertions <- t(x = as.matrix(x = expected.insertions))
  rownames(x = expected.insertions) <- "expected"

  # encode motif position as additional row in matrix
  motif.vec <- t(x = matrix(
    data = c(
      rep(x = 0, 1000),
      1,
      rep(x = 0, 1000)
    )
  )
  )
  rownames(x = motif.vec) <- "motif"

  # append
  norm.matrix <- rbind(norm.matrix, expected.insertions)
  norm.matrix <- rbind(norm.matrix, motif.vec)

  # store the normalized TSS matrix
  object <- suppressWarnings(SetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment",
    new.data = norm.matrix,
    key = "TSS"
  ))
  return(object)
}

#' @importFrom Rsamtools TabixFile
#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @importFrom stats ecdf
#' @importFrom Matrix rowSums
#' @importFrom Seurat DefaultAssay
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pblapply
TSSFast <- function(
  object,
  tss.positions,
  assay = NULL,
  process_n = 2000,
  verbose = TRUE
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))

  # extract fragments
  frags <- Fragments(object = object[[assay]])
  if (length(x = frags) == 0) {
    stop("Fragments file not set for assay ", assay)
  }

  # get regions
  upstream.flank <- Extend(
    x = tss.positions,
    upstream = 1000,
    downstream = -901,
    from.midpoint = TRUE
  )
  downstream.flank <- Extend(
    x = tss.positions,
    upstream = -901,
    downstream = 1000,
    from.midpoint = TRUE
  )
  centers <- Extend(
    x = tss.positions,
    upstream = 500,
    downstream = 500,
    from.midpoint = TRUE
  )

  # chunk ranges
  process_n <- SetIfNull(x = process_n, y = length(x = centers))
  nchunk <- ceiling(x = length(x = upstream.flank) / process_n)
  upstream.flank <- ChunkGRanges(
    granges = upstream.flank,
    nchunk = nchunk
  )
  downstream.flank <- ChunkGRanges(
    granges = downstream.flank,
    nchunk = nchunk
  )
  centers <- ChunkGRanges(
    granges = centers,
    nchunk = nchunk
  )

  # iterate over fragment files and parts of region
  if (verbose) {
    message("Extracting fragments at TSSs")
    pb <- txtProgressBar(
      min = 1,
      max = length(x = centers),
      style = 3,
      file = stderr()
    )
  }

  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }

  center.counts <- vector(mode = "numeric", length = ncol(x = object))
  flank.counts <- vector(mode = "numeric", length = ncol(x = object))
  for (i in seq_along(along.with = frags)) {
    # open fragment file
    tbx.path <- GetFragmentData(object = frags[[i]], slot = "path")
    cellmap <- GetFragmentData(object = frags[[i]], slot = "cells")
    if (is.null(x = cellmap)) {
      cellmap <- colnames(x = object)
      names(x = cellmap) <- cellmap
    } else {
      cellmap <- cellmap[intersect(names(x = cellmap), colnames(x = object))]
    }
    tbx <- TabixFile(file = tbx.path)
    # iterate over chunked ranges
    res <- mylapply(
      X = seq_along(along.with = centers),
      FUN = function(x) {
        extract_tss_counts(
          cellnames = colnames(x = object),
          region.centers = centers[[x]],
          upstream = upstream.flank[[x]],
          downstream = downstream.flank[[x]],
          tabix.file = tbx,
          cell.name.map = cellmap
        )
      }
    )

    # sum results from each chunk of granges
    cc <- lapply(X = res, FUN = `[[`, 1)
    fc <- lapply(X = res, FUN = `[[`, 2)
    center.counts <- center.counts + Reduce(f = `+`, x = cc)
    flank.counts <- flank.counts + Reduce(f = `+`, x = fc)
  }

  if (verbose) {
    message("\nComputing TSS enrichment score")
  }

  # take mean accessibility per base
  flank.mean <- flank.counts / 200
  flank.mean[flank.counts == 0] <- mean(x = flank.mean, na.rm = TRUE)

  center.norm <- center.counts / flank.mean

  # compute TSS enrichment score and add to object
  object$TSS.enrichment <- center.norm / 1001
  e.dist <- ecdf(x = object$TSS.enrichment)
  object$TSS.percentile <- round(
    x = e.dist(object$TSS.enrichment),
    digits = 2
  )
  return(object)
}

#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @importFrom Rsamtools seqnamesTabix
#' @importFrom Matrix rowSums
extract_tss_counts <- function(
  cellnames,
  region.centers,
  tabix.file,
  upstream,
  downstream,
  cell.name.map
) {
  tabix.file <- open(con = tabix.file)
  # initialize vectors
  fc <- vector(mode = "numeric", length = length(x = cellnames))
  names(x = fc) <- cellnames
  cc <- vector(mode = "numeric", length = length(x = cellnames))
  names(x = cc) <- cellnames

  # remove seqlevels not present in fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region.centers),
    y = seqnamesTabix(file = tabix.file)
  )
  if (length(x = common.seqlevels) == 0) {
    close(con = tabix.file)
    return(list(cc, fc))
  }
  uflanks.use <- keepSeqlevels(
    x = upstream,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  dflanks.use <- keepSeqlevels(
    x = downstream,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  centers.use <- keepSeqlevels(
    x = region.centers,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )

  # count integration events
  cuts.center <- SingleFileCutMatrix(
    cellmap = cell.name.map,
    tabix.file = tabix.file,
    region = centers.use,
    verbose = FALSE
  )
  counts.center <- rowSums(x = cuts.center)
  cuts.flank <- SingleFileCutMatrix(
    cellmap = cell.name.map,
    tabix.file = tabix.file,
    region = c(uflanks.use, dflanks.use),
    verbose = FALSE
  )
  counts.flank <- rowSums(x = cuts.flank)
  cc[names(x = counts.center)] <- cc + as.vector(x = counts.center)
  fc[names(x = counts.flank)] <- fc + as.vector(x = counts.flank)
  close(con = tabix.file)
  return(list(cc, fc))
}

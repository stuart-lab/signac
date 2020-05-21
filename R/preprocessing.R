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
#' @return Returns a \code{\link[Seurat]{Seurat}} object with
#' \code{\link[Seurat]{VariableFeatures}} set to the randomly sampled features.
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

#' Feature Matrix
#'
#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments A list of \code{\link{Fragment}} objects.
#' @param features A GRanges object containing a set of genomic intervals.
#' These will form the rows of the matrix, with each entry recording the number
#' of unique reads falling in the genomic region for each cell.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param chunk Number of chunks to use when processing the fragments file.
#' Fewer chunks may enable faster processing,
#'  but will use more memory.
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param verbose Display messages
#'
#' @export
#' @concept preprocessing
#' @concept utilities
#' @return Returns a sparse matrix
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(fpath)
#' FeatureMatrix(
#'   fragments = fragments,
#'   features = granges(atac_small)
#' )
FeatureMatrix <- function(
  fragments,
  features,
  cells = NULL,
  chunk = 50,
  sep = c("-", "-"),
  verbose = TRUE
) {
  if (!inherits(x = fragments, what = "list")) {
    if (inherits(x = fragments, what = "Fragment")) {
      fragments <- list(fragments)
    } else {
      stop("fragments should be a list of Fragment objects")
    }
  }
  # if cells is not NULL, iterate over all fragment objects
  # and find which objects contain cells that are requested
  if (!is.null(x = cells)) {
    obj.use <- c()
    for (i in seq_along(along.with = fragments)) {
      if (any(cells %in% Cells(x = fragments[[i]]))) {
        obj.use <- c(obj.use, i)
      }
    }
  } else {
    obj.use <- seq_along(along.with = fragments)
  }
  # create a matrix from each fragment file
  mat.list <- sapply(
    X = obj.use,
    FUN = function(x) {
      SingleFeatureMatrix(
        fragment = fragments[[x]],
        features = features,
        cells = cells,
        sep = sep,
        verbose = verbose,
        chunk = chunk
      )
    })
  # cbind all the matrices
  featmat <- do.call(what = cbind, args = mat.list)
  return(featmat)
}

# Run FeatureMatrix on a single Fragment object
# @inheritParams FeatureMatrix
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom Matrix sparseMatrix
#' @importMethodsFrom GenomicRanges intersect
#' @importFrom Rsamtools TabixFile seqnamesTabix
SingleFeatureMatrix <- function(
  fragment,
  features,
  cells = NULL,
  chunk = 50,
  sep = c("-", "-"),
  verbose = TRUE
) {
  fragment.path <- GetFragmentData(object = fragment, slot = "path")
  if (!is.null(cells)) {
    # only look for cells that are in the fragment file
    cells <- intersect(x = cells, y = Cells(x = fragment))
  }
  tbx <- TabixFile(file = fragment.path)
  features <- keepSeqlevels(
    x = features,
    value = intersect(
      x = seqnames(x = features),
      y = seqnamesTabix(file = tbx)
    ),
    pruning.mode = "coarse"
  )
  feature.list <- ChunkGRanges(
    granges = features,
    nchunk = chunk
  )
  if (verbose) {
    message("Extracting reads overlapping genomic regions")
  }
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  cells.in.regions <- mylapply(
    X = feature.list,
    FUN = GetCellsInRegion,
    tabix = tbx,
    cells = cells,
    sep = sep
  )
  if (verbose) {
    message("Constructing matrix")
  }
  cell.vector <- unlist(x = lapply(X = cells.in.regions, FUN = `[[`, 1))
  feature.vector <- unlist(x = lapply(X = cells.in.regions, FUN = `[[`, 2))
  all.cells <- unique(x = cell.vector)
  all.features <- unique(x = feature.vector)
  cell.lookup <- seq_along(along.with = all.cells)
  feature.lookup <- seq_along(along.with = all.features)
  names(x = cell.lookup) <- all.cells
  names(x = feature.lookup) <- all.features
  matrix.features <- feature.lookup[feature.vector]
  matrix.cells <- cell.lookup[cell.vector]
  featmat <- sparseMatrix(
    i = matrix.features,
    j = matrix.cells,
    x = rep(x = 1, length(x = cell.vector))
  )
  featmat <- as(Class = "dgCMatrix", object = featmat)
  rownames(x = featmat) <- names(x = feature.lookup)
  colnames(x = featmat) <- names(x = cell.lookup)
  # add zero columns for missing cells
  if (!is.null(x = cells)) {
    missing.cells <- setdiff(x = cells, y = colnames(x = featmat))
    if (!(length(x = missing.cells) == 0)) {
      null.mat <- sparseMatrix(
        i = c(),
        j = c(),
        dims = c(nrow(x = featmat), length(missing.cells))
      )
      rownames(x = null.mat) <- rownames(x = featmat)
      colnames(x = null.mat) <- missing.cells
      featmat <- cbind(featmat, null.mat)
    }
    featmat <- featmat[, cells]
  }
  # add zero rows for missing features
  all.features <- GRangesToString(grange = features, sep = sep)
  missing.features <- all.features[!(all.features %in% rownames(featmat))]
  if (length(x = missing.features) > 0) {
    null.mat <- sparseMatrix(
      i = c(),
      j = c(),
      dims = c(length(x = missing.features), ncol(x = featmat))
    )
    rownames(x = null.mat) <- missing.features
    featmat <- rbind(featmat, null.mat)
  }
  return(featmat)
}


#' @param assay Name of assay to use
#' @param min.cutoff Cutoff for feature to be included in the VariableFeatures
#' for the object. This can be a percentile specified as 'q' followed by the
#' minimum percentile, for example 'q5' to set the top 95\% most common features
#' as the VariableFeatures for the object. Alternatively, this can be an integer
#' specifying the minumum number of cells containing the feature for the feature
#' to be included in the set of VariableFeatures. For example, setting to 10
#' will include features in >10 cells in the set of VariableFeatures. If NULL,
#' include all features in VariableFeatures.
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
#' @param peak.assay Name of the assay containing a peak x cell matrix
#' @param bin.assay Name of the assay containing a bin x cell matrix
#' @param chromosome Which chromosome to use. Default is chromosome 1 ('chr1').
#' If NULL, use the whole genome.
#' @param verbose Display messages
#'
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData AddMetaData
#'
#' @export
#' @concept qc
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @examples
#' FRiP(object = atac_small, peak.assay = 'peaks', bin.assay = 'bins')
FRiP <- function(
  object,
  peak.assay,
  bin.assay,
  chromosome = "chr1",
  verbose = TRUE
) {
  if (verbose) {
    message("Calculating fraction of reads in peaks per cell")
  }
  peak.data <- GetAssayData(
    object = object, assay = peak.assay, slot = "counts"
  )
  bin.data <- GetAssayData(
    object = object, assay = bin.assay, slot = "counts"
  )
  if (!is.null(x = chromosome)) {
    peak.data <- peak.data[grepl(
      pattern = paste0("^", chromosome, "\\-|^", chromosome, ":"),
      x = rownames(x = peak.data)
    ), ]
    bin.data <- bin.data[grepl(
      pattern = paste0("^", chromosome, "\\-|^", chromosome, ":"),
      x = rownames(x = bin.data)
    ), ]
  }
  peak.counts <- colSums(x = peak.data)
  bin.counts <- colSums(x = bin.data)
  frip <- peak.counts / bin.counts
  object <- AddMetaData(object = object, metadata = frip, col.name = "FRiP")
  return(object)
}

#' Genome bin matrix
#'
#' Construct a bin x cell matrix from a fragments file.
#'
#' This function bins the genome and calls \code{\link{FeatureMatrix}} to
#' construct a bin x cell matrix.
#'
#' @param fragments Path to tabix-indexed fragments file
#' @param genome A vector of chromosome sizes for the genome. This is used to
#' construct the genome bin coordinates. The can be obtained by calling
#' \code{\link[GenomeInfoDb]{seqlengths}} on a
#' \code{\link[BSgenome]{BSgenome-class}} object.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param binsize Size of the genome bins to use
#' @param chunk Number of chunks to use when processing the fragments file.
#' Fewer chunks may enable faster processing, but will use more memory.
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param verbose Display messages
#'
#' @importFrom GenomicRanges tileGenome
#' @export
#' @concept preprocessing
#' @concept utilities
#' @return Returns a sparse matrix
#' @examples
#' genome <- 780007
#' names(genome) <- 'chr1'
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(fpath)
#' GenomeBinMatrix(
#'   fragments = fragments,
#'   genome = genome,
#'   binsize = 1000
#' )
GenomeBinMatrix <- function(
  fragments,
  genome,
  cells = NULL,
  binsize = 5000,
  chunk = 50,
  sep = c("-", "-"),
  verbose = TRUE
) {
  tiles <- tileGenome(
    seqlengths = genome,
    tilewidth = binsize,
    cut.last.tile.in.chrom = TRUE
  )
  binmat <- FeatureMatrix(
    fragments = fragments,
    features = tiles,
    cells = cells,
    chunk = chunk,
    sep = sep,
    verbose = verbose
  )
  return(binmat)
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
#' @param region Which region to use. Can be a GRanges region, a string, or a
#' vector of strings. Default is human chromosome 1.
#' @param min.threshold Lower bound for the mononucleosome size. Default is 147
#' @param max.threshold Upper bound for the mononucleosome size. Default is 294
#' @param verbose Display messages
#' @param ... Arguments passed to other functions
#'
#' @importFrom dplyr group_by summarize
#' @importFrom stats ecdf
#'
#' @return Returns a \code{\link[Seurat]{Seurat}} object with
#' added metadata for the ratio of mononucleosomal to nucleosome-free fragments
#' per cell, and the percentile rank of each ratio.
#' @export
#' @concept qc
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small)
#' )
#' NucleosomeSignal(object = atac_small)
NucleosomeSignal <- function(
  object,
  assay = NULL,
  region = "chr1-1-249250621",
  min.threshold = 147,
  max.threshold = 294,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  fragments.use <- MultiGetReadsInRegion(
    object = object,
    region = region,
    assay = assay,
    cells = colnames(x = object),
    verbose = verbose,
    ...
  )
  mn_ratio <- function(x) {
    mononucleosome <- sum(x[x > min.threshold & x < max.threshold])
    nucleosome_free <- sum(x[x <= min.threshold])
    return(mononucleosome / nucleosome_free)
  }
  if (verbose) {
    message("Computing ratio of mononucleosomal to nucleosome-free fragments")
  }
  fragments.use <- as.data.frame(x = fragments.use[, c("cell", "length")])
  fragments.use <- group_by(.data = fragments.use, cell)
  fragment.summary <- as.data.frame(
    x = summarize(
      .data = fragments.use,
      nucleosome_signal = mn_ratio(x = length)
    )
  )
  rownames(x = fragment.summary) <- fragment.summary$cell
  fragment.summary$cell <- NULL
  e.dist <- ecdf(x = fragment.summary$nucleosome_signal)
  fragment.summary$nucleosome_percentile <- round(
    x = e.dist(fragment.summary$nucleosome_signal),
    digits = 2
  )
  object <- AddMetaData(object = object, metadata = fragment.summary)
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
#'  (\url{https://doi.org/10.1101/460147}). This computes
#'  \eqn{\log(TF \times IDF)}.
#'  \item{2}: The TF-IDF implementation used by Cusanovich & Hill
#'  et al. 2018 (\url{https://doi.org/10.1016/j.cell.2018.06.052}). This
#'  computes \eqn{TF \times (\log(IDF))}.
#'  \item{3}: The log-TF method used by Andrew Hill (\url{http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/}).
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
  vals <- slot(object = object, name = "x")
  vals[is.na(x = vals)] <- 0
  slot(object = object, name = "x") <- vals
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
#' @param verbose Display messages
#' @importFrom Matrix rowMeans
#' @importFrom methods slot
#' @importFrom stats ecdf
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges start width strand
#'
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @export
#' @concept qc
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small)
#' )
#' TSSEnrichment(object = atac_small)
#' }
TSSEnrichment <- function(
  object,
  tss.positions = NULL,
  n = NULL,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (is.null(x = tss.positions)) {
    # work out TSS positions from gene annotations
    annotations <- Annotation(object = object[[assay]])
    tss.positions <- GetTSSPositions(ranges = annotations)
  }
  if (!is.null(x = n)) {
    if (n > length(x = tss.positions)) {
      n <- length(x = tss.positions)
    }
    tss.positions <- tss.positions[1:n, ]
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

  # compute mean read counts in 100 bp at eack flank for each cell
  # (200 bp total averaged)
  if (verbose) {
    message("Computing mean insertion frequency in flanking regions")
  }
  flanking.mean <- rowMeans(x = cutmatrix[, c(1:100, 1901:2001)])

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
  object$TSS.enrichment <- rowMeans(x = norm.matrix[, 501:1500], na.rm = TRUE)
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

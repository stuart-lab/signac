#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

#' @param verbose Display messages
#' @rdname BinarizeCounts
#' @importFrom methods is slot "slot<-"
#' @export
#' @examples
#' x <- matrix(data = sample(0:3, size = 25, replace = TRUE), ncol = 5)
#' BinarizeCounts(x)
BinarizeCounts.default <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  if (inherits(x = object, what = 'dgCMatrix')) {
    slot(object = object, name = 'x') <- rep.int(
      x = 1,
      times = length(
        x = slot(object = object, name = 'x')
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
#' @examples
#' BinarizeCounts(atac_small[['peaks']])
BinarizeCounts.Assay <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  data.matrix <- GetAssayData(object = object, slot = 'counts')
  object <- SetAssayData(
    object = object,
    slot = 'counts',
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

#' CreateMotifMatrix
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
#'   features = StringToGRanges(rownames(atac_small), sep = c(":", "-")),
#'   pwm = pwm,
#'   genome = BSgenome.Hsapiens.UCSC.hg19,
#'   sep = c(":", "-")
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
  if (!requireNamespace('motifmatchr', quietly = TRUE)) {
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
      motif.matrix <- as(Class = 'dgCMatrix', object = motif.matrix)
    }
  }
  rownames(motif.matrix) <- GRangesToString(grange = features, sep = sep)
  return(motif.matrix)
}

#' DownsampleFeatures
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

#' FeatureMatrix
#'
#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments Path to tabix-indexed fragments file
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
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom Matrix sparseMatrix
#' @importMethodsFrom GenomicRanges intersect
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @export
#' @return Returns a sparse matrix
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' FeatureMatrix(
#'   fragments = fpath,
#'   features = StringToGRanges(rownames(atac_small), sep = c(":", "-"))
#' )
FeatureMatrix <- function(
  fragments,
  features,
  cells = NULL,
  chunk = 50,
  sep = c('-', '-'),
  verbose = TRUE
) {
  tbx <- TabixFile(file = fragments)
  features <- keepSeqlevels(
    x = features,
    value = intersect(
      x = seqnames(x = features), y = seqnamesTabix(file = tbx)
    ),
    pruning.mode = "coarse"
  )
  feature.list <- ChunkGRanges(
    granges = features,
    nchunk = chunk
  )
  if (verbose) {
    message('Extracting reads overlapping genomic regions')
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
  featmat <- as(Class = 'dgCMatrix', object = featmat)
  rownames(x = featmat) <- names(x = feature.lookup)
  colnames(x = featmat) <- names(x = cell.lookup)
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
    return(featmat[, cells])
  } else {
    return(featmat)
  }
}

globalVariables(names = c('chr', 'start'), package = 'Signac')
#' FilterFragments
#'
#' Remove cells from a fragments file that are not present in a given list of
#' cells. Note that this reads the whole fragments file into memory, so may
#' require a lot of memory depending on the size of the fragments file.
#'
#' @param fragment.path Path to a tabix-indexed fragments file
#' @param cells A vector of cells to retain
#' @param output.path Name and path for output tabix file. A tabix index file
#' will also be created in the same location, with the .tbi file extension.
#' @param assume.sorted Assume sorted input and don't sort the filtered file.
#' Can save a lot of time, but indexing will fail if assumption is wrong.
#' @param compress Compress filtered fragments using bgzip (default TRUE)
#' @param index Index the filtered tabix file (default TRUE)
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[data.table]{fread}}
#'
#' @importFrom data.table fread fwrite
#' @importFrom Rsamtools indexTabix bgzip
#' @export
#' @return None
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' output.path = file.path(tempdir(), "filtered.tsv")
#'
#' FilterFragments(
#'   fragment.path = fpath,
#'   cells = colnames(atac_small),
#'   output.path = output.path
#' )
#' }
FilterFragments <- function(
  fragment.path,
  cells,
  output.path,
  assume.sorted = FALSE,
  compress = TRUE,
  index = TRUE,
  verbose = TRUE,
  ...
) {
  if (verbose) {
    message("Retaining ", length(x = cells), " cells")
    message("Reading fragments")
  }
  reads <- fread(
    file = fragment.path,
    col.names = c('chr', 'start', 'end', 'cell', 'count'),
    showProgress = verbose,
    ...
  )
  reads <- reads[reads$cell %in% cells, ]
  if (!assume.sorted) {
    if (verbose) {
      message("Sorting fragments")
    }
    reads <- reads[with(data = reads, expr = order(chr, start)), ]
  }
  if (verbose) {
    message("Writing output")
  }
  fwrite(
    x = reads,
    file = output.path,
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE,
    sep = '\t'
  )
  rm(reads)
  invisible(x = gc())
  if (compress) {
    if (verbose) {
      message("Compressing output")
    }
    outf <- bgzip(file = output.path)
    if (file.exists(outf)) {
      file.remove(output.path)
    }
    if (index) {
      if (verbose) {
        message("Building index")
      }
      index.file <- indexTabix(
        file = paste0(outf), format = 'bed', zeroBased = TRUE
      )
    }
  }
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
#' @examples
#' FindTopFeatures(object = atac_small[['peaks']][])
FindTopFeatures.default <- function(
  object,
  assay = NULL,
  min.cutoff = 'q5',
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
#' @examples
#' FindTopFeatures(object = atac_small[['peaks']])
FindTopFeatures.Assay <- function(
  object,
  assay = NULL,
  min.cutoff = 'q5',
  verbose = TRUE,
  ...
) {
  data.use <- GetAssayData(object = object, slot = 'counts')
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
      x = sub(pattern = "q",replacement = "", x = as.character(x = min.cutoff))
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
#' @method FindTopFeatures Seurat
#' @examples
#' FindTopFeatures(atac_small)
FindTopFeatures.Seurat <- function(
  object,
  assay = NULL,
  min.cutoff = 'q5',
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
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @examples
#' FRiP(object = atac_small, peak.assay = 'peaks', bin.assay = 'bins')
FRiP <- function(
  object,
  peak.assay,
  bin.assay,
  chromosome = 'chr1',
  verbose = TRUE
) {
  if (verbose) {
    message('Calculating fraction of reads in peaks per cell')
  }
  peak.data <- GetAssayData(
    object = object, assay = peak.assay, slot = 'counts'
  )
  bin.data <- GetAssayData(
    object = object, assay = bin.assay, slot = 'counts'
  )
  if (!is.null(x = chromosome)) {
    peak.data <- peak.data[grepl(
      pattern = paste0('^', chromosome, '\\-|^', chromosome, ':'),
      x = rownames(x = peak.data)
    ), ]
    bin.data <- bin.data[grepl(
      pattern = paste0('^', chromosome, '\\-|^', chromosome, ':'),
      x = rownames(x = bin.data)
    ), ]
  }
  peak.counts <- colSums(x = peak.data)
  bin.counts <- colSums(x = bin.data)
  frip <- peak.counts / bin.counts
  object <- AddMetaData(object = object, metadata = frip, col.name = 'FRiP')
  return(object)
}

#' GenomeBinMatrix
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
#' @return Returns a sparse matrix
#' @examples
#' gn <- 780007
#' names(gn) <- 'chr1'
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' GenomeBinMatrix(
#'   fragments = fpath,
#'   genome = gn,
#'   binsize = 1000,
#'   chunk = 1
#' )
GenomeBinMatrix <- function(
  fragments,
  genome,
  cells = NULL,
  binsize = 5000,
  chunk = 50,
  sep = c('-', '-'),
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

globalVariables(names = 'cell', package = 'Signac')
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
#' @param ... Additional arguments passed to \code{\link{GetReadsInRegion}}
#'
#' @importFrom dplyr group_by summarize
#' @importFrom stats ecdf
#'
#' @return Returns a \code{\link[Seurat]{Seurat}} object with
#' added metadata for the ratio of mononucleosomal to nucleosome-free fragments
#' per cell, and the percentile rank of each ratio.
#' @export
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' atac_small <- SetFragments(object = atac_small, file = fpath)
#' NucleosomeSignal(object = atac_small)
NucleosomeSignal <- function(
  object,
  assay = NULL,
  region = 'chr1-1-249250621',
  min.threshold = 147,
  max.threshold = 294,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  fragments.use <- GetReadsInRegion(
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
  fragments.use <- as.data.frame(x = fragments.use[, c('cell', 'length')])
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
#' @param sep A length-2 character vector containing the separators to be used
#' when constructing genomic coordinates from the regions. The first element is
#' used to separate the chromosome from the genomic coordinates, and the second
#' element used to separate the start and end coordinates.
#' @importFrom BiocGenerics width
#' @importMethodsFrom GenomicRanges width
#' @rdname RegionStats
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RegionStats(
#' object = rownames(atac_small),
#' genome = BSgenome.Hsapiens.UCSC.hg19, sep = c(":", "-")
#' )
#' }
RegionStats.default <- function(
  object,
  genome,
  sep = c('-', '-'),
  verbose = TRUE,
  ...
) {
  if (!requireNamespace('BSgenome', quietly = TRUE)) {
    stop("Please install BSgenome: BiocManager::install('BSgenome')")
  }
  if (!requireNamespace('Biostrings', quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }
  if (inherits(x = object, what = 'character')) {
    object <- StringToGRanges(regions = object, sep = sep)
  }
  sequence.length <- width(x = object)
  sequences <- BSgenome::getSeq(x = genome, names = object)
  gc <- Biostrings::letterFrequency(
    x = sequences, letters = 'CG'
  ) / sequence.length * 100
  colnames(gc) <- 'GC.percent'
  dinuc <- Biostrings::dinucleotideFrequency(sequences)
  sequence.stats <- cbind(dinuc, gc, sequence.length)
  rownames(sequence.stats) <- GRangesToString(grange = object, sep = sep)
  return(sequence.stats)
}

#' @rdname RegionStats
#' @method RegionStats Assay
#' @importFrom methods slot
#' @importFrom Seurat GetAssayData
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' RegionStats(
#' object = atac_small[['peaks']],
#' genome = BSgenome.Hsapiens.UCSC.hg19, sep = c(":", "-")
#' )
#' }
RegionStats.Assay <- function(
  object,
  genome,
  sep = c('-', '-'),
  verbose = TRUE,
  ...
) {
  regions <- rownames(x = object)
  feature.metadata <- RegionStats(
    object = regions,
    genome = genome,
    sep = sep,
    verbose = verbose,
    ...
  )
  meta.data <- GetAssayData(object = object, slot = 'meta.features')
  meta.data <- cbind(meta.data, feature.metadata)
  slot(object = object, name = 'meta.features') <- meta.data
  return(object)
}

#' @param assay Name of assay to use
#' @rdname RegionStats
#' @method RegionStats Seurat
#' @export
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
  sep = c('-', '-'),
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- RegionStats(
    object = assay.data,
    genome = genome,
    sep = sep,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' @param method Which TF-IDF implementation to use. Choice of:
#' \itemize{
#'  \item{1}: The LSI implementation used by Stuart & Butler et al. 2019
#'  (\url{https://doi.org/10.1101/460147}).
#'  \item{2}: The standard LSI implementation used by Cusanovich & Hill
#'  et al. 2018 (\url{https://doi.org/10.1016/j.cell.2018.06.052}).
#'  \item{3}: The log-TF method
#'  \item{4}: The 10x Genomics method (no TF normalization)
#' }
#' @param scale.factor Which scale factor to use. Default is 10000.
#' @param verbose Print progress
#' @rdname RunTFIDF
#' @importFrom Matrix colSums rowSums Diagonal tcrossprod
#' @importFrom methods is "slot<-" slot
#' @export
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
  if (method == 4) {
    tf <- object
  } else {
    tf <- tcrossprod(x = object, y = Diagonal(x = 1 / npeaks))
  }
  idf <- ncol(x = object) / rowSums(x = object)
  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    slot(object = tf, name = 'x') <- log1p(
      x = slot(object = tf, name = 'x') * scale.factor
    )
    idf <- log(1 + idf)
  }
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  if (method == 1) {
    slot(object = norm.data, name = 'x') <- log1p(
      x = slot(object = norm.data, name = 'x') * scale.factor
    )
  }
  colnames(x = norm.data) <- colnames(x = object)
  rownames(x = norm.data) <- rownames(x = object)
  return(norm.data)
}

#' @rdname RunTFIDF
#' @method RunTFIDF Assay
#' @export
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
    object = GetAssayData(object = object, slot = 'counts'),
    method = method,
    assay = assay,
    scale.factor = scale.factor,
    verbose = verbose,
    ...
  )
  new.data <- as(object = new.data, Class = 'dgCMatrix')
  object <- SetAssayData(
    object = object,
    slot = 'data',
    new.data = new.data
  )
  return(object)
}

#' @param assay Name of assay to use
#' @rdname RunTFIDF
#' @method RunTFIDF Seurat
#' @export
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
#' @param tss.positions A GRanges object containing the TSS positions
#' @param cells A vector of cells to include. If NULL (default), use all cells
#' in the object
#' @param verbose Display messages
#' @importFrom Matrix rowMeans
#' @importFrom methods slot
#' @importFrom stats ecdf
#'
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @export
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v75)
#' gene.ranges <- genes(EnsDb.Hsapiens.v75)
#' gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
#' tss.ranges <- GRanges(
#'   seqnames = seqnames(gene.ranges),
#'   ranges = IRanges(start = start(gene.ranges), width = 2),
#'   strand = strand(gene.ranges)
#' )
#' seqlevelsStyle(tss.ranges) <- 'UCSC'
#' tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
#'
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' atac_small <- SetFragments(object = atac_small, file = fpath)
#' TSSEnrichment(object = atac_small, tss.positions = tss.ranges[1:100])
#' }
TSSEnrichment <- function(
  object,
  tss.positions,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  cutmatrix <- CreateRegionPileupMatrix(
    object = object,
    regions = tss.positions,
    upstream = 1000,
    downstream = 1000,
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

  # store the normalized TSS matrix. For now put it in misc
  object <- AddToMisc(
    object = object,
    assay = assay,
    new.data = norm.matrix,
    save.as = 'TSS.enrichment.matrix'
  )
  return(object)
}

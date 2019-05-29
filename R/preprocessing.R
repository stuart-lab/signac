#' @param verbose Display messages
#' @rdname BinarizeCounts
#' @export
#'
BinarizeCounts.default <- function(
  object,
  assay = NULL,
  verbose = TRUE
) {
  if (class(object) == 'dgCMatrix') {
    object@x <- rep.int(x = 1, times = length(object@x))
  } else {
    object[object > 1] <- 1
  }
  return(object)
}

#' @rdname BinarizeCounts
#' @method BinarizeCounts Assay
#' @importFrom Seurat GetAssayData SetAssayData
#' @export
BinarizeCounts.Assay <- function(
  object,
  assay = NULL,
  verbose = TRUE
) {
  data.matrix <- GetAssayData(object = object, slot = 'counts')
  object <- SetAssayData(
    object = object,
    slot = 'counts',
    new.data = BinarizeCounts(object = data.matrix, assay = assay, verbose = verbose)
  )
  return(object)
}

#' @param assay Name of assay to use. Can be a list of assays,
#' and binarization will be applied to each.
#' @rdname BinarizeCounts
#' @method BinarizeCounts Seurat
#' @importFrom Seurat GetAssay
#' @export
BinarizeCounts.Seurat <- function(
  object,
  assay = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object)
  for (i in 1:length(x = assay)) {
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

#' Convert a peak matrix to a gene activity matrix
#'
#' This function will take in a peak matrix and an annotation file (gtf) and collapse the peak
#' matrix to a gene activity matrix. It makes the simplifying assumption that all counts in the gene
#' body plus X kb up and or downstream should be attributed to that gene.
#'
#' @param peak.matrix Matrix of peak counts
#' @param annotation.file Path to GTF annotation file
#' @param seq.levels Which seqlevels to keep (corresponds to chromosomes usually)
#' @param include.body Include the gene body?
#' @param upstream Number of bases upstream to consider
#' @param downstream Number of bases downstream to consider
#' @param verbose Print progress/messages
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame distanceToNearest
#' @importFrom rtracklayer import mcols
#' @importFrom GenomeInfoDb keepSeqlevels seqlevelsStyle
#' @importFrom SummarizedExperiment promoters
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix colSums
#'
#' @export
#'
CreateGeneActivityMatrix <- function(
  peak.matrix,
  annotation.file,
  seq.levels = c(1:22, "X", "Y"),
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE
) {
  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- makeGRangesFromDataFrame(df = peak.df)

  # get annotation file, select genes
  gtf <- import(con = annotation.file)
  gtf <- keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')
  seqlevelsStyle(gtf) <- "UCSC"
  gtf.genes <- gtf[gtf$type == 'gene']

  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[subjectHits(x = keep.overlaps)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- paste0(peak.ids$seqnames, ":", peak.ids$start, "-", peak.ids$end)
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')

  # collapse into expression matrix
  peak.matrix <- as(object = peak.matrix, Class = 'matrix')
  all.features <- unique(x = annotations$new_feature)

  if (PlanThreads() > 1) {
    mysapply <- future_sapply
  } else {
    mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  }
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){
    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(colSums(x = submat))
    } else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = 'dgCMatrix'))
}

#' CreateMotifMatrix
#'
#' Create a motif x feature matrix from a set of genomic ranges,
#' the genome, and a set of position weight matrices
#'
#' @param features A GRanges object containing a set of genomic features
#' @param pwm A PFMatrixList object containing position weight matrices to use
#' @param genome Any object compatible with the \code{genome} argument in \code{\link{matchMotifs}}
#' @param sep A length-2 character vector containing the separators to be used when constructing
#' matrix rownames from the GRanges
#' @param ... Additional arguments passed to \code{\link{matchMotifs}}
#'
#' @return Returns a sparse matrix
#' @importFrom motifmatchr matchMotifs motifMatches
#' @export
CreateMotifMatrix <- function(
  features,
  pwm,
  genome,
  sep = c("-", "-"),
  ...
) {
  motif_ix <- matchMotifs(pwms = PFMatrixList, subject = features, genome = genome)
  motif.matrix <- motifMatches(object = motif_ix)
  motif.matrix <- as(Class = 'dgCMatrix', object = motif.matrix)
  rownames(motif.matrix) <- GRangesToString(grange = features, sep = sep)
  return(motif.matrix)
}

#' FeatureMatrix
#'
#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments Path to tabix-indexed fragments file
#' @param features A GRanges object containing a set of genomic intervals. These will form the rows of the
#' matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param chunk Number of chunks to use when processing the fragments file. Fewer chunks may enable faster processing,
#'  but will use more memory.
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome
#' and coordinates, second separator is used to separate start and end coordinates.
#' @param verbose Display messages
#'
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pblapply
#' @importFrom Matrix sparseMatrix
#' @importFrom Rsamtools TabixFile seqnamesTabix
#'
#' @export
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
    value = seqnamesTabix(file = tbx),
    pruning.mode = "coarse"
  )
  feature.list <- ChunkGRanges(
    granges = features,
    nchunk = chunk
  )
  if (verbose) {
    message('Extracting reads overlapping genome bins')
  }
  if (PlanThreads() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  cells.in.regions <- mylapply(
    X = feature.list,
    FUN = GetCellsInRegion,
    tabix = tbx,
    cells = cells
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
  rownames(featmat) <- names(x = feature.lookup)
  colnames(featmat) <- names(x = cell.lookup)
  if (!is.null(x = cells)) {
    cells.accept <- intersect(cells, colnames(x = featmat))
    return(featmat[, cells.accept])
  } else {
    return(featmat)
  }
}

#' FilterFragments
#'
#' Remove cells from a fragments file that are not present in a given list of cells.
#' Note that this reads the whole fragments file into memory, so may require a lot of memory
#' depending on the size of the fragments file.
#'
#' @param fragment.path Path to a tabix-indexed fragments file
#' @param cells A vector of cells to retain
#' @param output.path Name and path for output tabix file. A tabix index file will also be created in the same location, with
#' the .tbi file extension.
#' @param assume.sorted Assume sorted input and don't sort the filtered file. Can save a lot of time, but indexing will
#' fail if assumption is wrong.
#' @param compress Compress filtered fragments using bgzip (default TRUE)
#' @param index Index the filtered tabix file (default TRUE)
#' @param verbose Display messages
#'
#' @importFrom data.table fread fwrite
#' @importFrom Rsamtools indexTabix bgzip
#'
#' @export
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
    showProgress = verbose
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
      index.file <- indexTabix(file = paste0(outf), format = 'bed', zeroBased = TRUE)
    }
  }
}

#' @importFrom Matrix rowSums
#' @importFrom stats ecdf
#' @rdname FindTopFeatures
#' @export
FindTopFeatures.default <- function(
  object,
  assay = NULL,
  min.cutoff = 'q5',
  verbose = TRUE,
  ...
) {
  featurecounts <- rowSums(object)
  e.dist <- ecdf(featurecounts)
  hvf.info <- data.frame(
    row.names = names(featurecounts),
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
  if (is.null(min.cutoff)) {
    VariableFeatures(object) <- rownames(hvf.info)
  } else if (is.numeric(min.cutoff)) {
    VariableFeatures(object) <- rownames(hvf.info[hvf.info$count > min.cutoff, ])
  } else {
    percentile.use <- as.numeric(x = sub(pattern = "q", replacement = "", x = as.character(x = min.cutoff)))/100
    VariableFeatures(object) <- rownames(hvf.info[hvf.info$percentile > percentile.use, ])
  }
  return(object)
}

#' @rdname FindTopFeatures
#' @importFrom Seurat DefaultAssay GetAssay
#' @export
#' @method FindTopFeatures Seurat
FindTopFeatures.Seurat <- function(
  object,
  assay = NULL,
  min.cutoff = 'q5',
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
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
#' @param verbose Display messages
#'
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData AddMetaData
#'
#' @export
#'
FRiP <- function(
  object,
  peak.assay,
  bin.assay,
  verbose = TRUE
) {
  if (verbose) {
    message('Calculating fraction of reads in peaks per cell')
  }
  peak.counts <- colSums(x = GetAssayData(object = object, assay = peak.assay, slot = 'counts'))
  bin.counts <-colSums(x = GetAssayData(object = object, assay = bin.assay, slot = 'counts'))
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
#' @param genome A vector of chromosome sizes for the genome. This is used to construct the
#' genome bin coordinates. The can be obtained by calling \code{\link[GenomeInfoDb]{seqlengths}} on
#' a \code{\link[BSgenome]{BSgenome-class}} object.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param binsize Size of the genome bins to use
#' @param chunk Number of chunks to use when processing the fragments file. Fewer chunks may enable faster processing,
#'  but will use more memory.
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome
#' and coordinates, second separator is used to separate start and end coordinates.
#' @param verbose Display messages
#'
#' @importFrom GenomicRanges tileGenome
#'
#' @export
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

#' Calculate periodogram score per cell
#'
#' @param object A Seurat object
#' @param fragments Path to an indexed fragment file
#'
#' @export
#'
Period <- function(
  object,
  fragments = NULL
) {
  # TODO
  return(object)
}

#' @param method Which TF-IDF implementation to use. Choice of:
#' \itemize{
#'  \item{1}: The LSI implementation used by Stuart & Butler et al. 2019 (\url{https://doi.org/10.1101/460147}).
#'  \item{2}: The standard LSI implementation used by Cusanovich & Hill et al. 2018 (\url{https://doi.org/10.1016/j.cell.2018.06.052}).
#'  \item{3}: The log-TF method
#'  \item{4}: The 10x Genomics method (no TF normalization)
#' }
#' @param scale.factor Which scale factor to use. Default is 10000.
#' @param verbose Print progress
#' @rdname RunTFIDF
#' @importFrom Matrix colSums rowSums t Diagonal
#' @export
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
#' mat_norm <- RunTFIDF(object = mat)
#'
RunTFIDF.default <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  if (class(x = object) == "data.frame") {
    object <- as.matrix(x = object)
  }
  if (class(x = object) != "dgCMatrix") {
    object <- as(object = object, Class = "dgCMatrix")
  }
  if (verbose) {
    message("Performing TF-IDF normalization")
  }
  npeaks <- colSums(x = object)
  if (method == 4) {
    tf <- object
  } else {
    tf <- t(x = t(x = object) / npeaks)
  }
  idf <- ncol(x = object) / rowSums(x = object)
  if (method == 2) {
    idf <- log(1 + idf)
  } else if (method == 3) {
    tf@x <- log1p(x = tf@x * scale.factor)
    idf <- log(1 + idf)
  }
  norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  if (method == 1) {
    norm.data@x <- log1p(x = norm.data@x * scale.factor)
  }
  colnames(x = norm.data) <- colnames(x = object)
  rownames(x = norm.data) <- rownames(x = object)
  return(norm.data)
}

#' @rdname RunTFIDF
#' @method RunTFIDF Assay
#' @export
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
RunTFIDF.Seurat <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
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

#' Set the fragments file path for creating plots
#'
#' Give path of indexed fragments file that goes with data in the object.
#' Checks for a valid path and an index file with the same name (.tbi) at the same path.
#' Stores the path under the tools slot for access by visualization functions.
#' One fragments file can be stored for each assay.
#'
#' @param object A Seurat object
#' @param file Path to indexed fragment file. See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @param assay Assay used to generate the fragments. If NULL, use the active assay.
#'
#' @export
#'
SetFragments <- function(
  object,
  file,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!(assay %in% names(x = slot(object = object, name = 'assays')))) {
    stop('Requested assay not present in object')
  }
  index.file <- paste0(file, '.tbi')
  if (all(file.exists(file, index.file))) {
    file <- normalizePath(path = file)
    current.tools <- slot(object = object, name = 'tools')
    current.tools$fragments[[assay]] <- file
    slot(object = object, name = 'tools') <- current.tools
    return(object)
  } else {
    stop('Requested file does not exist or is not indexed')
  }
}

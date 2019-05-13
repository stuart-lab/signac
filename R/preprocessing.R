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

#' @param assay Name of assay to use
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
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- BinarizeCounts(
    object = assay.data,
    assay = assay,
    verbose = verbose
  )
  object[[assay]] <- assay.data
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
  if (!PackageCheck('GenomicRanges', error = FALSE)) {
    stop("Please install GenomicRanges from Bioconductor.")
  }
  if (!PackageCheck('rtracklayer', error = FALSE)) {
    stop("Please install rtracklayer from Bioconductor.")
  }

  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)

  # get annotation file, select genes
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')
  GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
  gtf.genes <- gtf[gtf$type == 'gene']

  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
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
      return(Matrix::colSums(x = submat))
    } else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = 'dgCMatrix'))
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
#' mat_norm <- RunTFIDF(data = mat)
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

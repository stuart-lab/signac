#' @include generics.R
#'
NULL

#' Link peaks to genes
#'
#' Find peaks that are predictive of gene expression. Fits a generalized linear
#' model with elastic net regularization between peak accessibility and gene
#' expression, for all peaks within a given distance threshold (default 200 kb).
#' Returns a peak x gene sparse matrix, with each entry being the correlation
#' coefficient between the accessibility of the peak and expression of the gene.
#'
#' @param object A Seurat object
#' @param peak.assay Name of assay containing peak information
#' @param expression.assay Name of assay containing gene expression information
#' @param gene.coords GRanges object containing coordinates of genes in the
#' expression assay. If NULL, extract from gene annotations stored in the assay.
#' @param binary Binarize the peak counts (default FALSE)
#' @param alpha Alpha parameter for elastic net regularization. Controls balance
#' between LASSO (alpha=1) and ridge (alpha=0) regression.
#' @param distance Distance threshold for peaks to include in regression model
#' @param min.cells Minimum number of cells positive for the peak and gene
#' needed to include in the results.
#' @param expression.threshold Minimum value for a gene to be classified as
#' expressed. Only used for filtering genes from regression.
#' @param keep.all Retain all peaks in the output. Useful for knowing what peaks
#' had 0 coefficient and what were not tested due to expression or accessibility
#' thresholds.
#' @param genes.use Genes to test. If NULL, determine from expression assay.
#' @param verbose Display messages
#'
#' @importFrom Seurat GetAssayData
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#'
#' @return Returns a sparse matrix
#' @export
LinkPeaks <- function(
  object,
  peak.assay,
  expression.assay,
  gene.coords = NULL,
  binary = FALSE,
  alpha = 0.95,
  distance = 5e+05,
  min.cells = 10,
  expression.threshold = 0.1,
  genes.use = NULL,
  keep.all = FALSE,
  verbose = TRUE
) {
  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(
      ranges = Annotation(object = object[[peak.assay]])
    )
  }
  peak.data <- GetAssayData(
    object = object, assay = peak.assay, slot = 'counts'
  )
  expression.data <- GetAssayData(
    object = object, assay = expression.assay, slot = 'data'
  )
  peakcounts <- rowSums(x = peak.data > 0)
  genecounts <- rowSums(x = expression.data > expression.threshold)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  if (is.null(x = genes.use)) {
    expression.data <- expression.data[genes.keep, ]
  } else {
    genes.keep <- intersect(
      x = names(x = genes.keep[genes.keep]), y = genes.use
    )
    expression.data <- expression.data[genes.keep, ]
  }
  if (verbose) {
    message(
      "Testing ",
      nrow(x = expression.data),
      " genes and ",
      sum(peaks.keep),
      " peaks"
    )
  }
  if (binary) {
    peak.data <- BinarizeCounts(object = peak.data)
  }
  genes <- rownames(x = expression.data)
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
  peaks <- granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToGene(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
  genes.use <- colnames(x = peak_distance_matrix)
  coef.vec <- c()
  gene.vec <- c()
  if (verbose) {
    pb <- txtProgressBar(
      min = 1,
      max = length(x = genes.use),
      style = 3,
      file = stderr()
    )
  }
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }

  # run in parallel across genes
  res <- mylapply(
    X = seq_along(along.with = genes.use),
    FUN = function(i) {
      peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
      gene.expression <- expression.data[genes.use[[i]], ]
      if (sum(peak.use) < 2) {
        # no peaks close to gene
        if (verbose) setTxtProgressBar(pb = pb, value = i)
        next
      } else {
        peak.access <- t(x = peak.data[peak.use, ])
        cvfit <- cv.glmnet(
          x = peak.access,
          y = gene.expression,
          alpha = alpha,
          family = 'poisson'
        )
        lambda <- cvfit$lambda.1se
        coef.results <- coef(object = cvfit, s = lambda)
        coef.results <- coef.results[2:nrow(x = coef.results), ]
        if (keep.all) {
          # keep all genes so we know what was tested
          coef.results.filtered <- coef.results
        } else {
          coef.results.filtered <- coef.results[coef.results != 0]
        }
        if (length(x = coef.results.filtered) == 0) {
          if (verbose) setTxtProgressBar(pb = pb, value = i)
          next
        } else {
          gene.vec <- c(gene.vec, rep(i, length(x = coef.results.filtered)))
          coef.vec <- c(coef.vec, coef.results.filtered)
          if (verbose) setTxtProgressBar(pb = pb, value = i)
        }
      }
      return(list("gene" = gene.vec, "coef" = coef.vec))
    }
  )
  # combine results
  gene.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 1))
  coef.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 2))

  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(NULL)
  }
  peak.key <- seq_along(
    along.with = unique(x = names(x = coef.vec))
  )
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = coef.vec)],
    x = coef.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  return(coef.matrix)
}

#' Find peaks near genes
#'
#' Find peaks that are within a given distance threshold to each gene
#'
#' @param peaks A GRanges object containing peak coordinates
#' @param genes A GRanges object containing gene coordinates
#' @param distance Distance threshold. Peaks within this distance from the gene
#' will be recorded.
#' @param sep Separator for peak names when creating results matrix
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#'
#' @return Returns a sparse matrix
#'
DistanceToGene <- function(peaks, genes, distance = 200000, sep = c("-", "-")) {
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = genes, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

#' Find top links
#'
#' Filter a set of promoter-enhancer links and retain the top X\%
#' based on correlation score. Links below the percentile cutoff
#' will be set to zero. Note that the cutoff is defined across
#' all genes, so some genes may have no links after filtering.
#' The cutoff value is determined considering non-zero link
#' scores only (eg, top 50\% among non-zero scores).
#'
#' @param links A sparse matrix containing genes in the rows and enhancers in
#' the columns.
#' Values correspond to the link score.
#' @param cutoff Percentile cutoff for retaining links. 0.50 (default) retains
#' the top 50\% of links.
#' @return Returns a sparse matrix
#' @export
TopLinks <- function(links, cutoff = 0.25) {
  nonzero.values <- abs(x = links@x)
  cutoff.value <- quantile(x = nonzero.values, 1 - cutoff)
  links[abs(links) < cutoff.value] <- 0
  return(links)
}

# Return names of nonzero elements for each matrix row
#
# @param x A matrix
# @param keep.sparse Keep as a sparse matrix. Slower but uses less memory.
# @return Return list where each element correponds to a row
# of the input matrix, and contains a vector of the names of the
# nonzero elements for the row
GetLinkedElements <- function(x, keep.sparse = FALSE) {
  if (keep.sparse) {
    results <- list()
    for (i in seq_along(along.with = rownames(x = x))) {
      results[[rownames(x = x)[[i]]]] <- names(x = which(abs(x[i, ]) > 0))
    }
    return(results)
  } else {
    return(
      apply(
        X = x,
        MARGIN = 1,
        FUN = function(y) names(x = which(abs(x = y) > 0))
        )
      )
  }
}

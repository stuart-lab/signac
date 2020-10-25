#' @include generics.R
#'
NULL

#' Get peaks linked to genes
#'
#' Find peaks linked to a given set of genes
#'
#' @param object A Seurat object
#' @param features A list of genes to find linked peaks for
#' @param assay Name of assay to use. If NULL, use the default assay
#' @param min.abs.score Minimum absolute value of the link score for a link to
#' be returned
#' @export
#' @concept links
#' @seealso GetLinkedGenes
GetLinkedPeaks <- function(
  object,
  features,
  assay = NULL,
  min.abs.score = 0
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  lnk <- Links(object = object[[assay]])
  if (length(x = lnk) == 0) {
    stop("No links present in assay. Run LinkPeaks first.")
  }
  lnk.keep <- lnk[(abs(x = lnk$score) > min.abs.score) & lnk$gene %in% features]
  return(unique(x = lnk.keep$peak))
}

#' Get genes linked to peaks
#'
#' Find genes linked to a given set of peaks
#'
#' @param object A Seurat object
#' @param features A list of peaks to find linked genes for
#' @param assay Name of assay to use. If NULL, use the default assay
#' @param min.abs.score Minimum absolute value of the link score for a link to
#' be returned
#' @export
#' @concept links
#' @seealso GetLinkedPeaks
GetLinkedGenes <- function(
  object,
  features,
  assay = NULL,
  min.abs.score = 0
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  lnk <- Links(object = object[[assay]])
  if (length(x = lnk) == 0) {
    stop("No links present in assay. Run LinkPeaks first.")
  }
  lnk.keep <- lnk[(abs(x = lnk$score) > min.abs.score) & lnk$peak %in% features]
  return(unique(x = lnk.keep$gene))
}

#' Cicero connections to links
#'
#' Convert the output of Cicero connections to a set of genomic ranges where
#' the start and end coordinates of the range are the midpoints of the linked
#' elements. Only elements on the same chromosome are included in the output.
#'
#' @param conns A dataframe containing co-accessible elements. This would
#' usually be the output of \code{\link[cicero]{run_cicero}} or
#' \code{\link[cicero]{assemble_connections}}. Specifically, this should be a
#' dataframe where the first column contains the genomic coordinates of the
#' first element in the linked pair of elements, with chromosome, start, end
#' coordinates separated by "-" characters. The second column should be the
#' second element in the linked pair, formatted in the same way as the first
#' column. A third column should contain the co-accessibility scores.
#' @param ccans This is optional, but if supplied should be a dataframe
#' containing the cis-co-accessibility network (CCAN) information generated
#' by \code{\link[cicero]{generate_ccans}}. Specifically, this should be a
#' dataframe containing the name of the peak in the first column, and the
#' CCAN that it belongs to in the second column.
#' @param threshold Threshold for retaining a coaccessible site. Links with
#' a value less than or equal to this threshold will be discarded.
#'
#' @export
#' @importFrom GenomicRanges start end makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stringi stri_split_fixed
#'
#' @concept links
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
ConnectionsToLinks <- function(conns, ccans = NULL, threshold = 0) {
  # add chromosome information
  chr1 <- stri_split_fixed(str = conns$Peak1, pattern = "-")
  conns$chr1 <- unlist(x = chr1)[3 * (seq_along(along.with = chr1)) - 2]
  chr2 <- stri_split_fixed(str = conns$Peak2, pattern = "-")
  conns$chr2 <- unlist(x = chr2)[3 * (seq_along(along.with = chr2)) - 2]

  # filter out trans-chr links
  conns <- conns[conns$chr1 == conns$chr2, ]

  # filter based on threshold
  conns <- conns[!is.na(conns$coaccess), ]
  conns <- conns[conns$coaccess > threshold, ]

  # add group information
  if (!is.null(x = ccans)) {
    ccan.lookup <- ccans$CCAN
    names(x = ccan.lookup) <- ccans$Peak
    groups <- as.vector(x = ifelse(
      test = is.na(x = ccan.lookup[conns$Peak1]),
      yes = ccan.lookup[conns$Peak2],
      no = ccan.lookup[conns$Peak1]
    )
    )
    conns$group <- groups
  } else {
    conns$group <- NA
  }

  # extract genomic regions
  coords.1 <- StringToGRanges(regions = conns$Peak1, sep = c("-", "-"))
  coords.2 <- StringToGRanges(regions = conns$Peak2, sep = c("-", "-"))
  chr <- seqnames(x = coords.1)

  # find midpoints
  midpoints.1 <- start(x = coords.1) + (width(x = coords.1) / 2)
  midpoints.2 <- start(x = coords.2) + (width(x = coords.2) / 2)

  startpos <- ifelse(
    test = midpoints.1 > midpoints.2,
    yes = midpoints.2,
    no = midpoints.1
  )
  endpos <- ifelse(
    test = midpoints.1 > midpoints.2,
    yes = midpoints.1,
    no = midpoints.2
  )

  link.df <- data.frame(chromosome = chr,
                        start = startpos,
                        end = endpos,
                        score = conns$coaccess,
                        group = conns$group)

  # convert to genomic ranges
  links <- makeGRangesFromDataFrame(df = link.df, keep.extra.columns = TRUE)
  return(links)
}

#' Link peaks to genes
#'
#' Find peaks that are correlated with gene expression. Fits a generalized
#' linear model with elastic net regularization between peak accessibility and
#' gene expression for all peaks within a given distance threshold from the gene
#' TSS.
#'
#' @param object A Seurat object
#' @param peak.assay Name of assay containing peak information
#' @param expression.assay Name of assay containing gene expression information
#' @param expression.slot Name of slot to pull expression data from
#' @param gene.coords GRanges object containing coordinates of genes in the
#' expression assay. If NULL, extract from gene annotations stored in the assay.
#' @param binary Binarize the peak counts (default FALSE)
#' @param alpha Alpha parameter for elastic net regularization. Controls balance
#' between LASSO (alpha=1) and ridge (alpha=0) regularization.
#' @param distance Distance threshold for peaks to include in regression model
#' @param min.cells Minimum number of cells positive for the peak and gene
#' needed to include in the results.
#' @param expression.threshold Minimum value for a gene to be classified as
#' expressed. Only used for filtering genes from regression.
#' @param keep.all Retain all peaks in the output. Useful for knowing what peaks
#' had 0 coefficient and what were not tested due to expression or accessibility
#' thresholds.
#' @param genes.use Genes to test. If NULL, determine from expression assay.
#' @param family Response type (for gene expression values). Valid values are
#' "gaussian", "poisson", "binomial", or a \code{\link[stats]{glm}}-family
#' object
#' @param null Compute null distribution of coefficients using randomly selected
#' peaks from other chromosomes.
#' @param n_sample Number of peaks to sample at random when computing the null
#' distribution.
#' @param verbose Display messages
#'
#' @importFrom Seurat GetAssayData
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats coef
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#'
#' @return Returns a Seurat object with the \code{Links} information set, with
#' each recorded link being the correlation coefficient between the
#' accessibility of the peak and expression of the gene.
#'
#' @export
#' @concept links
LinkPeaks <- function(
  object,
  peak.assay,
  expression.assay,
  expression.slot = "data",
  gene.coords = NULL,
  binary = FALSE,
  alpha = 0.95,
  distance = 5e+05,
  min.cells = 10,
  expression.threshold = 0.1,
  family = "poisson",
  genes.use = NULL,
  keep.all = FALSE,
  null = FALSE,
  n_sample = 300,
  verbose = TRUE
) {
  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(
      ranges = Annotation(object = object[[peak.assay]])
    )
  }
  meta.features <- GetAssayData(
    object = object, assay = peak.assay, slot = "meta.features"
  )
  peak.data <- GetAssayData(
    object = object, assay = peak.assay, slot = 'counts'
  )
  expression.data <- GetAssayData(
    object = object, assay = expression.assay, slot = expression.slot
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
  peak_distance_matrix <- DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)

  coef.vec <- c()
  gene.vec <- c()
  pval.vec <- c()
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
      gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))

      if (sum(peak.use) < 2) {
        # no peaks close to gene
        return(list("gene" = NULL, "coef" = NULL, "pval" = NULL))
      } else {
        peak.access <- t(x = peak.data[peak.use, ])
        cvfit <- cv.glmnet(
          x = peak.access,
          y = gene.expression,
          alpha = alpha,
          family = family
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
          return(list("gene" = NULL, "coef" = NULL, "pval" = NULL))
        } else {
          if (null) {
            # compute null distribution with same lambda
            # sample from peaks on a different chromosome to the gene
            trans.peaks <- all.peaks[
              !grepl(pattern = paste0("^", gene.chrom), x = all.peaks)
            ]
            peaks.test <- colnames(x = peak.access)
            meta.use <- meta.features[trans.peaks, ]
            meta.use <- rbind(meta.use, meta.features[peaks.test, ])
            # sample matching total accessibility and GC
            peak.random <- MatchRegionStats(
              meta.feature = meta.use,
              regions = peaks.test,
              features.match = c("GC.percent", "count"),
              n = n_sample,
              verbose = FALSE
            )
            peak.access <- t(x = peak.data[peak.random, ])
            cvfit <- glmnet(
              x = peak.access,
              y = gene.expression,
              alpha = alpha,
              family = family
            )
            coef.results.rand <- coef(object = cvfit, s = lambda)
            coef.results.rand <- coef.results.rand[
              2:nrow(x = coef.results.rand),
            ]
            # number of random peaks with coefficient at least as extreme as obs
            pvals <- sapply(
              X = coef.results.filtered,
              FUN = function(x) {
                if (x < 0) {
                  sum(coef.results.rand < x)
                } else {
                  sum(coef.results.rand > x)
                }
              })
            pvals <- pvals / n_sample
            pval.vec <- c(pval.vec, pvals)
          }
          gene.vec <- c(gene.vec, rep(i, length(x = coef.results.filtered)))
          coef.vec <- c(coef.vec, coef.results.filtered)
        }
      }
      return(list("gene" = gene.vec, "coef" = coef.vec, "pval" = pval.vec))
    }
  )
  # combine results
  gene.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 1))
  coef.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 2))
  if (null) {
    pval.vec <- do.call(what = c, args = sapply(X = res, FUN = `[[`, 3))
  }

  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
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
  links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  if (null) {
    # add pvalues
    pv.matrix <- sparseMatrix(
      i = gene.vec,
      j = peak.key[names(x = pval.vec)],
      x = pval.vec,
      dims = c(length(x = genes.use), max(peak.key))
    )
    rownames(x = pv.matrix) <- genes.use
    colnames(x = pv.matrix) <- names(x = peak.key)
    pv.lnk <- LinksToGRanges(linkmat = pv.matrix, gene.coords = gene.coords.use)
    links$pvalue <- pv.lnk$score
  }
  Links(object = object[[peak.assay]]) <- links
  return(object)
}

### Not exported ###

# Link matrix to granges
#
# Create set of genomic ranges from a sparse matrix containing links
#
# @param linkmat A sparse matrix with genes in the rows and peaks in the
# columns
# @param gene.coords Genomic coordinates for each gene
# @return Returns a GRanges object
#' @importFrom GenomicRanges resize start width GRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics sort
LinksToGRanges <- function(linkmat, gene.coords, sep = c("-", "-")) {
  # get TSS for each gene
  tss <- resize(gene.coords, width = 1, fix = 'start')
  gene.idx <- sapply(
    X = rownames(x = linkmat),
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx]

  # get midpoint of each peak
  peak.ranges <- StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)

  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "dgTMatrix")

  # create dataframe
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )

  # work out start and end coords
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL

  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}


# Find peaks near genes
#
# Find peaks that are within a given distance threshold to each gene
#
# @param peaks A GRanges object containing peak coordinates
# @param genes A GRanges object containing gene coordinates
# @param distance Distance threshold. Peaks within this distance from the gene
# will be recorded.
# @param sep Separator for peak names when creating results matrix
#
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges resize
#
# @return Returns a sparse matrix
DistanceToTSS <- function(peaks, genes, distance = 200000, sep = c("-", "-")) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
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

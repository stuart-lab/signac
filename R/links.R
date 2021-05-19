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
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
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
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
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
#' See the Cicero package for more information:
#' \url{https://bioconductor.org/packages/cicero/}
#'
#' @param conns A dataframe containing co-accessible elements. This would
#' usually be the output of \code{run_cicero} or
#' \code{assemble_connections}. Specifically, this should be a
#' dataframe where the first column contains the genomic coordinates of the
#' first element in the linked pair of elements, with chromosome, start, end
#' coordinates separated by "-" characters. The second column should be the
#' second element in the linked pair, formatted in the same way as the first
#' column. A third column should contain the co-accessibility scores.
#' @param ccans This is optional, but if supplied should be a dataframe
#' containing the cis-co-accessibility network (CCAN) information generated
#' by \code{generate_ccans}. Specifically, this should be a
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
#' Find peaks that are correlated with the expression of nearby genes.
#' For each gene, this function computes the correlation coefficient between
#' the gene expression and accessibility of each peak within a given distance
#' from the gene TSS, and computes an expected correlation coefficient for each
#' peak given the GC content, accessibility, and length of the peak. The expected
#' coefficient values for the peak are then used to compute a z-score and p-value.
#'
#' This function was inspired by the method originally described by SHARE-seq
#' (Sai Ma et al. 2020, Cell). Please consider citing the original SHARE-seq
#' work if using this function: \doi{10.1016/j.cell.2020.09.056}
#'
#' @param object A Seurat object
#' @param peak.assay Name of assay containing peak information
#' @param expression.assay Name of assay containing gene expression information
#' @param expression.slot Name of slot to pull expression data from
#' @param gene.coords GRanges object containing coordinates of genes in the
#' expression assay. If NULL, extract from gene annotations stored in the assay.
#' @param distance Distance threshold for peaks to include in regression model
#' @param min.distance Minimum distance between peak and TSS to include in
#' regression model. If NULL (default), no minimum distance is used.
#' @param min.cells Minimum number of cells positive for the peak and gene
#' needed to include in the results.
#' @param genes.use Genes to test. If NULL, determine from expression assay.
#' @param method Which correlation coefficient to compute. Can be "pearson"
#' (default), "spearman", or "kendall".
#' @param n_sample Number of peaks to sample at random when computing the null
#' distribution.
#' @param pvalue_cutoff Minimum p-value required to retain a link. Links with a
#' p-value equal or greater than this value will be removed from the output.
#' @param score_cutoff Minimum absolute value correlation coefficient for a link
#' to be retained
#' @param verbose Display messages
#'
#' @importFrom Seurat GetAssayData
#' @importFrom stats pnorm sd
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom qlcMatrix corSparse
#' @importMethodsFrom Matrix t
#'
#' @return Returns a Seurat object with the \code{Links} information set. This is
#' a \code{\link[GenomicRanges]{granges}} object accessible via the \code{\link{Links}}
#' function, with the following information:
#' \itemize{
#'   \item{score: the correlation coefficient between the accessibility of the
#'   peak and expression of the gene}
#'   \item{zscore: the z-score of the correlation coefficient, computed based on
#'   the distribution of correlation coefficients from a set of background peaks}
#'   \item{pvalue: the p-value associated with the z-score for the link}
#'   \item{gene: name of the linked gene}
#'   \item{peak: name of the linked peak}
#' }
#'
#' @export
#' @concept links
LinkPeaks <- function(
  object,
  peak.assay,
  expression.assay,
  expression.slot = "data",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  method = "pearson",
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  verbose = TRUE
) {
  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("min.distance should be a numeric value")
    }
    if (min.distance < 0) {
      warning("Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    } else if (min.distance == 0) {
      min.distance <- NULL
    }
  }

  if (is.null(x = gene.coords)) {
    gene.coords <- CollapseToLongestTranscript(
      ranges = Annotation(object = object[[peak.assay]])
    )
  }
  meta.features <- GetAssayData(
    object = object, assay = peak.assay, slot = "meta.features"
  )
  features.match <- c("GC.percent", "count")
  if (!("GC.percent" %in% colnames(x = meta.features))) {
    stop("GC content per peak has not been computed.\n",
         "Run RegionsStats before calling this function.")
  }
  peak.data <- GetAssayData(
    object = object, assay = peak.assay, slot = 'counts'
  )
  if (!("count" %in% colnames(x = meta.features))) {
    # compute total count
    hvf.info <- FindTopFeatures(object = peak.data)
    hvf.info <- hvf.info[rownames(x = meta.features), ]
    meta.features <- cbind(meta.features, hvf.info)
  }
  expression.data <- GetAssayData(
    object = object, assay = expression.assay, slot = expression.slot
  )
  peakcounts <- meta.features[rownames(x = peak.data), "count"]
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  if (is.null(x = genes.use)) {
    expression.data <- expression.data[genes.keep, ]
  } else {
    genes.keep <- intersect(
      x = names(x = genes.keep[genes.keep]), y = genes.use
    )
    expression.data <- expression.data[genes.keep, , drop = FALSE]
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
  genes <- rownames(x = expression.data)
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
  peaks <- granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
  if (!is.null(x = min.distance)) {
    peak_distance_matrix_min <- DistanceToTSS(
      peaks = peaks,
      genes = gene.coords.use,
      distance = min.distance
    )
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }
  if (sum(peak_distance_matrix) == 0) {
    stop("No peaks fall within distance threshold\n",
         "Have you set the proper genome and seqlevelsStyle for ",
         peak.assay,
         " assay?")
  }
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)

  peak.data <- t(x = peak.data)

  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
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
      gene.expression <- t(x = expression.data[genes.use[[i]], , drop = FALSE])
      gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))

      if (sum(peak.use) < 2) {
        # no peaks close to gene
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      } else {
        peak.access <- peak.data[, peak.use, drop = FALSE]
        coef.result <- corSparse(
          X = peak.access,
          Y = gene.expression
        )
        rownames(x = coef.result) <- colnames(x = peak.access)
        coef.result <- coef.result[abs(x = coef.result) > score_cutoff, , drop = FALSE]

        if (nrow(x = coef.result) == 0) {
          return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
        } else {

          # select peaks at random with matching GC content and accessibility
          # sample from peaks on a different chromosome to the gene
          peaks.test <- rownames(x = coef.result)
          trans.peaks <- all.peaks[
            !grepl(pattern = paste0("^", gene.chrom), x = all.peaks)
          ]
          meta.use <- meta.features[trans.peaks, ]
          pk.use <- meta.features[peaks.test, ]
          bg.peaks <- lapply(
            X = seq_len(length.out = nrow(x = pk.use)),
            FUN = function(x) {
              MatchRegionStats(
                meta.feature = meta.use,
                query.feature = pk.use[x, , drop = FALSE],
                features.match = c("GC.percent", "count", "sequence.length"),
                n = n_sample,
                verbose = FALSE
              )
            }
          )
          # run background correlations
          bg.access <- peak.data[, unlist(x = bg.peaks), drop = FALSE]
          bg.coef <- corSparse(
            X = bg.access,
            Y = gene.expression
          )
          rownames(bg.coef) <- colnames(bg.access)
          zscores <- vector(mode = "numeric", length = length(x = peaks.test))
          for (j in seq_along(along.with = peaks.test)) {
            coef.use <- bg.coef[(((j - 1) * n_sample) + 1):(j * n_sample), ]
            z <- (coef.result[j] - mean(x = coef.use)) / sd(x = coef.use)
            zscores[[j]] <- z
          }
          names(x = coef.result) <- peaks.test
          names(x = zscores) <- peaks.test
          zscore.vec <- c(zscore.vec, zscores)
          gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
          coef.vec <- c(coef.vec, coef.result)
        }
        gc(verbose = FALSE)
        pval.vec <- pnorm(q = -abs(x = zscore.vec))
        links.keep <- pval.vec < pvalue_cutoff
        if (sum(x = links.keep) == 0) {
          return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
        } else {
          gene.vec <- gene.vec[links.keep]
          coef.vec <- coef.vec[links.keep]
          zscore.vec <- zscore.vec[links.keep]
          return(list("gene" = gene.vec, "coef" = coef.vec, "zscore" = zscore.vec))
        }
      }
    }
  )
  # combine results
  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 3))

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
  # add zscores
  z.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = zscore.vec)],
    x = zscore.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
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

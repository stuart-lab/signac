#' @include generics.R
#'
NULL

#' @param assay Name of assay to use. If NULL, use the default assay
#' @importFrom SeuratObject DefaultAssay
#' @method GetLinkedPeaks Seurat
#' @concept links
#' @rdname GetLinkedPeaks
#' @export
GetLinkedPeaks.Seurat <- function(
    object,
    features,
    assay = NULL,
    min.abs.score = 0,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  links <- GetLinkedPeaks(
    object = object[[assay]],
    features = features,
    min.abs.score = min.abs.score,
    ...
  )
  return(links)
}

#' @export
#' @method GetLinkedPeaks Assay5
#' @concept links
#' @rdname GetLinkedPeaks
GetLinkedPeaks.Assay5 <- function(
    object,
    ...
) {
  stop("GetLinkedPeaks requires a GRangesAssay object")
}

#' @param features A list of genes to find linked peaks for
#' @param min.abs.score Minimum absolute value of the link score for a link to
#' be returned
#' @export
#' @method GetLinkedPeaks GRangesAssay
#' @concept links
#' @rdname GetLinkedPeaks
GetLinkedPeaks.GRangesAssay <- function(
  object,
  features,
  min.abs.score = 0,
  ...
) {
  lnk <- Links(object = object)
  if (length(x = lnk) == 0) {
    stop("No links present in assay. Run LinkPeaks first.")
  }
  lnk.keep <- lnk[(abs(x = lnk$score) > min.abs.score) & lnk$gene %in% features]
  return(unique(x = lnk.keep$peak))
}

#' @param assay Name of assay to use. If NULL, use the default assay
#' @importFrom SeuratObject DefaultAssay
#' @method GetLinkedGenes Seurat
#' @concept links
#' @rdname GetLinkedGenes
#' @export
GetLinkedGenes.Seurat <- function(
    object,
    features,
    assay = NULL,
    min.abs.score = 0,
    ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  links <- GetLinkedGenes(
    object = object[[assay]],
    features = features,
    min.abs.score = min.abs.score,
    ...
  )
  return(links)
}

#' @export
#' @method GetLinkedGenes Assay5
#' @concept links
#' @rdname GetLinkedGenes
GetLinkedGenes.Assay5 <- function(
    object,
    ...
) {
  stop("GetLinkedGenes requires a GRangesAssay object")
}

#' @param features A list of peaks to find linked genes for
#' @param min.abs.score Minimum absolute value of the link score for a link to
#' be returned
#' @export
#' @method GetLinkedGenes GRangesAssay
#' @concept links
#' @rdname GetLinkedGenes
GetLinkedGenes.GRangesAssay <- function(
    object,
    features,
    min.abs.score = 0,
    ...
) {
  lnk <- Links(object = object)
  if (length(x = lnk) == 0) {
    stop("No links present in assay. Run LinkPeaks first.")
  }
  lnk.keep <- lnk[(abs(x = lnk$score) > min.abs.score) & lnk$peak %in% features]
  return(unique(x = lnk.keep$gene))
}

#' Cicero connections to links
#'
#' Convert the output of Cicero connections to a
#' `InteractionSet::GInteractions()` object.
#'
#' See the Cicero package for more information:
#' <https://bioconductor.org/packages/cicero/>
#'
#' @param conns A dataframe containing co-accessible elements. This would
#' usually be the output of `run_cicero` or
#' `assemble_connections`. Specifically, this should be a
#' dataframe where the first column contains the genomic coordinates of the
#' first element in the linked pair of elements, with chromosome, start, end
#' coordinates separated by "-" characters. The second column should be the
#' second element in the linked pair, formatted in the same way as the first
#' column. A third column should contain the co-accessibility scores.
#' @param ccans This is optional, but if supplied should be a dataframe
#' containing the cis-co-accessibility network (CCAN) information generated
#' by `generate_ccans`. Specifically, this should be a
#' dataframe containing the name of the peak in the first column, and the
#' CCAN that it belongs to in the second column.
#' @param threshold Threshold for retaining a coaccessible site. Links with
#' a value less than or equal to this threshold will be discarded.
#'
#' @export
#' @importFrom InteractionSet GInteractions
#'
#' @concept links
#' @return Returns a [InteractionSet::GInteractions()] object
ConnectionsToLinks <- function(
  conns,
  ccans = NULL,
  threshold = 0
) {
  
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
  
  # filter based on threshold
  conns <- conns[!is.na(conns$coaccess), ]
  conns <- conns[conns$coaccess > threshold, ]
  
  # create ginteractions
  gi <- GInteractions(GRanges(conns$Peak1), GRanges(conns$Peak2))
  gi$score <- conns$coaccess
  gi$group <- conns$group

  return(gi)
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
#' @param peak.layer Name of layer to pull chromatin data from
#' @param expression.assay Name of assay containing gene expression information
#' @param expression.layer Name of layer to pull expression data from
#' @param method Correlation method to use. One of "pearson" or "spearman"
#' @param key Key to use when storing link information in the assay
#' @param gene.coords GRanges object containing coordinates of genes in the
#' expression assay. If NULL, extract from gene annotations stored in the assay.
#' @param distance Distance threshold for peaks to include in regression model
#' @param min.distance Minimum distance between peak and TSS to include in
#' regression model. If NULL (default), no minimum distance is used.
#' @param min.cells Minimum number of cells positive for the peak and gene
#' needed to include in the results.
#' @param genes.use Genes to test. If NULL, determine from expression assay.
#' @param n_sample Number of peaks to sample at random when computing the null
#' distribution.
#' @param pvalue_cutoff Minimum p-value required to retain a link. Links with a
#' p-value equal or greater than this value will be removed from the output.
#' @param score_cutoff Minimum absolute value correlation coefficient for a link
#' to be retained
#' @param gene.id Set to TRUE if genes in the expression assay are named
#' using gene IDs rather than gene names.
#' @param verbose Display messages
#' @param peak.slot Deprecated (use `peak.layer`)
#' @param expression.slot Deprecated (used `expression.layer`)
#'
#' @importFrom SeuratObject GetAssayData Layers
#' @importFrom stats pnorm sd
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom lifecycle is_present deprecated deprecate_warn
#' @importMethodsFrom Matrix t
#'
#' @return Returns a Seurat object with the `Links` information set. This is
#' a [GenomicRanges::granges()] object accessible via the [Links()]
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
  peak.layer = "counts",
  expression.layer = "data",
  method = "pearson",
  key = "linkpeaks",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE,
  peak.slot = deprecated(),
  expression.slot = deprecated()
) {
  if (!inherits(x = object[[peak.assay]], what = "GRangesAssay")) {
    stop("The requested assay is not a GRangesAssay")
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
  if (is_present(arg = expression.slot)) {
    deprecate_warn(
      when = "1.16.0",
      what = "LinkPeaks(expression.slot)",
      with = "LinkPeaks(expression.layer)"
    )
    expression.layer <- expression.slot
  }
  if (is_present(arg = peak.slot)) {
    deprecate_warn(
      when = "1.16.0",
      what = "LinkPeaks(peak.slot)",
      with = "LinkPeaks(peak.layer)"
    )
    peak.layer <- peak.slot
  }
  features.match <- c("GC.percent", "count", "sequence.length")
  if (method == "pearson") {
    cor_method <- corSparse
  } else if (method == "spearman") {
    cor_method <- SparseSpearmanCor
  } else {
    stop("method can be one of 'pearson' or 'spearman'.")
  }

  if (is.null(x = gene.coords)) {
    annot <- Annotation(object = object[[peak.assay]])
    if (is.null(x = annot)) {
      stop("Gene annotations not found")
    }
    gene.coords <- CollapseToLongestTranscript(
      ranges = annot
    )
  }
  meta.features <- object[[peak.assay]][[]]
  if (!(all(
    c("GC.percent", "sequence.length") %in% colnames(x = meta.features)
    ))) {
    stop("DNA sequence information for each peak has not been computed.\n",
         "Run RegionStats before calling this function.")
  }
  if (!("count" %in% colnames(x = meta.features))) {
    data.use <- GetAssayData(object = object[[peak.assay]], layer = peak.layer)
    hvf.info <- FindTopFeatures(object = data.use, verbose = FALSE)
    hvf.info <- hvf.info[rownames(meta.features), , drop = FALSE]
    meta.features <- cbind(meta.features, hvf.info)
  }
  if (!(peak.layer %in% Layers(object = object[[peak.assay]]))) {
    stop("Requested peak layer not found")
  }
  peak.data <- LayerData(object = object[[peak.assay]], layer = peak.layer)
  if (!(expression.layer %in% Layers(object = object[[expression.assay]]))) {
    stop("Requested expression layer not found")
  }
  expression.data <- LayerData(
    object = object[[expression.assay]], layer = expression.layer
  )
  peakcounts <- rowSums(x = peak.data > 0)
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  if (!is.null(x = genes.use)) {
    genes.keep <- intersect(
      x = names(x = genes.keep[genes.keep]), y = genes.use
    )
  }
  expression.data <- expression.data[genes.keep, , drop = FALSE]
  genes <- rownames(x = expression.data)
  if (gene.id) {
    gene.coords.use <- gene.coords[gene.coords$gene_id %in% genes,]
    gene.coords.use$gene_name <- gene.coords.use$gene_id
  } else {
    gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
  }
  if (length(x = gene.coords.use) == 0) {
    stop("Could not find gene coordinates for requested genes")
  }
  if (length(x = gene.coords.use) < nrow(x = expression.data)) {
    message("Found gene coordinates for ",
            length(x = gene.coords.use),
            " genes")
  }
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
  if (verbose) {
    message(
      "Testing ",
      nrow(x = expression.data),
      " genes and ",
      sum(rowSums(x = peak_distance_matrix) > 0),
      " peaks"
    )
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
        coef.result <- cor_method(
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
            !grepl(pattern = paste0("^", gene.chrom, "-"), x = all.peaks)
          ]
          meta.use <- meta.features[trans.peaks, ]
          pk.use <- meta.features[peaks.test, ]
          bg.peaks <- lapply(
            X = seq_len(length.out = nrow(x = pk.use)),
            FUN = function(x) {
              MatchRegionStats(
                meta.feature = meta.use,
                query.feature = pk.use[x, , drop = FALSE],
                features.match = features.match,
                n = n_sample,
                verbose = FALSE
              )
            }
          )
          # run background correlations
          bg.access <- peak.data[, unlist(x = bg.peaks), drop = FALSE]
          bg.coef <- cor_method(
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
  # links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  links <- LinksToGInteractions(
    linkmat = coef.matrix,
    gene.coords = gene.coords.use
  )
  # add zscores
  z.matrix <- sparseMatrix(
    i = gene.vec,
    j = peak.key[names(x = zscore.vec)],
    x = zscore.vec,
    dims = c(length(x = genes.use), max(peak.key))
  )
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  # z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
  z.lnk <- LinksToGInteractions(
    linkmat = z.matrix,
    gene.coords = gene.coords.use
  )
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  Links(object = object[[peak.assay]], key = key) <- links
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
#' @importFrom GenomicRanges resize start width makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics sort
LinksToGRanges <- function(linkmat, gene.coords) {
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
  peak.ranges <- GRanges(colnames(x = linkmat))
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)

  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "TsparseMatrix")

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

LinksToGInteractions <- function(linkmat, gene.coords) {
  x <- as(object = linkmat, Class = 'TsparseMatrix')
  peak.coords <- GRanges(colnames(x = linkmat))
  gene.coords$strand <- strand(x = gene.coords)
  region1 <- peak.coords[x@j+1]
  region2 <- gene.coords[x@i+1]
  gi <- GInteractions(region1, region2)
  gi$score <- x@x
  return(gi)
}

# Find peaks near genes
#
# Find peaks that are within a given distance threshold to each gene
#
# @param peaks A GRanges object containing peak coordinates
# @param genes A GRanges object containing gene coordinates
# @param distance Distance threshold. Peaks within this distance from the gene
# will be recorded.
#
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges resize
#
# @return Returns a sparse matrix
DistanceToTSS <- function(
  peaks,
  genes,
  distance = 200000
  ) {
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
  rownames(x = hit_matrix) <- as.character(x = peaks)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

LinkPeaks_BPCells <- function(
    object,
    peak.matrix.assay,
    peak.assay,
    expression.assay,
    peak.slot = "counts",
    expression.slot = "data",
    method = "pearson",
    gene.coords = NULL,
    distance = 5e+05,
    min.distance = NULL,
    min.cells = 10,
    genes.use = NULL,
    n_sample = 200,
    pvalue_cutoff = 0.05,
    score_cutoff = 0.05,
    gene.id = FALSE,
    verbose = TRUE
) {
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
  meta.features <- GetAssayData(
    object = object, assay = peak.assay, layer = "meta.features"
  )
  if (!(all(
    c("GC.percent", "sequence.length") %in% colnames(x = meta.features)
  ))) {
    stop("DNA sequence information for each peak has not been computed.\n",
         "Run RegionsStats before calling this function.")
  }
  if (!("count" %in% colnames(x = meta.features))) {
    data.use <- GetAssayData(object = object[[peak.matrix.assay]], layer = "counts")
    rownames(data.use) <- gsub(':', '-',   rownames(data.use) )
    hvf.info <- FindTopFeatures(object = data.use, verbose = FALSE)
    hvf.info <- hvf.info[rownames(meta.features), , drop = FALSE]
    meta.features <- cbind(meta.features, hvf.info)
  }
  peak.data <- GetAssayData(
    object = object, assay = peak.matrix.assay, layer = peak.slot
  )
  
  rownames(peak.data) <- gsub(':', '-',   rownames(peak.data) )
  expression.data <- GetAssayData(
    object = object, assay = expression.assay, layer = expression.slot
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



GetPromoterPeaksMatrix <- function(
    fragment_file,
    gene,
    upstream = 2000,
    downstream = 2000 
    ) {
  library(EnsDb.Hsapiens.v86)
  gene_info <- genes(EnsDb.Hsapiens.v86, filter = GeneNameFilter(gene))
  seqlevelsStyle(gene_info) <- "UCSC" 
  promoter_region <- promoters(x = gene_info, upstream = upstream, downstream = downstream)
  peak_mat <- BPCells::peak_matrix(frag, promoter_region, mode="insertions")
  return(peak_mat)
}

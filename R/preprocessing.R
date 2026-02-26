#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

#' @export
#' @concept qc
#' @rdname ATACqc
ATACqc.default <- function(
  object,
  annotations,
  fragtk.path = NULL,
  outdir = tempdir(),
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  # find fragtk
  fragtk.path <- fragtk.path %||% unname(obj = Sys.which(names = "fragtk"))
  if (nchar(x = fragtk.path) == 0) {
    stop(
      "fragtk not found. Please install fragtk:",
      "https://crates.io/crates/fragtk"
    )
  }
  if (!file.exists(fragtk.path)) {
    stop("fragtk executable does not exist at supplied path")
  }
  if (!dir.exists(paths = outdir)) {
    stop("Requested output directory does not exist")
  }

  object <- normalizePath(path = object, mustWork = TRUE)

  # get tss positions from annotations
  tss <- GetTSSPositions(ranges = annotations)

  # temp files
  tss.path <- tempfile(pattern = "signac_fragtk_tss", tmpdir = outdir)
  out.path <- tempfile(pattern = "signac_fragtk_qc", tmpdir = outdir)

  # write tss
  write.table(
    x = as.data.frame(x = tss),
    file = tss.path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # call fragtk qc
  cmd <- paste0(
    fragtk.path,
    " qc --fragments ",
    object,
    " --bed ",
    tss.path,
    " --outfile ",
    out.path
  )

  system(
    command = cmd,
    wait = TRUE,
    ignore.stderr = !verbose,
    ignore.stdout = !verbose
  )

  # load results
  md <- read.table(file = out.path, header = TRUE, row.names = 1, sep = "\t")

  # remove temp files
  if (cleanup) {
    files.to.remove <- c(tss.path, out.path)
    for (i in files.to.remove) {
      if (file.exists(i)) {
        file.remove(i)
      }
    }
  }

  return(md)
}

#' @export
#' @concept qc
#' @rdname ATACqc
#' @method ATACqc ChromatinAssay5
#' @importFrom utils write.table
ATACqc.ChromatinAssay5 <- function(
  object,
  annotations = NULL,
  fragtk.path = NULL,
  outdir = tempdir(),
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  annotations <- annotations %||% Annotation(object = object)
  if (is.null(x = annotations)) {
    stop("Gene annotation information not set")
  }
  frags <- Fragments(object = object)

  results <- data.frame()

  for (i in seq_along(along.with = frags)) {
    fragments <- GetFragmentData(object = frags[[i]], slot = "file.path")
    if (verbose) {
      message("Processing ", fragments)
    }

    md <- ATACqc(
      object = fragments,
      fragtk.path = fragtk.path,
      annotations = annotations,
      outdir = outdir,
      cleanup = cleanup,
      verbose = verbose,
      ...
    )

    # convert cell names
    cellconvert <- GetFragmentData(object = frags[[i]], slot = "cells")
    cc <- names(x = cellconvert)
    names(x = cc) <- cellconvert
    # in case some cells are missing
    cellconvert <- cellconvert[cellconvert %in% rownames(x = md)]
    md <- md[cellconvert, ]
    rownames(x = md) <- cc[rownames(x = md)]

    # concat across fragment files
    results <- rbind(results, md)
  }
  return(results)
}

#' @param assay Name of assay to use. If NULL, use the default assay.
#' @rdname ATACqc
#' @method ATACqc Seurat
#' @importFrom SeuratObject AddMetaData DefaultAssay
#' @export
#' @concept qc
ATACqc.Seurat <- function(
  object,
  assay = NULL,
  annotations = NULL,
  fragtk.path = NULL,
  outdir = tempdir(),
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  md <- ATACqc(
    object = object[[assay]],
    fragtk.path = fragtk.path,
    annotations = annotations,
    outdir = outdir,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  object <- AddMetaData(object = object, metadata = md)
  return(object)
}

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
  if (inherits(x = object, what = "CsparseMatrix")) {
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
#' @importFrom SeuratObject LayerData SetAssayData
#' @export
#' @concept preprocessing
BinarizeCounts.Assay <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  data.matrix <- LayerData(object = object, layer = "counts")
  object <- SetAssayData(
    object = object,
    layer = "counts",
    new.data = BinarizeCounts(
      object = data.matrix, assay = assay, verbose = verbose
    )
  )
  return(object)
}

#' @rdname BinarizeCounts
#' @method BinarizeCounts StdAssay
#' @importFrom SeuratObject LayerData LayerData<-
#' @export
#' @concept preprocessing
#' @examples
#' BinarizeCounts(atac_small[["peaks"]])
BinarizeCounts.StdAssay <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  data.matrix <- LayerData(object = object, layer = "counts")
  LayerData(object, layer = "counts") <- BinarizeCounts(
    object = data.matrix, assay = assay, verbose = verbose
  )
  return(object)
}

#' @param assay Name of assay to use. Can be a list of assays,
#' and binarization will be applied to each.
#' @rdname BinarizeCounts
#' @method BinarizeCounts Seurat
#' @importFrom SeuratObject DefaultAssay
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
  assay <- assay %||% DefaultAssay(object = object)
  for (i in seq_along(along.with = assay)) {
    assay.data <- object[[assay[[i]]]]
    assay.data <- BinarizeCounts(
      object = assay.data,
      assay = assay[[i]],
      verbose = verbose
    )
    object[[assay[[i]]]] <- assay.data
  }
  return(object)
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
#' @importFrom SeuratObject DefaultAssay "VariableFeatures<-"
#' @return Returns a [SeuratObject::Seurat()] object with
#' [SeuratObject::VariableFeatures()] set to the randomly sampled features.
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
  assay <- assay %||% DefaultAssay(object = object)
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
#' minimum percentile, for example 'q5' to set the top 95% most common features
#' as the VariableFeatures for the object. Alternatively, this can be an integer
#' specifying the minimum number of counts for the feature
#' to be included in the set of VariableFeatures. For example, setting to 10
#' will include features with >10 total counts in the set of VariableFeatures.
#' If `NULL`, include all features in VariableFeatures.
#' @param verbose Display messages
#'
#' @importFrom Matrix rowSums
#' @importFrom stats ecdf
#' @rdname FindTopFeatures
#' @export
#' @concept preprocessing
#' @examples
#' FindTopFeatures(object = atac_small[["peaks"]]["data"])
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
  return(hvf.info)
}

#' @param layer Name of layer to use
#' @rdname FindTopFeatures
#' @importFrom SeuratObject LayerData VariableFeatures
#' @importFrom utils packageVersion
#' @export
#' @method FindTopFeatures Assay5
#' @concept preprocessing
#' @examples
#' FindTopFeatures(object = atac_small[["peaks"]])
FindTopFeatures.Assay5 <- function(
  object,
  assay = NULL,
  layer = "counts",
  min.cutoff = "q5",
  key = "topfeatures",
  verbose = TRUE,
  ...
) {
  layer <- Layers(object = object, search = layer)
  for (i in seq_along(along.with = layer)) {
    if (isTRUE(x = verbose)) {
      message("Finding variable features for layer ", layer[i])
    }
    data.use <- LayerData(object = object, layer = layer[i], fast = TRUE)
    hvf <- FindTopFeatures(
      object = data.use,
      assay = assay,
      min.cutoff = min.cutoff,
      verbose = verbose,
      ...
    )
    rownames(x = hvf) <- Features(x = object, layer = layer[i])
    if (i == 1) {
      hvf.use <- hvf
    } else {
      # sum feature counts across layers
      hvf.use[rownames(x = hvf), ] <- hvf.use[rownames(x = hvf), ] + hvf
    }
  }
  # re-compute percentile due to multiple layers
  e.dist <- ecdf(x = hvf.use$count)
  hvf.use$percentile <- e.dist(hvf.use$count)
  hvf.use$rank <- rank(x = hvf.use$percentile)
  hvf.use$variable <- FALSE

  if (is.null(x = min.cutoff)) {
    hvf.use$variable <- TRUE
  } else if (is.numeric(x = min.cutoff)) {
    hvf.use[hvf.use[, 1] > min.cutoff, "variable"] <- TRUE
  } else {
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = min.cutoff))
    ) / 100
    hvf.use[hvf.use[, 2] > percentile.use, "variable"] <- TRUE
  }
  colnames(x = hvf.use) <- paste(
    "vf",
    key,
    layer[i],
    colnames(x = hvf.use),
    sep = "_"
  )
  object[["var.features"]] <- NULL
  object[["var.features.rank"]] <- NULL
  object[[names(x = hvf.use)]] <- hvf.use
  return(object)
}

#' @rdname FindTopFeatures
#' @importFrom utils packageVersion
#' @export
#' @method FindTopFeatures StdAssay
#' @concept preprocessing
#' @examples
#' FindTopFeatures(object = atac_small[["peaks"]])
FindTopFeatures.StdAssay <- function(
  object,
  assay = NULL,
  layer = "counts",
  min.cutoff = "q5",
  key = "topfeatures",
  verbose = TRUE,
  ...
) {
  FindTopFeatures.Assay5(
    object = object,
    assay = assay,
    layer = layer,
    min.cutoff = min.cutoff,
    key = key,
    verbose = verbose,
    ...
  )
}

#' @param key Key to use when storing the highly variable feature information in
#' the assay.
#' @rdname FindTopFeatures
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept preprocessing
#' @method FindTopFeatures Seurat
#' @examples
#' FindTopFeatures(atac_small)
FindTopFeatures.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  min.cutoff = "q5",
  key = "topfeatures",
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  assay.data <- object[[assay]]
  assay.data <- FindTopFeatures(
    object = assay.data,
    assay = assay,
    layer = layer,
    min.cutoff = min.cutoff,
    key = key,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' @rdname FitMeanVar
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param layer Name of layer to use. If NULL, use the default layer(s).
#' @param nfeatures Number of features to selected as top variable features.
#' @param loess.span `span` parameter passed to the [stats::loess()] function
#' @param min.cutoff Minimum number of counts for a feature to be eligible for
#' variable feature selection.
#' @param weight.mean How much to weight the ranking of features according to
#' their mean. Setting `weight.mean=0` will rank features according to their
#' residual variance only.
#' @param bins Number of bins to use when downsampling features across the range
#' of mean count values.
#' @param sample_per_bin Number of features to select per mean count bin in
#' feature downsampling step.
#' @param key Key to use when storing the highly variable feature information in
#' the assay.
#' @param verbose Display messages.
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept preprocessing
#' @method FitMeanVar Seurat
#' @examples
#' \dontrun{
#' FitMeanVar(atac_small)
#' }
FitMeanVar.Seurat <- function(
  object,
  assay = NULL,
  layer = NULL,
  nfeatures = 20000,
  loess.span = 0.1,
  min.cutoff = 10,
  weight.mean = 0.5,
  bins = 1000,
  sample_per_bin = 50,
  key = "dsLoess",
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  assay.data <- object[[assay]]
  assay.data <- FitMeanVar(
    object = assay.data,
    layer = layer,
    loess.span = loess.span,
    nfeatures = nfeatures,
    min.cutoff = min.cutoff,
    weight.mean = weight.mean,
    bins = bins,
    sample_per_bin = sample_per_bin,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' @rdname FitMeanVar
#' @importFrom SeuratObject Layers LayerData Features VariableFeatures
#' VariableFeatures<-
#' @export
#' @concept preprocessing
#' @method FitMeanVar Assay5
FitMeanVar.Assay5 <- function(
  object,
  layer = NULL,
  loess.span = 0.1,
  nfeatures = 20000,
  min.cutoff = 10,
  weight.mean = 0,
  bins = 1000,
  sample_per_bin = 50,
  key = "dsLoess",
  verbose = TRUE,
  ...
) {
  layer <- Layers(object = object, search = layer)
  feature.ranks <- list()
  for (i in seq_along(along.with = layer)) {
    if (isTRUE(x = verbose)) {
      message("Finding variable features for layer ", layer[i])
    }
    data <- LayerData(object = object, layer = layer[i], fast = TRUE)
    hvf <- FitMeanVar(
      object = data,
      loess.span = loess.span,
      min.cutoff = min.cutoff,
      weight.mean = weight.mean,
      bins = bins,
      nfeatures = nfeatures,
      sample_per_bin = sample_per_bin,
      verbose = verbose,
      ...
    )
    colnames(x = hvf) <- paste(
      "vf",
      key,
      layer[i],
      colnames(x = hvf),
      sep = "_"
    )
    rownames(x = hvf) <- Features(x = object, layer = layer[i])
    object[["var.features"]] <- NULL
    object[["var.features.rank"]] <- NULL
    object[[names(x = hvf)]] <- hvf
    feature.ranks[[i]] <- setNames(
      object = hvf[[paste("vf", key, layer[i], "rank", sep = "_")]],
      nm = rownames(x = hvf)
    )
  }
  # sum ranks
  feature.ranks <- unlist(x = feature.ranks)
  feature.ranks <- tapply(
    X = feature.ranks, INDEX = names(x = feature.ranks), FUN = sum
  )
  feature.ranks <- sort(x = feature.ranks, decreasing = FALSE)
  if (!is.na(x = nfeatures)) {
    top_features <- head(x = names(x = feature.ranks), n = nfeatures)
    VariableFeatures(object = object) <- top_features
  }
  return(object)
}

#' @rdname FitMeanVar
#' @param random.seed Random seed to set for sampling.
#' @importFrom SeuratObject DefaultAssay
#' @importFrom sparseMatrixStats rowVars
#' @importFrom stats loess predict
#' @export
#' @concept preprocessing
#' @method FitMeanVar default
FitMeanVar.default <- function(
  object,
  nfeatures = 20000,
  loess.span = 0.1,
  min.cutoff = 10,
  weight.mean = 0,
  bins = 1000,
  sample_per_bin = 50,
  random.seed = 1234,
  verbose = TRUE,
  ...
) {
  rs <- rowSums(x = object)

  if (is.character(x = min.cutoff)) {
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(min.cutoff))
    ) / 100
    count.thresh <- quantile(x = rs, probs = percentile.use, na.rm = TRUE)[[1]]
  } else {
    count.thresh <- min.cutoff
  }
  if (verbose) {
    message("Retained ", nrow(x = object), " features after count filtering")
  }
  if (nrow(x = object) == 0) {
    stop("No features remain after filtering by min.cutoff")
  }
  df <- data.frame(
    mean = rowMeans(x = object),
    variance = rowVars(x = object),
    variance.expected = 0,
    total.counts = rs
  )

  # run on a dataframe containing mean, variance, total.counts
  df <- FitMeanVar(
    object = df,
    loess.span = loess.span,
    weight.mean = weight.mean,
    bins = bins,
    sample_per_bin = sample_per_bin,
    random.seed = random.seed,
    verbose = verbose,
    ...
  )

  df$rank[df$total.counts < count.thresh] <- NA
  vf <- head(
    x = order(df$rank, decreasing = FALSE),
    n = nfeatures
  )
  df$variable <- FALSE
  df$variable[vf] <- TRUE
  return(df)
}

#' @rdname FitMeanVar
#' @importFrom sparseMatrixStats rowVars
#' @importFrom stats loess predict
#' @export
#' @concept preprocessing
#' @method FitMeanVar data.frame
FitMeanVar.data.frame <- function(
  object,
  loess.span = 0.1,
  weight.mean = 0,
  bins = 1000,
  sample_per_bin = 50,
  random.seed = 1234,
  verbose = TRUE,
  ...
) {
  if (!all(c("mean", "variance") %in% colnames(x = object))) {
    stop("Mean and variance information must be stored in the input dataframe")
  }
  set.seed(random.seed)
  object$log_mean <- log1p(x = object$mean)
  breaks <- seq(
    min(object$log_mean, na.rm = TRUE),
    max(object$log_mean, na.rm = TRUE),
    length.out = bins + 1
  )
  object$bin <- findInterval(
    x = object$log_mean,
    vec = breaks,
    rightmost.closed = TRUE
  )
  sampled_df <- do.call(
    what = rbind,
    args = lapply(X = split(object, object$bin), FUN = function(subset) {
      if (nrow(subset) > sample_per_bin) {
        subset <- subset[sample(
          x = nrow(x = subset),
          size = sample_per_bin,
          replace = FALSE
        ), ]
      }
      return(subset)
    })
  )
  loess_fit <- loess(
    formula = log1p(x = variance) ~ log_mean,
    data = sampled_df,
    span = loess.span
  )
  object$log.variance.expected <- predict(
    object = loess_fit, newdata = object$log_mean
  )
  object$variance.expected <- expm1(x = object$log.variance.expected)
  object$variance.residual <- log1p(x = object$variance) -
    object$log.variance.expected
  object$variance.residual[is.na(x = object$variance.residual)] <- 0

  object$residual.rank <- rank(
    x = -object$variance.residual, ties.method = "average"
  )
  object$mean.rank <- rank(x = -object$mean, ties.method = "average")
  object$rank <- (weight.mean * object$mean.rank) +
    ((1 - weight.mean) * object$residual.rank)
  return(object)
}

#' @param assay Name of assay to use
#' @param min.counts Minimum number of counts for feature to be eligible for
#' variable features
#' @param ncell.batch Number of cells to process in each batch. Higher number
#' increases speed but uses more memory.
#' @param nfeatures Number of top features to set as the variable features
#' @param theta Theta value for analytic Pearson residual calculation
#' @param verbose Display messages
#'
#' @importFrom Matrix rowMeans
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @rdname PearsonResidualVar
#' @export
#' @concept preprocessing
#' @examples
#' PearsonResidualVar(object = atac_small[["peaks"]]["counts"])
PearsonResidualVar.default <- function(
  object,
  assay = NULL,
  nfeatures = 20000,
  min.counts = 100,
  ncell.batch = 100,
  theta = 10,
  weight.mean = 0,
  verbose = TRUE,
  ...
) {
  # compute Pearson residual variance for each feature
  # low-memory implementation that does not construct the entire matrix of
  # cell x feature Pearson residuals

  # (X - μ) / sqrt(μ + μ²/θ)
  # clip at sqrt(n_cell)

  N <- ncol(x = object)
  clip_threshold <- sqrt(x = N)

  feature_means <- rowMeans(x = object)
  nonzero_mean <- feature_means > 0

  if (verbose) {
    message(
      "Retaining ",
      sum(nonzero_mean),
      " features with mean greater than zero"
    )
  }
  rn <- rownames(x = object)
  rcount <- rowSums(x = object)
  object <- object[nonzero_mean, ]
  feature_means.nonzero <- feature_means[nonzero_mean]

  denominator <- sqrt(
    feature_means.nonzero +
      ((feature_means.nonzero * feature_means.nonzero) / theta)
  )

  # iterate over the values for each feature,
  # compute the pearson residual variance
  resid_sums <- vector(mode = "numeric", length = nrow(x = object))
  resid_sum_square <- vector(mode = "numeric", length = nrow(x = object))
  nbatch <- ceiling(N / ncell.batch)
  if (nbatch <= 1) verbose <- FALSE

  if (verbose) pb <- txtProgressBar(min = 1, max = nbatch, style = 3)
  for (i in seq_len(length.out = nbatch)) {
    cells.interval.start <- 1 + ((i - 1) * ncell.batch)
    cells.interval.end <- min(N, (i * ncell.batch))

    resid <- as.matrix(
      x = (object[, cells.interval.start:cells.interval.end] -
             feature_means.nonzero) / denominator
    )
    resid[resid > clip_threshold] <- clip_threshold
    resid[resid < -clip_threshold] <- -clip_threshold
    rs <- rowSums(x = resid)
    resid_sums <- resid_sums + rs
    resid_sum_square <- resid_sum_square + rowSums(resid^2)
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
  }

  # Variance = [sum_of_squares - n * mean^2] / n
  resid_mean <- resid_sums / N
  pearson_residual_variance <- (resid_sum_square - (N * resid_mean^2)) / N
  resid.all <- rep(x = 0, length(x = feature_means))
  resid.all[nonzero_mean] <- pearson_residual_variance

  res_rank <- rank(x = -resid.all, ties.method = "average")
  mean_rank <- rank(x = -feature_means, ties.method = "average")

  # construct dataframe
  hvf.info <- data.frame(
    row.names = rn,
    count = rcount,
    mean = feature_means,
    ResidualVariance = resid.all,
    rank = (weight.mean * mean_rank) + ((1 - weight.mean) * res_rank)
  )

  return(hvf.info)
}

#' @rdname PearsonResidualVar
#' @importFrom SeuratObject VariableFeatures LayerData Layers
#' @importFrom utils packageVersion
#' @export
#' @method PearsonResidualVar Assay5
#' @concept preprocessing
#' @examples
#' PearsonResidualVar(object = atac_small[["peaks"]])
PearsonResidualVar.Assay5 <- function(
  object,
  assay = NULL,
  nfeatures = 20000,
  theta = 10,
  min.counts = 100,
  weight.mean = 0,
  ncell.batch = 100,
  key = "pearson",
  verbose = TRUE,
  ...
) {
  layer <- Layers(object = object, search = "counts")
  feature.ranks <- list()

  for (i in seq_along(along.with = layer)) {
    if (isTRUE(x = verbose)) {
      message("Finding variable features for layer ", layer[i])
    }
    data.use <- LayerData(object = object, layer = layer[i], fast = TRUE)
    hvf <- PearsonResidualVar(
      object = data.use,
      assay = assay,
      min.counts = min.counts,
      theta = theta,
      weight.mean = weight.mean,
      ncell.batch = ncell.batch,
      verbose = verbose,
      ...
    )
    colnames(x = hvf) <- paste(
      "vf",
      key,
      layer[i],
      colnames(x = hvf),
      sep = "_"
    )
    rownames(x = hvf) <- Features(x = object, layer = layer[i])
    object[["var.features"]] <- NULL
    object[["var.features.rank"]] <- NULL
    object[[names(x = hvf)]] <- hvf
    feature.ranks[[i]] <- setNames(
      object = hvf[[paste("vf", key, layer[i], "rank", sep = "_")]],
      nm = rownames(x = hvf)
    )
  }

  # sum ranks
  feature.ranks <- unlist(x = feature.ranks)
  feature.ranks <- tapply(
    X = feature.ranks, INDEX = names(x = feature.ranks), FUN = sum
  )
  feature.ranks <- sort(x = feature.ranks, decreasing = FALSE)
  if (!is.na(x = nfeatures)) {
    top_features <- head(x = names(x = feature.ranks), n = nfeatures)
    VariableFeatures(object = object) <- top_features
  }
  return(object)
}

#' @rdname PearsonResidualVar
#' @importFrom utils packageVersion
#' @export
#' @method PearsonResidualVar StdAssay
#' @concept preprocessing
#' @examples
#' PearsonResidualVar(object = atac_small[["peaks"]])
PearsonResidualVar.StdAssay <- function(
  object,
  assay = NULL,
  min.counts = 100,
  weight.mean = 0,
  theta = 10,
  ncell.batch = 100,
  key = "pearson",
  verbose = TRUE,
  ...
) {
  PearsonResidualVar.Assay5(
    object = object,
    assay = assay,
    min.counts = min.counts,
    weight.mean = weight.mean,
    theta = theta,
    ncell.batch = ncell.batch,
    key = key,
    verbose = verbose,
    ...
  )
}

#' @rdname PearsonResidualVar
#' @param weight.mean Weighting to apply to the feature mean relative to the
#' Pearson residual variance for ranking features. `weight.mean=0` will
#' rank features based on the Pearson residual variance only.
#' @param key Key to use when storing the highly variable feature information in
#' the assay.
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept preprocessing
#' @method PearsonResidualVar Seurat
#' @examples
#' PearsonResidualVar(atac_small)
PearsonResidualVar.Seurat <- function(
  object,
  assay = NULL,
  min.counts = 100,
  weight.mean = 0.5,
  theta = 10,
  ncell.batch = 100,
  key = "pearson",
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  assay.data <- object[[assay]]
  assay.data <- PearsonResidualVar(
    object = assay.data,
    assay = assay,
    min.counts = min.counts,
    weight.mean = weight.mean,
    ncell.batch = ncell.batch,
    theta = theta,
    key = key,
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
#' [CountFragments()] function.
#' @param col.name Name of column in metadata to store the FRiP information.
#' @param verbose Display messages
#'
#' @importFrom Matrix colSums
#' @importFrom SeuratObject LayerData AddMetaData
#'
#' @export
#' @concept qc
#' @return Returns a [SeuratObject::Seurat()] object
#' @examples
#' FRiP(object = atac_small, assay = "peaks", total.fragments = "fragments")
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
  peak.data <- LayerData(object = object, assay = assay, layer = "counts")
  total_fragments_cell <- object[[]][[total.fragments]]
  peak.counts <- colSums(x = peak.data)
  frip <- peak.counts / total_fragments_cell
  object <- AddMetaData(object = object, metadata = frip, col.name = col.name)
  return(object)
}

#' @param genome A `BSgenome` object or any other object supported by
#' `getSeq`. Do `showMethods("getSeq")` to get the list of all
#' supported object types.
#' @param verbose Display messages
#'
#' @importMethodsFrom GenomicRanges width
#' @rdname RegionStats
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' RegionStats(
#'   object = rownames(atac_small),
#'   genome = BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
RegionStats.default <- function(
  object,
  genome,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop("Please install BSgenome: BiocManager::install('BSgenome')")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }
  sequence.length <- width(x = object)
  common.seq <- intersect(x = seqlevels(x = object), y = seqlevels(x = genome))
  if (length(x = common.seq) < length(x = seqlevels(x = object))) {
    warning("Not all seqlevels present in supplied genome", immediate. = TRUE)
  }
  seq.keep <- as.character(x = seqnames(x = object)) %in% common.seq
  enum <- seq_along(along.with = seq.keep)
  object <- object[seq.keep]
  object <- keepSeqlevels(
    x = object, value = common.seq, pruning.mode = "coarse"
  )
  sequences <- Biostrings::getSeq(x = genome, object)
  gc <- Biostrings::letterFrequency(
    x = sequences, letters = "CG"
  ) / sequence.length[seq.keep] * 100
  colnames(gc) <- "GC.percent"
  dinuc <- Biostrings::dinucleotideFrequency(sequences)
  sequence.stats <- cbind(dinuc, gc)
  # fill missing seqnames with NA
  nadf <- as.data.frame(
    x = matrix(nrow = sum(!seq.keep), ncol = ncol(x = sequence.stats))
  )
  colnames(x = nadf) <- colnames(x = sequence.stats)
  rownames(x = nadf) <- enum[!seq.keep]
  rownames(x = sequence.stats) <- enum[seq.keep]
  sequence.stats <- rbind(sequence.stats, nadf)
  sequence.stats <- sequence.stats[enum, ]
  sequence.stats <- cbind(sequence.stats, sequence.length)
  return(sequence.stats)
}

#' @rdname RegionStats
#' @method RegionStats GRangesAssay
#' @importFrom methods slot
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' RegionStats(
#'   object = atac_small[["peaks"]],
#'   genome = BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
RegionStats.GRangesAssay <- function(
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
  meta.data <- object[[]]
  feature.metadata <- feature.metadata[rownames(x = meta.data), ]
  meta.data <- cbind(meta.data, feature.metadata)
  object[[names(x = meta.data)]] <- NULL
  object[[names(x = meta.data)]] <- meta.data
  return(object)
}

#' @param assay Name of assay to use
#' @rdname RegionStats
#' @method RegionStats Seurat
#' @export
#' @concept motifs
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' RegionStats(
#'   object = atac_small,
#'   assay = "bins",
#'   genome = BSgenome.Hsapiens.UCSC.hg38
#' )
#' }
RegionStats.Seurat <- function(
  object,
  genome,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- object[[assay]]
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
#' @param idf A precomputed IDF vector to use. If NULL, compute based on the
#' input data matrix.
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
  idf = NULL,
  verbose = TRUE,
  ...
) {
  if (inherits(x = object, what = "data.frame")) {
    object <- as.matrix(x = object)
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
    if (inherits(x = object, what = "IterableMatrix")) {
      tf <- BPCells::multiply_cols(mat = object, vec = 1 / npeaks)
    } else {
      tf <- tcrossprod(x = object, y = Diagonal(x = 1 / npeaks))
    }
  }
  if (!is.null(x = idf)) {
    precomputed_idf <- TRUE
    if (!inherits(x = idf, what = "numeric")) {
      stop("idf parameter must be a numeric vector")
    }
    if (length(x = idf) != nrow(x = object)) {
      stop(
        "Length of supplied IDF vector does not match",
        " number of rows in input matrix"
      )
    }
    if (any(idf == 0)) {
      stop("Supplied IDF values cannot be zero")
    }
    if (verbose) {
      message("Using precomputed IDF vector")
    }
  } else {
    precomputed_idf <- FALSE
    rsums <- rowSums(x = object)
    if (any(rsums == 0)) {
      warning("Some features contain 0 total counts")
    }
    idf <- ncol(x = object) / rsums
  }

  if (method == 2) {
    if (!precomputed_idf) {
      idf <- log(1 + idf)
    }
  } else if (method == 3) {
    slot(object = tf, name = "x") <- log1p(
      x = slot(object = tf, name = "x") * scale.factor
    )
    if (!precomputed_idf) {
      idf <- log(1 + idf)
    }
  }
  if (inherits(x = object, what = "IterableMatrix")) {
    norm.data <- BPCells::multiply_rows(mat = tf, vec = idf)
  } else {
    norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
  }
  if (method == 1) {
    if (inherits(x = norm.data, what = "IterableMatrix")) {
      norm.data <- log1p(scale.factor * norm.data)
    } else {
      slot(object = norm.data, name = "x") <- log1p(
        x = slot(object = norm.data, name = "x") * scale.factor
      )
    }
  }
  colnames(x = norm.data) <- colnames(x = object)
  rownames(x = norm.data) <- rownames(x = object)
  if (inherits(x = norm.data, what = "sparseMatrix")) {
    # set NA values to 0
    vals <- slot(object = norm.data, name = "x")
    vals[is.na(x = vals)] <- 0
    slot(object = norm.data, name = "x") <- vals
  }

  return(norm.data)
}

#' @rdname RunTFIDF
#' @method RunTFIDF Assay5
#' @export
#' @concept preprocessing
#' @importFrom SeuratObject LayerData Layers
#' @examples
#' RunTFIDF(atac_small[["peaks"]])
RunTFIDF.Assay5 <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  idf = NULL,
  layer = "counts",
  save = "data",
  verbose = TRUE,
  ...
) {
  olayer <- layer <- unique(x = layer)
  layer <- Layers(object = object, search = layer)
  if (length(x = save) != length(x = layer)) {
    save <- make.unique(names = gsub(
      pattern = olayer,
      replacement = save,
      x = layer
    ))
  }
  for (i in seq_along(along.with = layer)) {
    l <- layer[i]
    if (verbose) {
      message("Processing layer: ", l)
    }
    LayerData(
      object = object,
      layer = save[i],
      features = Features(x = object, layer = l),
      cells = Cells(x = object, layer = l)
    ) <- RunTFIDF(
      object = LayerData(object = object, layer = l, fast = NA),
      method = method,
      scale.factor = scale.factor,
      idf = idf,
      verbose = verbose,
      ...
    )
  }
  return(object)
}

#' @rdname RunTFIDF
#' @method RunTFIDF StdAssay
#' @export
#' @concept preprocessing
#' @examples
#' RunTFIDF(atac_small[["peaks"]])
RunTFIDF.StdAssay <- function(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 1e4,
  idf = NULL,
  layer = "counts",
  save = "data",
  verbose = TRUE,
  ...
) {
  RunTFIDF.Assay5(
    object = object,
    assay = assay,
    method = method,
    scale.factor = scale.factor,
    idf = idf,
    layer = layer,
    save = save,
    verbose = verbose,
    ...
  )
}

#' @param assay Name of assay to use.
#' @param layer Name of layer to use.
#' @param save Name of layer to save results in.
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
  idf = NULL,
  layer = "counts",
  save = "data",
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  assay.data <- object[[assay]]
  assay.data <- RunTFIDF(
    object = assay.data,
    assay = assay,
    method = method,
    scale.factor = scale.factor,
    idf = idf,
    layer = layer,
    save = save,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

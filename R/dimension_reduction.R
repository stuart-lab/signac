#' @include generics.R
#'
NULL

#' Calculate the Jaccard index between two matrices
#'
#' Finds the Jaccard similarity between rows of the two matrices. Note that
#' the matrices must be binary, and any rows with zero total counts will result
#' in an NaN entry that could cause problems in downstream analyses.
#'
#' This will calculate the raw Jaccard index, without normalizing for the
#' expected similarity between cells due to differences in sequencing depth.
#'
#' @param x The first matrix
#' @param y The second matrix
#'
#' @importFrom Matrix tcrossprod rowSums
#' @return Returns a matrix
#'
#' @export
#' @concept dimension_reduction
#' @examples
#' x <- matrix(data = sample(c(0, 1), size = 25, replace = TRUE), ncol = 5)
#' Jaccard(x = x, y = x)
Jaccard <- function(x, y) {
  if (any(x > 1) || any(y > 1)) {
    warning("Matrices contain values greater than 1.
            Please binarize matrices before running Jaccard")
  }
  intersection <- tcrossprod(x = x, y = y)
  union.counts.x <- rowSums(x = x)
  union.counts.y <- rowSums(x = y)
  A <- matrix(
    data = rep(x = union.counts.x, ncol(x = intersection)),
    ncol = ncol(x = intersection)
  )
  B <- matrix(
    data = rep(x = union.counts.y, nrow(x = intersection)),
    ncol = nrow(x = intersection)
  )
  jaccard.matrix <- as.matrix(x = intersection / ((A + t(B)) - intersection))
  return(jaccard.matrix)
}

#' @param assay Which assay to use. If NULL, use the default assay
#' @param n Number of singular values to compute
#' @param reduction.key Key for dimension reduction object
#' @param scale.max Clipping value for cell embeddings.
#' Default (NULL) is no clipping.
#' @param scale.embeddings Scale cell embeddings within each component to
#' mean 0 and SD 1 (default TRUE).
#' @param tol Tolerance (tol) parameter for [RSpectra::svds()].
#' Larger values speed up convergence due to greater amount of allowed error.
#' @param pca Run PCA. Setting this option to TRUE will perform implicit scaling
#' and centering of the input matrix to enable memory-efficient computation of
#' the principal components.
#' @param verbose Print messages
#'
#' @importFrom RSpectra svds
#' @importFrom SeuratObject CreateDimReducObject
#' @importMethodsFrom Matrix t
#'
#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @examples
#' x <- matrix(data = rnorm(100), ncol = 10)
#' RunSVD(x)
RunSVD.default <- function(
  object,
  assay = NULL,
  n = 50,
  scale.embeddings = !pca,
  pca = FALSE,
  reduction.key = ifelse(pca, "PCA_", "LSI_"),
  scale.max = NULL,
  verbose = TRUE,
  tol = 1e-05,
  ...
) {
  if (is.null(x = rownames(x = object))) {
    rownames(x = object) <- seq_len(length.out = nrow(x = object))
  }
  if (is.null(x = colnames(x = object))) {
    colnames(x = object) <- seq_len(length.out = ncol(x = object))
  }
  n <- min(n, (ncol(x = object) - 1))

  opts <- list("tol" = tol)
  if (pca) {
    if (verbose) {
      message("Running PCA")
    }
    # data needs to be standardized
    if (inherits(x = object, what = "IterableMatrix")) {
      # scale and center BPCells matrix
      # opts params in BPCells:::svds.IterableMatrix not implemented
      s <- BPCells::matrix_stats(object, row_stats = "variance")
      r_means <- s$row_stats["mean", ]
      r_vars <- s$row_stats["variance", ]
      object <- (object - r_means) / r_vars
    } else {
      opts <- c(opts, list("center" = TRUE, "scale" = TRUE))
    }
  } else {
    if (verbose) {
      message("Running SVD")
    }
  }

  components <- svds(A = t(x = object), k = n, opts = opts)
  feature.loadings <- components$v
  sdev <- components$d / sqrt(x = max(1, nrow(x = object) - 1))
  if (pca) {
    # weight by eigenvalues
    cell.embeddings <- components$u %*% diag(components$d)
  } else {
    cell.embeddings <- components$u
  }
  if (scale.embeddings) {
    if (verbose) {
      message("Scaling cell embeddings")
    }
    embed.mean <- apply(X = cell.embeddings, MARGIN = 2, FUN = mean)
    embed.sd <- apply(X = cell.embeddings, MARGIN = 2, FUN = sd)
    norm.embeddings <- t((t(cell.embeddings) - embed.mean) / embed.sd)
    if (!is.null(x = scale.max)) {
      norm.embeddings[norm.embeddings > scale.max] <- scale.max
      norm.embeddings[norm.embeddings < -scale.max] <- -scale.max
    }
  } else {
    norm.embeddings <- cell.embeddings
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(
    reduction.key, seq_len(length.out = n)
  )
  rownames(x = norm.embeddings) <- colnames(x = object)
  colnames(x = norm.embeddings) <- paste0(
    reduction.key, seq_len(length.out = n)
  )
  reduction.data <- CreateDimReducObject(
    embeddings = norm.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = components
  )
  return(reduction.data)
}

# from SeuratObject (not exported)
#' @importFrom SeuratObject VariableFeatures LayerData
#' @importFrom stats var
PrepDR5 <- function(
  object,
  features = NULL,
  layer = "scale.data",
  verbose = TRUE
) {
  layer <- layer[1L]
  olayer <- layer
  layer <- Layers(object = object, search = layer)
  if (is.null(layer)) {
    stop(
      paste0(
        "No layer matching pattern '",
        olayer,
        "' not found. Please run ScaleData and retry"
      )
    )
  }
  data.use <- LayerData(object = object, layer = layer)
  features <- features %||% VariableFeatures(object = object)
  if (!length(x = features)) {
    stop(
      "No variable features, run FindVariableFeatures()",
      " or provide a vector of features",
      call. = FALSE
    )
  }
  if (is(data.use, "IterableMatrix")) {
    features.var <- BPCells::matrix_stats(
      matrix = data.use, row_stats = "variance"
    )$row_stats["variance", ]
  } else {
    features.var <- sparseMatrixStats::rowVars(x = data.use)
  }
  features.keep <- features[features.var > 0]
  if (!length(x = features.keep)) {
    stop("None of the requested features have any variance", call. = FALSE)
  } else if (length(x = features.keep) < length(x = features)) {
    exclude <- setdiff(x = features, y = features.keep)
    if (isTRUE(x = verbose)) {
      warning(
        "The following ",
        length(x = exclude),
        " features requested have zero variance;",
        " running reduction without them: ",
        paste(exclude, collapse = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  features.use <- features[features %in% rownames(data.use)]
  if (!isTRUE(all.equal(features, features.use))) {
    missing_features <- setdiff(features, features.use)
    if (length(missing_features) > 0) {
      warning_message <- paste("The following features were not available: ",
        paste(missing_features, collapse = ", "),
        ".",
        sep = ""
      )
      warning(warning_message, immediate. = TRUE)
    }
  }
  data.use <- data.use[features.use, ]
  return(data.use)
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @method RunSVD Assay5
#' @examples
#' RunSVD(atac_small[["peaks"]], features = rownames(atac_small))
RunSVD.Assay5 <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  pca = FALSE,
  n = 50,
  reduction.key = ifelse(pca, "PCA_", "LSI_"),
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  data.use <- PrepDR5(
    object = object,
    features = features,
    layer = layer,
    verbose = verbose
  )
  reduction.data <- RunSVD(
    object = data.use,
    assay = assay,
    features = features,
    pca = pca,
    n = n,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @method RunSVD Assay
RunSVD.Assay <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  pca = FALSE,
  n = 50,
  reduction.key = ifelse(pca, "PCA_", "LSI_"),
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  RunSVD.Assay5(
    object = object,
    assay = assay,
    features = features,
    layer = layer,
    n = n,
    pca = pca,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @method RunSVD StdAssay
#' @examples
#' RunSVD(atac_small[["peaks"]], features = rownames(atac_small))
RunSVD.StdAssay <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  pca = FALSE,
  n = 50,
  reduction.key = ifelse(pca, "PCA_", "LSI_"),
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  RunSVD.Assay5(
    object = object,
    assay = assay,
    features = features,
    layer = layer,
    n = n,
    pca = pca,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
}

#' @param reduction.name Name for stored dimension reduction object.
#' @param layer Name of layer to use.
#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @importFrom SeuratObject DefaultAssay
#' @examples
#' RunSVD(atac_small, features = rownames(atac_small))
#' @method RunSVD Seurat
RunSVD.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  layer = "data",
  n = 50,
  pca = FALSE,
  reduction.key = ifelse(pca, "PCA_", "LSI_"),
  reduction.name = ifelse(pca, "pca", "lsi"),
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- object[[assay]]
  reduction.data <- RunSVD(
    object = assay.data,
    assay = assay,
    features = features,
    layer = layer,
    n = n,
    pca = pca,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}

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
  if (any(x > 1) | any(y > 1)) {
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
#' @param irlba.work work parameter for \code{\link[irlba]{irlba}}.
#' Working subspace dimension, larger values can speed convergence at the
#' cost of more memory use.
#' @param tol Tolerance (tol) parameter for \code{\link[irlba]{irlba}}. Larger
#' values speed up convergence due to greater amount of allowed error.
#' @param verbose Print messages
#'
#' @importFrom irlba irlba
#' @importFrom stats sd
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
  scale.embeddings = TRUE,
  reduction.key = "LSI_",
  scale.max = NULL,
  verbose = TRUE,
  irlba.work = n * 3,
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
  if (verbose) {
    message("Running SVD")
  }
  
  if (inherits(x = object, what = 'matrix')) {
    svd.function <- irlba
  } else if (inherits(x = object, what = 'sparseMatrix')) {
    svd.function <- irlba
  } else if (inherits(x = object, what = 'IterableMatrix')) {
    svd.function <- function(A, nv, ...) BPCells::svds(A=A, k = nv)
  } else {
    stop("Unknown matrix format")
  }
  
  components <- svd.function(A = t(x = object), nv = n, work = irlba.work, tol = tol)
  feature.loadings <- components$v
  sdev <- components$d / sqrt(x = max(1, nrow(x = object) - 1))
  cell.embeddings <- components$u
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
#' @importFrom stats var
PrepDR5 <- function(object, features = NULL, layer = 'scale.data', verbose = TRUE) {
  layer <- layer[1L]
  olayer <- layer
  layer <- Layers(object = object, search = layer)
  if (is.null(layer)) {
    stop(paste0("No layer matching pattern '", olayer, "' not found. Please run ScaleData and retry"))
  }
  data.use <- LayerData(object = object, layer = layer)
  features <- features %||% VariableFeatures(object = object)
  if (!length(x = features)) {
    stop("No variable features, run FindVariableFeatures() or provide a vector of features", call. = FALSE)
  }
  if (is(data.use, "IterableMatrix")) {
    features.var <- BPCells::matrix_stats(matrix=data.use, row_stats="variance")$row_stats["variance",]
  } else {
    features.var <- apply(X = data.use, MARGIN = 1L, FUN = var)
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
        " features requested have zero variance; running reduction without them: ",
        paste(exclude, collapse = ', '),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  features.use <- features[features %in% rownames(data.use)]
  if(!isTRUE(all.equal(features, features.use))) {
    missing_features <- setdiff(features, features.use)
    if(length(missing_features) > 0) {
      warning_message <- paste("The following features were not available: ",
                               paste(missing_features, collapse = ", "),
                               ".", sep = "")
      warning(warning_message, immediate. = TRUE)
    }
  }
  data.use <- data.use[features.use, ]
  return(data.use)
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname RunSVD
#' @importFrom SeuratObject VariableFeatures GetAssayData
#' @export
#' @concept dimension_reduction
#' @method RunSVD Assay5
#' @examples
#' \dontrun{
#' RunSVD(atac_small[['peaks']])
#' }
RunSVD.Assay5 <- function(
  object,
  assay = NULL,
  layer = "data",
  features = NULL,
  n = 50,
  reduction.key = "LSI_",
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
  # features <- SetIfNull(x = features, y = VariableFeatures(object = object))
  # data.use <- GetAssayData(
  #   object = object,
  #   layer = "data"
  # )[features, ]
  reduction.data <- RunSVD(
    object = data.use,
    assay = assay,
    features = features,
    n = n,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname RunSVD
#' @importFrom SeuratObject VariableFeatures GetAssayData
#' @export
#' @concept dimension_reduction
#' @method RunSVD StdAssay
#' @examples
#' \dontrun{
#' RunSVD(atac_small[['peaks']])
#' }
RunSVD.StdAssay <- function(
    object,
    assay = NULL,
    layer = "data",
    features = NULL,
    n = 50,
    reduction.key = "LSI_",
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
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
}

#' @param reduction.name Name for stored dimension reduction object.
#' Default 'svd'
#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @examples
#' \dontrun{
#' RunSVD(atac_small)
#' }
#' @method RunSVD Seurat
RunSVD.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  layer = "data",
  n = 50,
  reduction.key = "LSI_",
  reduction.name = "lsi",
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  assay.data <- object[[assay]]
  reduction.data <- RunSVD(
    object = assay.data,
    assay = assay,
    features = features,
    layer = layer,
    n = n,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}

#' @param assay Which assay to use. If NULL, use the default assay
#' @param n Number of singular values to compute
#' @param weight.by.var Weight PCs by variance explained
#' @param reduction.key Key for dimension reduction object
#' @param irlba.work work parameter for \code{\link[irlba]{irlba}}.
#' Working subspace dimension, larger values can speed convergence at the
#' cost of more memory use.
#' @param tol Tolerance (tol) parameter for \code{\link[irlba]{irlba}}. Larger
#' values speed up convergence due to greater amount of allowed error.
#' @param verbose Print messages
#'
#' @importFrom irlba irlba
#' @importFrom stats sd
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom sparseMatrixStats rowVars
#' @importMethodsFrom Matrix t
#'
#' @rdname SparsePCA
#' @export
#' @concept dimension_reduction
#' @examples
#' x <- matrix(data = rnorm(100), ncol = 10)
#' SparsePCA(x)
SparsePCA.default <- function(
    object,
    assay = NULL,
    n = 50,
    weight.by.var = TRUE,
    reduction.key = "PCA_",
    irlba.work = n * 3,
    tol = 1e-05,
    verbose = TRUE,
    ...
) {
  if (is.null(x = rownames(x = object))) {
    rownames(x = object) <- seq_len(length.out = nrow(x = object))
  }
  if (is.null(x = colnames(x = object))) {
    colnames(x = object) <- seq_len(length.out = ncol(x = object))
  }
  n <- min(n, (ncol(x = object) - 1))
  if (verbose) {
    message("Running PCA")
  }
  
  d_rowmeans <- rowMeans(x = object)
  d_sd <- sqrt(x = sparseMatrixStats::rowVars(x = object))
  nz_var <- d_sd > 0
  if (verbose) {
    if (sum(nz_var) != nrow(x = object)) {
      message("Retaining ", sum(nz_var), " features with non-zero variance")
    }
  }
  object <- object[nz_var, ]
  d_rowmeans <- d_rowmeans[nz_var]
  d_sd <- d_sd[nz_var]
  pcs <- irlba::irlba(
    A = t(x = object),
    scale = d_sd,
    center = d_rowmeans,
    nv = n,
    work = irlba.work,
    tol = tol
  )
  if (weight.by.var) {
    emb <- pcs$u %*% diag(pcs$d)
  }
  loadings <- pcs$v
  rownames(x = loadings) <- rownames(x = object)
  colnames(x = loadings) <- paste0(reduction.key, 1:n)
  rownames(x = emb) <- colnames(x = object)
  colnames(x = emb) <- colnames(x = loadings)
  sdev <- pcs$d/sqrt(x = max(1, ncol(x = object) - 1))
  dr <- CreateDimReducObject(
    embeddings = emb,
    loadings = loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key
  )
  return(dr)
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname SparsePCA
#' @export
#' @concept dimension_reduction
#' @method SparsePCA StdAssay
#' @examples
#' \dontrun{
#' SparsePCA(atac_small[['peaks']])
#' }
SparsePCA.StdAssay <- function(
    object,
    assay = NULL,
    features = NULL,
    n = 50,
    weight.by.var = TRUE,
    irlba.work = n * 3,
    tol = 1e-05,
    reduction.key = "PCA_",
    verbose = TRUE,
    ...
) {
  SparsePCA.Assay(
    object = object,
    assay = assay,
    features = features,
    n = n,
    weight.by.var = weight.by.var,
    irlba.work = irlba.work,
    tol = tol,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname SparsePCA
#' @importFrom SeuratObject VariableFeatures GetAssayData
#' @export
#' @concept dimension_reduction
#' @method SparsePCA Assay
#' @examples
#' \dontrun{
#' SparsePCA(atac_small[['peaks']])
#' }
SparsePCA.Assay <- function(
    object,
    assay = NULL,
    features = NULL,
    n = 50,
    weight.by.var = TRUE,
    irlba.work = n * 3,
    tol = 1e-05,
    reduction.key = "PCA_",
    verbose = TRUE,
    ...
) {
  features <- SetIfNull(x = features, y = VariableFeatures(object = object))
  data.use <- GetAssayData(
    object = object,
    layer = "data"
  )[features, ]
  reduction.data <- SparsePCA(
    object = data.use,
    assay = assay,
    features = features,
    n = n,
    weight.by.var = weight.by.var,
    irlba.work = irlba.work,
    tol = tol,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  return(reduction.data)
}

#' @param reduction.name Name for stored dimension reduction object.
#' Default 'pca'
#' @rdname SparsePCA
#' @export
#' @concept dimension_reduction
#' @examples
#' \dontrun{
#' SparsePCA(atac_small)
#' }
#' @method SparsePCA Seurat
SparsePCA.Seurat <- function(
    object,
    assay = NULL,
    features = NULL,
    n = 50,
    weight.by.var = TRUE,
    irlba.work = n * 3,
    tol = 1e-05,
    reduction.key = "PCA_",
    reduction.name = "pca",
    verbose = TRUE,
    ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  assay.data <- object[[assay]]
  reduction.data <- SparsePCA(
    object = assay.data,
    assay = assay,
    features = features,
    n = n,
    weight.by.var = weight.by.var,
    irlba.work = irlba.work,
    tol = tol,
    reduction.key = reduction.key,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}
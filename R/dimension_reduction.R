#' @include generics.R
#'
NULL

#' Calculate the Jaccard index between two matrices
#'
#' Finds the Jaccard similarity between rows of the two matricies. Note that the matrices must be binary,
#' and any rows with zero total counts will result in an NaN entry that could cause problems in downstream analyses.
#'
#' This will calculate the raw Jaccard index, without normalizing for the expected similarity
#' between cells due to differences in sequencing depth.
#'
#' @param x The first matrix
#' @param y The second matrix
#'
#' @importFrom Matrix tcrossprod rowSums
#' @return Returns a matrix
#'
#' @export
#' @examples
#' x <- matrix(data = sample(c(0, 1), size = 25, replace = TRUE), ncol = 5)
#' Jaccard(x = x, y = x)
Jaccard <- function(x, y) {
  if (any(x > 1) | any(y > 1)) {
    warning("Matrices contain values greater than 1. Please binarize matrices before running Jaccard")
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
#' @param standardize.embeddings Scale and center cell embeddings
#' @param scale.max Clipping value for cell embeddings. Default (NULL) is no clipping.
#' @param seed.use Set a random seed. By default, no seed is set.
#' @param verbose Print messages
#'
#' @importFrom irlba irlba
#' @importFrom stats sd
#' @importFrom Seurat CreateDimReducObject
#'
#' @rdname RunSVD
#' @export
#' @examples
#' x <- matrix(data = rnorm(100), ncol = 10)
#' RunSVD(x)
RunSVD.default <- function(
  object,
  assay = NULL,
  n = 50,
  reduction.key = 'SVD_',
  standardize.embeddings = TRUE,
  scale.max = NULL,
  seed.use = NULL,
  verbose = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  n <- min(n, (ncol(x = object) - 1))
  if (verbose) {
    message("Running SVD")
  }
  components <- irlba(A = t(object), nv = n)
  feature.loadings <- components$v
  sdev <- components$d / sqrt(x = max(1, nrow(x = object) - 1))
  cell.embeddings <- components$u
  if (standardize.embeddings) {
    if (verbose) {
      message('Scaling cell embeddings')
    }
    embed.mean <- apply(X = cell.embeddings, MARGIN = 1, FUN = mean)
    embed.sd <- apply(X = cell.embeddings, MARGIN = 1, FUN = sd)
    cell.embeddings <- (cell.embeddings - embed.mean) / embed.sd
  }
  if (!is.null(x = scale.max)) {
    cell.embeddings[cell.embeddings > scale.max] <- scale.max
    cell.embeddings[cell.embeddings < -scale.max] <- -scale.max
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:n)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- paste0(reduction.key, 1:n)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key
  )
  return(reduction.data)
}

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname RunSVD
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
#' @method RunSVD Assay
RunSVD.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = 'SVD_',
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  features <- features %||% VariableFeatures(object = object)
  data.use <- GetAssayData(
    object = object,
    slot = 'data'
  )[features, ]
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

#' @param reduction.name Name for stored dimension reduction object. Default 'lsi'
#'
#' @rdname RunSVD
#'
#' @export
#' @method RunSVD Seurat
RunSVD.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = 'SVD_',
  reduction.name = 'svd',
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  reduction.data <- RunSVD(
    object = assay.data,
    assay = assay,
    features = features,
    n = n,
    reduction.key = reduction.key,
    scale.max = scale.max,
    verbose = verbose,
    ...
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}

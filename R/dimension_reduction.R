#' @include generics.R
#'
NULL

#' Calculate the Jaccard index between two matrices
#'
#' Finds the Jaccard similarity between rows of the two matricies. Note that
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
#' @param verbose Print messages
#'
#' @importFrom irlba irlba
#' @importFrom stats sd
#' @importFrom Seurat CreateDimReducObject
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
  ...
) {
  n <- min(n, (ncol(x = object) - 1))
  if (verbose) {
    message("Running SVD")
  }
  components <- irlba(A = t(x = object), nv = n, work = irlba.work)
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

#' @param features Which features to use. If NULL, use variable features
#'
#' @rdname RunSVD
#' @importFrom Seurat VariableFeatures GetAssayData
#' @export
#' @concept dimension_reduction
#' @method RunSVD Assay
#' @examples
#' RunSVD(atac_small[['peaks']])
RunSVD.Assay <- function(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "LSI_",
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  features <- SetIfNull(x = features, y = VariableFeatures(object = object))
  data.use <- GetAssayData(
    object = object,
    slot = "data"
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

#' @param reduction.name Name for stored dimension reduction object.
#' Default 'svd'
#' @rdname RunSVD
#' @export
#' @concept dimension_reduction
#' @examples
#' RunSVD(atac_small)
#' @method RunSVD Seurat
RunSVD.Seurat <- function(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "LSI_",
  reduction.name = "lsi",
  scale.max = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
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

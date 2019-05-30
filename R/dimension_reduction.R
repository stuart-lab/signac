#' @param assay Which assay to use. If NULL, use the default assay
#' @param n Number of singular values to compute
#' @param reduction.key Key for dimension reduction object
#' @param scale.max Clipping value for cell embeddings. Default (NULL) is no clipping.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param verbose Print messages
#'
#' @importFrom irlba irlba
#'
#' @rdname RunSVD
#' @export
RunSVD.default <- function(
  object,
  assay = NULL,
  n = 50,
  reduction.key = 'SVD_',
  scale.max = NULL,
  seed.use = 42,
  verbose = TRUE,
  ...
) {
  if (!is.null(seed.use)) {
    set.seed(seed = seed.use)
  }
  n <- min(n, ncol(x = object) - 1)
  if (verbose) {
    message("Running SVD")
  }
  components <- irlba(A = t(object), nv = n)
  feature.loadings <- components$v
  sdev <- components$d / sqrt(max(1, nrow(x = object) - 1))
  cell.embeddings <- components$u
  if (verbose) {
    message('Scaling cell embeddings')
  }
  embed.mean <- apply(X = cell.embeddings, MARGIN = 1, FUN = mean)
  embed.sd <- apply(X = cell.embeddings, MARGIN = 1, FUN = sd)
  norm.embeddings <- (cell.embeddings - embed.mean) / embed.sd
  if (!is.null(x = scale.max)) {
    norm.embeddings[norm.embeddings > scale.max] <- scale.max
    norm.embeddings[norm.embeddings < -scale.max] <- -scale.max
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:n)
  rownames(x = norm.embeddings) <- colnames(x = object)
  colnames(x = norm.embeddings) <- paste0(reduction.key, 1:n)
  reduction.data <- CreateDimReducObject(
    embeddings = norm.embeddings,
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
  features <- features %||% VariableFeatures(object)
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
  assay <- assay %||% DefaultAssay(object)
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

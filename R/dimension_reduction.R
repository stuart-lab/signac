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

#' @param verbose Display messages
#' @param graph.name Name of the neighbor graph to use ('nn' or 'snn')
#' @rdname RunMotifTSNE
#' @method RunMotifTSNE Motif
#' @importFrom Seurat RunTSNE
#' @export
RunMotifTSNE.Motif <- function(
  object,
  verbose = TRUE,
  graph.name = 'nn',
  ...
) {
  neighbor.graph <- GetMotifData(object = object, slot = 'neighbors')
  if (!(graph.name %in% names(x = neighbor.graph))) {
    stop("Requested neighbor graph is not present")
  }
  tsne.obj <- RunTSNE(
    object = as.matrix(x = neighbor.graph[[graph.name]]),
    is_distance = TRUE,
    assay = 'Motif',
    ...
  )
  reductions <- GetMotifData(object = object, slot = 'reductions')
  reductions$tSNE <- tsne.obj
  object <- SetMotifData(object = object, slot = 'reductions', new.data = reductions)
  return(object)
}

#' @rdname RunMotifTSNE
#' @method RunMotifTSNE Assay
#' @export
RunMotifTSNE.Assay <- function(
  object,
  verbose = TRUE,
  ...
) {
  motif.obj <- GetMotifObject(object = object)
  motif.obj <- RunMotifTSNE(object = motif.obj, verbose = verbose, ...)
  object <- AddMotifObject(object = object, motif.object = motif.obj, verbose = FALSE)
  return(object)
}

#' @param assay Name of assay to use
#' @rdname RunMotifTSNE
#' @method RunMotifTSNE Seurat
#' @export
RunMotifTSNE.Seurat <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- RunMotifTSNE(
    object = assay.data,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' @param verbose Display messages
#' @param graph.name Name of the neighbor graph to use ('nn' or 'snn')
#' @importFrom Seurat RunUMAP
#' @rdname RunMotifUMAP
#' @method RunMotifUMAP Motif
#' @export
RunMotifUMAP.Motif <- function(
  object,
  verbose = TRUE,
  graph.name = 'nn',
  ...
) {
  neighbor.graph <- GetMotifData(object = object, slot = 'neighbors')
  if (!(graph.name %in% names(x = neighbor.graph))) {
    stop("Requested neighbor graph is not present")
  }
  umap.obj <- RunUMAP(
    object = as.matrix(x = neighbor.graph[[graph.name]]),
    assay = 'Motif',
    ...
  )
  reductions <- GetMotifData(object = object, slot = 'reductions')
  reductions$UMAP <- umap.obj
  object <- SetMotifData(object = object, slot = 'reductions', new.data = reductions)
  return(object)
}

#' @rdname RunMotifUMAP
#' @method RunMotifUMAP Assay
#' @export
RunMotifUMAP.Assay <- function(
  object,
  verbose = TRUE,
  ...
) {
  motif.obj <- GetMotifObject(object = object)
  motif.obj <- RunMotifUMAP(object = motif.obj, verbose = verbose, ...)
  object <- AddMotifObject(object = object, motif.object = motif.obj, verbose = FALSE)
  return(object)
}

#' @param assay Name of assay to use
#' @rdname RunMotifUMAP
#' @method RunMotifUMAP Seurat
#' @export
RunMotifUMAP.Seurat <- function(
  object,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  assay.data <- GetAssay(object = object, assay = assay)
  assay.data <- RunMotifUMAP(
    object = assay.data,
    verbose = verbose,
    ...
  )
  object[[assay]] <- assay.data
  return(object)
}

#' @param assay Which assay to use. If NULL, use the default assay
#' @param n Number of singular values to compute
#' @param reduction.key Key for dimension reduction object
#' @param scale.max Clipping value for cell embeddings. Default (NULL) is no clipping.
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting
#' NULL will not set a seed.
#' @param verbose Print messages
#'
#' @importFrom irlba irlba
#' @importFrom stats sd
#' @importFrom Seurat CreateDimReducObject
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

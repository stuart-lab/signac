#' @importFrom methods setClass
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

## Class definitions

#' The Motif class
#'
#' The Motif class stores DNA motif information
#'
#' @slot data A sparse, binary, motif x feature matrix. Rows correspond to motif IDs, columns correspond to genomic features (peaks or bins).
#' Entries in the matrix should be 1 if the genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A list of position weight matrices for each motif
#' @slot meta.data A dataframe for storage of additional information related to each motif. This could include the names of proteins that bind the motif.
#'
#' @name Motif-class
#' @rdname Motif-class
#' @exportClass Motif
#'
Motif <- setClass(
  Class = 'Motif',
  slots = list(
    data = 'dgCMatrix',
    pwm = 'ANY',
    meta.data = 'ANY'
  )
)

## Functions

#' @param motif.object An object of class Motif
#' @rdname AddMotifObject
#' @method AddMotifObject Assay
#' @export
AddMotifObject.Assay <- function(
  object,
  motif.object
) {
  misc.data <- slot(object = object, name = 'misc') %||% list()
  if ('motif' %in% names(x = misc.data)) {
    warning('Overwriting existing motif object in assay')
  }
  if (!all(rownames(x = object) == colnames(x = motif.object))) {
    keep.features <- intersect(x = colnames(x = motif.object), y = rownames(x = object))
    if (length(x = keep.features) == 0) {
      stop('No features in common between the Assay and Motif objects')
    } else {
      warning('Features do not match in Assay and Motif object. Subsetting the Motif object.')
      motif.object <- motif.obj[, keep.features]
    }
  }
  misc.data[['motif']] <- motif.object
  slot(object = object, name = 'misc') <- misc.data
  return(object)
}

#' @param assay Name of assay to store motif object in
#' @rdname AddMotifObject
#' @method AddMotifObject Seurat
#' @export
AddMotifObject.Seurat <- function(
  object,
  motif.object,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- AddMotifObject(
    object = GetAssay(
      object = object,
      assay = assay
    ),
    motif.object = motif.object)
}

#' CreateMotifObject
#'
#' Create an object of class \code{Motif}
#'
#' @param data A
CreateMotifObject <- function(
  data,
  pwm = NULL,
  meta.data = NULL
) {
  if (!(class(x = data) %in% c('matrix', 'dgCMatrix'))) {
    stop('Data must be matrix or sparse matrix class. Supplied ', class(data))
  }
  if (class(x = data) == 'matrix') {
    data <- as(Class = 'dgCMatrix', object = data)
  }
  if (!is.null(x = pwm)) {
    if (!all(names(x = pwm) == rownames(x = data))) {
      stop('Motif names in data matrix and PWM list are inconsistent')
    }
  } else {
    pwm <- list()
  }
  if (!is.null(x = meta.data)) {
    if (!all(rownames(x = meta.data) == rownames(x = data))) {
      stop('Motif names in data matrix and metadata are inconsistent')
    }
  } else {
    meta.data <- data.frame()
  }
  motif.obj <- new(
    Class = 'Motif',
    data = data,
    pwm = pwm,
    meta.data = meta.data
  )
  return(motif.obj)
}

#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
GetMotifData.Motif <- function(object, slot = 'data', ...) {
  return(slot(object = object, name = slot))
}

#' @rdname GetMotifData
#' @method GetMotifData Assay
#' @export
GetMotifData.Assay <- function(object, slot = 'data', ...) {
  misc.data <- slot(object = object, name = 'misc')
  if ('motif' %in% names(x = misc.data)) {
    return(GetMotifData(object = misc.data[['motif']], slot = slot, ...))
  } else {
    stop('Motif object not present in assay')
  }
}

#' @param assay Which assay to use. Default is the current active assay
#' @rdname GetMotifData
#' @method GetMotifData Seurat
#' @export
GetMotifData.Seurat <- function(object, assay = NULL, slot = 'data', ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(GetMotifData(
    object = GetAssay(object = object, assay = assay)
  ))
}

#' @param new.data New data to add
#' @rdname SetMotifData
#' @method SetMotifData Motif
#' @export
SetMotifData.Motif <- function(object, slot, new.data, ...) {
  slots.use <- c('data', 'pwm', 'meta.data')
  if (!(slot %in% slots.use)) {
    stop('slot must be one of ', paste(slots.use, collapse = ', '), call. = FALSE)
  }
  if (slot == 'data') {
    if (class(x = new.data) == 'matrix') {
      new.data <- as(Class = 'dgCMatrix', object = new.data)
    }
  }
  # TODO check that new data is compatible with existing slots
  # rownames of data must match rownames of meta.data and names of pwm, if not empty
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param data motif matrix to add. Should be matrix or sparse matrix class
#' @rdname SetMotifData
#' @export
#' @method SetMotifData Assay
#' @import Matrix
SetMotifData.Assay <- function(object, slot, new.data, ...) {
  if (!(class(x = new.data) %in% c('matrix', 'dgCMatrix'))) {
    stop('Data must be matrix or sparse matrix class. Supplied ', class(x = new.data))
  }
  if (!all(rownames(x = object) == colnames(x = new.data))) {
    stop('Features do not match existing assay data. Column names in motif matrix should match row names in assay data')
  }
  if (class(x = new.data) == 'matrix') {
    new.data <- as(Class = 'dgCMatrix', object = new.data)
  }
  misc.data <- slot(object = object, name = 'misc') %||% list()
  if (class(misc.data) != 'list') {
    stop('misc slot already occupied and would be overwritten.
         This can be avoided by converting the data in misc to a list,
         so that additional data can be added')
  }
  if (!('motif' %in% names(x = misc.data))) {
    stop('Motif object not present in assay')
  }
  misc.data[['motif']] <- SetMotifData(object = misc.data[['motif']], slot = slot, new.data = new.data)
  slot(object = object, name = 'misc') <- misc.data
  return(object)
}

#' @param assay Name of assay whose data should be set
#' @rdname SetMotifData
#' @export
#' @method SetMotifData Seurat
SetMotifData.Seurat <- function(object, data, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- SetMotifData(object = object[[assay]], data = data, ...)
  return(object)
}

#' @method subset Motif
#' @aliases subset
#' @rdname subset.Motif
#' @seealso \code{\link[base]{subset}}
#' @export
subset.Motif <- function(x, features = NULL, motifs = NULL, ...) {
  features <- features %||% colnames(x = x)
  motifs <- motifs %||% rownames(x = x)
  new.data <- GetMotifData(object = x, slot = 'data')[motifs, features]
  new.pwm <- GetMotifData(object = x, slot = 'pwm')[motifs]
  new.meta <- GetMotifData(object = x, slot = 'meta.data')[motifs, ]
  new.motif <- new(
    Class = 'Motif',
    data = new.data,
    pwm = new.pwm,
    meta.data = new.meta
  )
  return(new.motif)
}

#' @inheritParams subset.Motif
#'
#' @rdname subset.Motif
#' @export
#' @method [ Motif
#'
"[.Motif" <- function(x, i, j, ...) {
  if (missing(x = i) && missing(x = j)) {
    return(x)
  }
  if (missing(x = i)) {
    i <- NULL
  } else if (missing(x = j)) {
    j <- colnames(x = x)
  }
  if (is.logical(x = i)) {
    if (length(i) != nrow(x = x)) {
      stop("Incorrect number of logical values provided to subset features")
    }
    i <- rownames(x = x)[i]
  }
  if (is.logical(x = j)) {
    if (length(j) != ncol(x = x)) {
      stop("Incorrect number of logical values provided to subset cells")
    }
    j <- colnames(x = x)[j]
  }
  if (is.numeric(x = i)) {
    i <- rownames(x = x)[i]
  }
  if (is.numeric(x = j)) {
    j <- colnames(x = x)[j]
  }
  return(subset.Motif(x = x, features = j, motifs = i, ...))
}

## S4 methods

setMethod(
  f = 'show',
  signature = 'Motif',
  definition = function(object) {
    cat('A Motif object containing', nrow(x = slot(object = object, name = "data")),
        "motifs in", ncol(x = slot(object = object, name = "data")), "regions\n")
  }
)

#' @method dimnames Motif
#' @export
dimnames.Motif <- function(x) {
  return(dimnames(x = GetMotifData(object = x)))
}

#' @method dim Motif
#' @export
dim.Motif <- function(x) {
  return(dim(x = GetMotifData(object = x)))
}

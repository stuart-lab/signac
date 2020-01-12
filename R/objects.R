#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod is slot slot<- new as slotNames
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

## Class definitions

#' The Motif class
#'
#' The Motif class stores DNA motif information
#'
#' @slot data A sparse, binary, feature x motif matrix. Columns
#' correspond to motif IDs, rows correspond to genomic features
#' (peaks or bins). Entries in the matrix should be 1 if the
#' genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A named list of position weight matrices
#' @slot motif.names A list containing the name of each motif
#' @slot meta.data A dataframe for storage of additional
#' information related to each motif. This could include the
#' names of proteins that bind the motif.
#'
#' @name Motif-class
#' @rdname Motif-class
#' @exportClass Motif
#'
Motif <- setClass(
  Class = 'Motif',
  slots = list(
    data = 'dgCMatrix',
    pwm = 'list',
    motif.names = 'list',
    meta.data = 'data.frame'
  )
)

## Functions

#' @param motif.object An object of class Motif
#' @param verbose Display messages
#' @param ... Additional arguments
#' @rdname AddMotifObject
#' @method AddMotifObject Assay
#' @export
#' @examples
#' obj <- GetMotifObject(atac_small[['peaks']])
#' atac_small[['peaks']] <- AddMotifObject(object = atac_small[['peaks']], motif.object = obj)
AddMotifObject.Assay <- function(
  object,
  motif.object,
  verbose = TRUE,
  ...
) {
  misc.data <- SetIfNull(x = slot(object = object, name = 'misc'), y = list())
  if ('motif' %in% names(x = misc.data) & verbose) {
    warning('Overwriting existing motif object in assay')
  }
  if (!all(rownames(x = object) == rownames(x = motif.object))) {
    keep.features <- intersect(x = rownames(x = motif.object), y = rownames(x = object))
    if (length(x = keep.features) == 0) {
      stop('No features in common between the Assay and Motif objects')
    } else {
      warning('Features do not match in Assay and Motif object. Subsetting the Motif object.')
      motif.object <- motif.object[keep.features, ]
    }
  }
  misc.data[['motif']] <- motif.object
  slot(object = object, name = 'misc') <- misc.data
  return(object)
}

#' @param assay Name of assay to store motif object in
#' @rdname AddMotifObject
#' @importFrom Seurat DefaultAssay
#' @method AddMotifObject Seurat
#' @export
#' @examples
#' obj <- GetMotifObject(object = atac_small)
#' atac_small[['peaks']] <- AddMotifObject(object = atac_small, motif.object = obj)
AddMotifObject.Seurat <- function(
  object,
  motif.object,
  assay = NULL,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
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
#' @param data A motif x region matrix
#' @param pwm A named list of position weight matrices or position frequency
#' matrices matching the motif names in \code{data}.
#' Can be of class PFMatrixList.
#' @param motif.names A named list of motif names. List element names
#' must match the names given in \code{pwm}. If NULL, use the names from the list
#' of position weight or position frequency matrices. This can be used to set
#' a alternative common name for the motif. If a PFMatrixList is passed to
#' \code{pwm}, it will pull the motif name from the PFMatrixList.
#' @param meta.data A data.frame containing metadata
#' @export
#' @return Returns a \code{\link{Motif}} object
#' @examples
#' motif.matrix <- matrix(data = sample(c(0,1), size = 100, replace = TRUE), ncol = 5)
#' motif <- CreateMotifObject(data = motif.matrix)
CreateMotifObject <- function(
  data = NULL,
  pwm = NULL,
  motif.names = NULL,
  meta.data = NULL
) {
  data <- SetIfNull(x = data, y = new(Class = 'dgCMatrix'))
  meta.data <- SetIfNull(x = meta.data, y = data.frame())
  if (!(inherits(x = data, what = 'matrix') | inherits(x = data, what = 'dgCMatrix'))) {
    stop('Data must be matrix or sparse matrix class. Supplied ', class(x = data))
  }
  if (inherits(x = data, what = 'matrix')) {
    data <- as(Class = 'dgCMatrix', object = data)
  }
  if ((nrow(x = data) > 0) & (length(x = pwm) > 0)) {
    if (!all(names(x = pwm) == colnames(x = data))) {
      stop('Motif names in data matrix and PWM list are inconsistent')
    }
  }
  if ((nrow(x = data) > 0) & (nrow(x = meta.data) > 0)) {
    if (!all(rownames(x = meta.data) == rownames(x = data))) {
      stop('Motif names in data matrix and metadata are inconsistent')
    }
  }
  if (inherits(x = pwm, what = 'list')) {
    if (is.null(names(x = pwm))) {
      stop("PWM must be a named list")
    }
  }
  if (!is.null(x = motif.names)) {
    if (length(x = motif.names) != length(x = pwm)) {
      stop("Number of motif names supplied does not match the number of motifs")
    }
  }
  if (inherits(x = pwm, what = "PFMatrixList") | inherits(x = pwm, what = "PWMatrixList")) {
    pwm.converted <- lapply(X = as.list(x = pwm), FUN = PFMatrixToList)
    pwm <- lapply(X = pwm.converted, FUN = "[[", 1)
    motif.names <- lapply(X = pwm.converted, FUN = "[[", 2)
  }
  pwm <- SetIfNull(x = pwm, y = list())
  if (is.null(x = motif.names)) {
    motif.names <- as.list(x = names(x = pwm))
    names(motif.names) <- names(x = pwm)
  }
  motif.obj <- new(
    Class = 'Motif',
    data = data,
    pwm = pwm,
    motif.names = motif.names,
    meta.data = meta.data
  )
  return(motif.obj)
}

#' @rdname GetMotifObject
#' @method GetMotifObject Assay
#' @export
#' @examples
#' GetMotifObject(object = atac_small[['peaks']])
GetMotifObject.Assay <- function(object, ...) {
  misc.data <- slot(object = object, name = 'misc')
  if ('motif' %in% names(x = misc.data)) {
    return(misc.data[['motif']])
  } else {
    stop('Motif object not present in assay')
  }
}

#' @param assay Which assay to use. Default is the current active assay
#' @rdname GetMotifObject
#' @importFrom Seurat DefaultAssay GetAssay
#' @method GetMotifObject Seurat
#' @export
#' @examples
#' GetMotifObject(object = atac_small)
GetMotifObject.Seurat <- function(object, assay = NULL, ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  return(GetMotifObject(
    object = GetAssay(object = object, assay = assay),
    ...
  ))
}

#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
#' @examples
#' motif.obj <- GetMotifObject(object = atac_small[['peaks']])
#' GetMotifData(object = motif.obj)
GetMotifData.Motif <- function(object, slot = 'data', ...) {
  return(slot(object = object, name = slot))
}

#' @rdname GetMotifData
#' @method GetMotifData Assay
#' @export
#' @examples
#' GetMotifData(object = atac_small[['peaks']])
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
#' @importFrom Seurat DefaultAssay GetAssay
#' @export
#' @examples
#' GetMotifData(object = atac_small)
GetMotifData.Seurat <- function(object, assay = NULL, slot = 'data', ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  return(GetMotifData(
    object = GetAssay(object = object, assay = assay),
    slot = slot,
    ...
  ))
}

#' @rdname SetMotifData
#' @method SetMotifData Motif
#' @export
#' @examples
#' motif.obj <- GetMotifObject(object = atac_small)
#' SetMotifData(object = motif.obj, slot = 'data', new.data = matrix())
SetMotifData.Motif <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop('slot must be one of ', paste(slotNames(x = object), collapse = ', '), call. = FALSE)
  }
  if (slot == 'data') {
    if (inherits(x = new.data, what = 'matrix')) {
      new.data <- as(Class = 'dgCMatrix', object = new.data)
    }
  }
  # TODO check that new data is compatible with existing slots
  # rownames of data must match rownames of meta.data and names of pwm, if not empty
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param new.data motif matrix to add. Should be matrix or sparse matrix class
#' @param slot Name of slot to use
#' @rdname SetMotifData
#' @export
#' @method SetMotifData Assay
#' @examples
#' SetMotifData(object = atac_small[['peaks']], slot = 'data', new.data = matrix())
SetMotifData.Assay <- function(object, slot, new.data, ...) {
  if (slot == 'data') {
    if (!(inherits(x = new.data, what = 'matrix') | inherits(x = new.data, what = 'dgCMatrix'))) {
      stop('Data must be matrix or sparse matrix class. Supplied ', class(x = new.data))
    }
    if (!all(rownames(x = object) == rownames(x = new.data))) {
      stop('Features do not match existing assay data. Column names in motif matrix should match row names in assay data')
    }
    if (inherits(x = new.data, what = 'matrix')) {
      new.data <- as(Class = 'dgCMatrix', object = new.data)
    }
  }
  misc.data <- SetIfNull(x = slot(object = object, name = 'misc'), y = list())
  if (!inherits(x = misc.data, what = 'list')) {
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
#' @importFrom Seurat DefaultAssay
#' @export
#' @method SetMotifData Seurat
#' @examples
#' motif.matrix <- GetMotifData(object = atac_small)
#' SetMotifData(object = atac_small, assay = 'peaks', slot = 'data', new.data = motif.matrix)
SetMotifData.Seurat <- function(object, assay = NULL, ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  object[[assay]] <- SetMotifData(object = object[[assay]], ...)
  return(object)
}

#' Return a subset of a Motif object
#'
#' @param x A Motif object
#' @param features Which features to retain
#' @param motifs Which motifs to retain
#' @param ... Arguments passed to other methods
#'
#' @aliases subset
#' @rdname subset.Motif
#' @method subset Motif
#'
#' @seealso \code{\link[base]{subset}}
#' @return Returns a subsetted \code{\link{Motif}} object
#' @export
#' @examples
#' motif.obj <- GetMotifObject(object = atac_small)
#' subset(x = motif.obj, features = head(rownames(motif.obj), 10))
subset.Motif <- function(x, features = NULL, motifs = NULL, ...) {
  features <- SetIfNull(x = features, y = rownames(x = x))
  motifs <- SetIfNull(x = motifs, y = colnames(x = x))
  new.data <- GetMotifData(object = x, slot = 'data')[features, motifs]
  new.pwm <- GetMotifData(object = x, slot = 'pwm')[motifs]
  new.names <- GetMotifData(object = x, slot = 'motif.names')[motifs]
  new.meta <- GetMotifData(object = x, slot = 'meta.data')[motifs, ]
  new.motif <- new(
    Class = 'Motif',
    data = new.data,
    pwm = new.pwm,
    motif.names = new.names,
    meta.data = new.meta
  )
  return(new.motif)
}

#' @param i Which columns to retain
#' @param j Which rows to retain
#'
#' @rdname subset.Motif
#' @export
#' @method [ Motif
#' @examples
#' motif.obj <- GetMotifObject(atac_small)
#' motif.obj[1:10,1:5]
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
  return(subset.Motif(x = x, features = i, motifs = j, ...))
}

## S4 methods

setMethod(
  f = 'show',
  signature = 'Motif',
  definition = function(object) {
    cat('A Motif object containing', ncol(x = slot(object = object, name = "data")),
        "motifs in", nrow(x = slot(object = object, name = "data")), "regions\n")
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

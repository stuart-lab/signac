#' @rdname GetMotifData
#' @method GetMotifData Assay
#' @export
GetMotifData.Assay <- function(object, ...) {
  misc.data <- slot(object = object, name = 'misc')
  if ('motif' %in% names(misc.data)) {
    return(misc.data[['motif']])
  } else {
    stop('Motif data not set')
  }
}

#' @param assay Which assay to use. Default is the current active assay
#' @rdname GetMotifData
#' @method GetMotifData Seurat
#' @export
GetMotifData.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(GetMotifData(
    object = GetAssay(object = object, assay = assay)
  ))
}

#' @param data motif matrix to add. Should be matrix or sparse matrix class
#' @rdname SetMotifData
#' @export
#' @method SetMotifData Assay
#' @import Matrix
SetMotifData.Assay <- function(object, data, ...) {
  if (!(class(data) %in% c('matrix', 'dgCMatrix'))) {
    stop('Data must be matrix or sparse matrix class. Supplied ', class(data))
  }
  if (!all(rownames(x = object) == colnames(x = data))) {
    stop('Features do not match existing assay data. Column names in motif matrix should match row names in assay data')
  }
  if (class(data) == 'matrix') {
    data <- as(Class = 'dgCMatrix', object = data)
  }
  misc.data <- slot(object = object, name = 'misc') %||% list()
  if (class(misc.data) != 'list') {
    stop('misc slot already occupied and would be overwritten.
         This can be avoided by converting the data in misc to a list,
         so that additional data can be added')
  }
  misc.data[['motif']] <- data
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

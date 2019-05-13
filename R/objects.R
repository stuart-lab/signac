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
#' @importFrom Matrix dgCMatrix-class
SetMotifData.Assay <- function(object, data, ...) {
  if (!(class(data) %in% c('matrix', 'dgCMatrix'))) {
    stop('Data must be matrix or sparse matrix class. Supplied ', class(data))
  }
  if (!all(rownames(object = object) == colnames(data))) {
    stop('Features do not match existing assay data. Column names in motif matrix should match row names in assay data')
  }
  if (class(data) == 'matrix') {
    data <- as(Class = 'dgCMatrix', object = data)
  }
  misc.data <- slot(object = object, name = 'misc')
  misc.data[['motif']] <- data
  slot(object = object, name = 'misc') <- misc.data
  return(object)
}

#' @param assay Name of assay whose data should be set
#' @rdname SetMotifData
#' @export
#' @motif SetMotifData Seurat
SetMotifData.Seurat <- function(object, data, assay, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- SetAssayData(object = object[[assay]], data = data, ...)
  return(object)
}

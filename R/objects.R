#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod is slot slot<- new as slotNames
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom TFBSTools PWMatrixList PFMatrixList
#'
NULL

## Class definitions

setClassUnion(name = 'AnyPWM', c("list", "PWMatrixList", "PFMatrixList"))
setClassUnion(name = 'AnyMatrix', c("matrix", "dgCMatrix"))

#' The Motif class
#'
#' The Motif class stores DNA motif information
#'
#' @slot data A sparse, binary, feature x motif matrix. Columns correspond to motif IDs, rows correspond to genomic features (peaks or bins).
#' Entries in the matrix should be 1 if the genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A list of position weight matrices for each motif
#' @slot neighbors A list of nearest neighbors graphs
#' @slot reductions A list of \code{\link[Seurat]{DimReduc}} objects
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
    pwm = 'AnyPWM',
    neighbors = 'list',
    reductions = 'list',
    meta.data = 'data.frame'
  )
)

#' The ChromatinAssay class
#'
#' The ChramatinAssay object is an extended \code{\link[Seurat]{Assay}}
#' for the storage and analysis of chromatin-based single-cell data.
#'
#' @slot ranges A \code{\link[GenomicRanges]{GRanges}} object describing the
#' genomic location of features in the object
#' @slot motifs A \code{\link{Motif}} object
#' @slot fragments A character vector containing the path/s to tabix-indexed fragment file/s
#' for the cells in the assay. See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @slot genome Name of the genome used
#' @slot annotation An object containing genomic annotations
#' @slot bias A matrix containing Tn5 integration bias information (frequency of Tn5 integration at
#' different kmers)
#' @slot positionEnrichment A named list of matrices containing positional enrichment scores for Tn5 integration
#' (for example, enrichment at the TSS)
#'
#' @name ChromatinAssay-class
#' @rdname ChromatinAssay-class
#' @exportClass ChromatinAssay
#'
ChromatinAssay <- setClass(
  Class = 'ChromatinAssay',
  contains = 'Assay',
  slots = list(
    'ranges' = 'GRanges',
    'motifs' = 'ANY',
    'fragments' = 'ANY',
    'genome' = 'ANY',
    'annotation' = 'ANY',
    'bias' = 'ANY',
    'positionEnrichment' = 'list'
  )
)

#' @param ranges A GRanges object
#' @param genome Name of genome used
#' @param annotation Genomic annotation
#' @param motifs A Motif object
#' @param fragments Path to fragments file
#' @param sep Charaters used to separate the chromosome, start, and end coordinates
#' in the row names of the data matrix
#'
#' @rdname as.ChromatinAssay
#' @export
#' @method as.ChromatinAssay Assay
#'
as.ChromatinAssay.Assay <- function(
  x,
  ranges = NULL,
  genome = NULL,
  annotation = NULL,
  motifs = NULL,
  fragments = NULL,
  sep = c("-", "-"),
  ...
) {
  new.assay <- as(object = x, Class = 'ChromatinAssay')
  ranges <- ranges %||% StringToGRanges(regions = rownames(x = x), sep = sep)
  new.assay <- SetAssayData(object = new.assay, slot = 'ranges', new.data = ranges)
  if (!is.null(x = fragments)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = 'fragments',
      new.data = fragments
    )
  }
  if (!is.null(x = genome)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = 'genome',
      new.data = genome
    )
  }
  if (!is.null(x = annotation)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = 'annotation',
      new.data = annotation
    )
  }
  if (!is.null(x = motifs)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = 'motifs',
      new.data = motifs
    )
  }
  return(new.assay)
}

setAs(
  from = 'Assay',
  to = 'ChromatinAssay',
  def = function(from) {
    object.list <- sapply(
      X = slotNames(x = from),
      FUN = slot,
      object = from,
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    object.list <- c(
      list(
        'Class' = 'ChromatinAssay'
      ),
      object.list
    )
    return(do.call(what = 'new', args = object.list))
  }
)

#' Create Seurat object with a chromatin assay
#'
#' This is a wrapper for \code{\link[Seurat]{CreateSeuratObject}} to create
#' an object containing a \code{\link{ChromatinAssay}} rather than a
#' standard Seurat \code{\link[Seurat]{Assay}}.
#'
#' @param counts Unnormalized data (raw counts). Should be a sparse matrix.
#' @param assay Name of the assay corresponding to the initial input data. Default is "ATAC".
#' @param genome Name of genome used
#' @param ranges A \code{\link[GenomicRanges]{GRanges}} object containing the genomic
#' position of each row of the counts matrix
#' @param project Sets the project name for the object
#' @param fragments A character vector containing path/s to tabix-indexed fragment file/s
#' for cells in the object
#' @param annotation A \code{\link[GenomicRanges]{GRanges}} object containing
#' genomic annotations for the genome used
#' @param motifs A \code{\link{Motif}} object
#' @param sep Charaters used to separate the chromosome, start, and end coordinates
#' in the row names of the data matrix
#' @param ... Parameters passed to \code{\link[Seurat]{CreateSeuratObject}}
#'
#' @importFrom Seurat CreateSeuratObject
#' @export
CreateSignacObject <- function(
  counts,
  assay = 'ATAC',
  project = 'SignacProject',
  ranges = NULL,
  fragments = NULL,
  annotation = NULL,
  genome = NULL,
  motifs = NULL,
  sep = c("-", "-"),
  ...
) {
  ranges <- ranges %||% StringToGRanges(regions = rownames(x = counts), sep = sep)
  seurat.obj <- CreateSeuratObject(
    counts = counts,
    project = project,
    assay = 'temp',
    ...
  )
  seurat.obj[[assay]] <- as.ChromatinAssay(
    x = seurat.obj[['temp']],
    ranges = ranges,
    fragments = fragments,
    annotation = annotation,
    genome = genome,
    motifs = motifs
  )
  DefaultAssay(object = seurat.obj) <- assay
  seurat.obj[['temp']] <- NULL
  return(seurat.obj)
}

## Functions

#' CreateMotifObject
#'
#' Create an object of class \code{Motif}
#'
#' @param data A motif x region matrix
#' @param pwm A list of position weight matrices or position frequency matrices matching the motif names in \code{data}.
#' Can be of class PFMatrixList
#' @param neighbors Neighbor data
#' @param reductions Dimension reduction data
#' @param meta.data A data.frame containing metadata
#'
#' @export
#' @examples
#' motif.matrix <- matrix(data = sample(c(0,1), size = 100, replace = TRUE), ncol = 5)
#' motif <- CreateMotifObject(data = motif.matrix)
CreateMotifObject <- function(
  data = NULL,
  pwm = NULL,
  neighbors = NULL,
  reductions = NULL,
  meta.data = NULL
) {
  data <- data %||% new(Class = 'dgCMatrix')
  pwm <- pwm %||% list()
  neighbors <- neighbors %||% list()
  reductions <- reductions %||% list()
  meta.data <- meta.data %||% data.frame()
  if (!(class(x = data) %in% c('matrix', 'dgCMatrix'))) {
    stop('Data must be matrix or sparse matrix class. Supplied ', class(x = data))
  }
  if (is(object = data, class2 = 'matrix')) {
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
  motif.obj <- new(
    Class = 'Motif',
    data = data,
    pwm = pwm,
    neighbors = neighbors,
    reductions = reductions,
    meta.data = meta.data
  )
  return(motif.obj)
}

#' @importFrom Seurat GetAssayData
#' @method GetAssayData ChromatinAssay
#' @export
GetAssayData.ChromatinAssay <- function(object, slot = 'data', assay = NULL, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop('slot must be one of ', paste(slotNames(x = object), collapse = ', '), call. = FALSE)
  }
  return(slot(object = object, name = slot))
}

#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
GetMotifData.Motif <- function(object, slot = 'data', ...) {
  return(slot(object = object, name = slot))
}

#' @importFrom Seurat GetAssayData
#' @rdname GetMotifData
#' @method GetMotifData ChromatinAssay
#' @export
GetMotifData.ChromatinAssay <- function(object, slot = 'data', ...) {
  motif.obj <- GetAssayData(object = object, slot = 'motifs')
  if (is.null(x = motif.obj)) {
    stop("Motif object not present in assay")
  } else {
    return(GetMotifData(object = motif.obj, slot = slot, ...))
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
  assay <- assay %||% DefaultAssay(object = object)
  return(GetMotifData(
    object = GetAssay(object = object, assay = assay),
    slot = slot,
    ...
  ))
}

#' @importFrom Seurat SetAssayData
#' @importFrom GenomeInfoDb genome
#' @method SetAssayData ChromatinAssay
#' @export
SetAssayData.ChromatinAssay <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop('slot must be one of ', paste(slotNames(x = object), collapse = ', '), call. = FALSE)
  }
  if (slot == 'counts') {
    stop("Modifying counts slot is not enabled")
  } else if (slot == 'data') {
    if (!(is(object = new.data, class2 = 'AnyMatrix'))) {
      stop("Data must be a matrix or sparseMatrix")
    }
    if (nrow(x = object) != nrow(x = new.data)) {
      stop('Number of rows in provided matrix does not match the number of rows in the object')
    }
    if (ncol(x = object) != ncol(x = new.data)) {
      stop("Number of columns in the provided matrix does not match the number of cells in the object")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == 'genome') {
    if (!is(object = new.data, class2 = 'character')) {
      stop("Genome must be a character class object")
    }
    # TODO check that genome matches the genome for granges and annotation
    slot(object = object, name = slot) <- new.data
  } else if (slot == 'fragments') {
    index.file <- paste0(new.data, ".tbi")
    if (all(file.exists(new.data, index.file))) {
      file <- normalizePath(path = new.data)
      slot(object = object, name = slot) <- new.data
    }
    else {
      stop("Requested file does not exist or is not indexed")
    }
  } else if (slot == 'annotation') {
    if (!is(object = new.data, class2 = 'GRanges')) {
      stop("Must provide a GRanges object")
    }
    current.genome <- genome(x = object)
    annotation.genome <- unique(x = genome(x = new.data))
    if (!is.null(x = current.genome) & !is.na(x = annotation.genome) & (current.genome != annotation.genome)) {
      stop("Annotation genome does not match genome of the object")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == 'bias') {
    if(!is(object = new.data, class2 = 'AnyMatrix')) {
      stop("Bias must be provided as a matrix or sparseMatrix")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == 'positionEnrichment') {
    if(!is(object = new.data, class2 = 'AnyMatrix')) {
      stop("Position enrichment must be provided as a matrix or sparseMatrix")
    }
    args <- list(...)
    if (!('key' %in% names(x = args))) {
      stop("Must supply a key when adding positionEnrichment data")
    } else {
      key <- args$key
    }
    current.pos <- slot(object = object, name = slot)
    current.pos[[key]] <- new.data
    slot(object = object, name = slot) <- current.pos
  } else if (slot == 'ranges') {
    if (!is(object = new.data, class2 = 'GRanges')) {
      stop("Must provide a GRanges object")
    } else if(length(x = new.data) != nrow(x = object)) {
      stop("Number of ranges provided is not equal to the number of features in the assay")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == 'motifs') {
    if (!inherits(x = new.data, what = 'Motif')) {
      stop("Must provide a Motif class object")
    }
    if (!all(rownames(x = object) == rownames(x = new.data))) {
      keep.features <- intersect(x = rownames(x = new.data),
                                 y = rownames(x = object))
      if (length(x = keep.features) == 0) {
        stop("No features in common between the ChromatinAssay and Motif objects")
      }
      else {
        warning("Features do not match in ChromatinAssay and Motif object. Subsetting the Motif object.")
        new.data <- new.data[keep.features, ]
      }
    }
    slot(object = object, name = slot) <- new.data
  }
  return(object)
}

#' @rdname SetMotifData
#' @method SetMotifData Motif
#' @export
SetMotifData.Motif <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop('slot must be one of ', paste(slotNames(x = object), collapse = ', '), call. = FALSE)
  }
  if (slot == 'data') {
    if (is(object = new.data, class2 = 'matrix')) {
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
#' @importFrom Seurat GetAssayData SetAssayData
#' @rdname SetMotifData
#' @export
#' @method SetMotifData ChromatinAssay
#' @import Matrix
SetMotifData.ChromatinAssay <- function(object, slot, new.data, ...) {
  if (slot == 'data') {
    if (!(class(x = new.data) %in% c('matrix', 'dgCMatrix'))) {
      stop('Data must be matrix or sparse matrix class. Supplied ', class(x = new.data))
    }
    if (!all(rownames(x = object) == rownames(x = new.data))) {
      stop('Features do not match existing assay data. Column names in motif matrix should match row names in assay data')
    }
    if (is(object = new.data, class2 = 'matrix')) {
      new.data <- as(Class = 'dgCMatrix', object = new.data)
    }
  }
  motif.obj <- GetAssayData(object = object, slot = 'motifs')
  if (is.null(x = motif.obj)) {
    stop("Motif object not present in assay")
  } else {
    motif.obj <- SetMotifData(object = motif.obj, slot = slot, new.data = new.data)
    object <- SetAssayData(object = object, slot = 'motifs', new.data = motif.obj)
    return(object)
  }
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
  assay <- assay %||% DefaultAssay(object = object)
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
#' @return Returns a subsetted Motif object
#' @export
#'
subset.Motif <- function(x, features = NULL, motifs = NULL, ...) {
  features <- features %||% rownames(x = x)
  motifs <- motifs %||% colnames(x = x)
  new.data <- GetMotifData(object = x, slot = 'data')[features, motifs]
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

# TODO define subset.ChromatinAssay
# needs to subet the genomic ranges and motif object

# TODO define merge.ChromatinAssay
# should be coordinate-aware and genome-aware
# genome must match in order to merge

#' @inheritParams subset.Motif
#' @param i Which columns to retain
#' @param j Which rows to retain
#'
#' @rdname subset.Motif
#' @export
#' @method [ Motif
#' @examples
#' motif.obj <- GetMotifObject(atac_small)
#' motif.obj[1:10,1:10]
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

#' @importFrom GenomeInfoDb genome
setMethod(
  f = 'show',
  signature = 'ChromatinAssay',
  definition = function(object) {
    cat('ChromatinAssay data with', nrow(x = object), 'features for', ncol(x = object), 'cells\n')
    cat('Variable features:', length(x = VariableFeatures(object = object)), '\n')
    cat('Genome:', genome(x = object), '\n')
  }
)

#' @importFrom GenomicRanges granges
setMethod(
  f = "granges",
  signature = "ChromatinAssay",
  definition = function(x, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for ChromatinAssay objects ",
           "does not support the 'use.mcols' argument")
    }
    slot(object = x, name = 'ranges')
  }
)

#' @importFrom GenomicRanges granges
#' @importFrom Seurat DefaultAssay
setMethod(
  f = "granges",
  signature = "Seurat",
  definition = function(x, use.mcols = FALSE, ...) {
    if (!identical(x = use.mcols, y = FALSE)) {
      stop("\"granges\" method for Seurat objects ",
           "does not support the 'use.mcols' argument")
    }
    assay <- DefaultAssay(object = x)
    granges(x = x[[assay]])
  }
)

#' @importFrom GenomeInfoDb genome
setMethod(
  f = "genome",
  signature = "ChromatinAssay",
  definition = function(x) {
    slot(object = x, name = 'genome')
  }
)

#' @importFrom GenomeInfoDb genome
#' @importFrom Seurat DefaultAssay
setMethod(
  f = "genome",
  signature = "Seurat",
  definition = function(x) {
    assay <- DefaultAssay(object = x)
    genome(x = x[[assay]])
  }
)

#' @rdname Annotation
#' @method Annotation ChromatinAssay
#' @export
Annotation.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = 'annotation'))
}

#' @param object A Seurat object or ChromatinAssay object
#' @importFrom Seurat DefaultAssay
#' @rdname Annotation
#' @method Annotation Seurat
#' @export
Annotation.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Annotation(object = object[[assay]]))
}

#' @rdname Motifs
#' @method Motifs ChromatinAssay
#' @export
Motifs.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = 'motifs'))
}

#' @param object A Seurat object
#' @rdname Motifs
#' @importFrom Seurat DefaultAssay
#' @method Motifs Seurat
#' @export
Motifs.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Motifs(object = object[[assay]]))
}

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

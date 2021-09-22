#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass setClassUnion setMethod is slot slot<- new as
#' slotNames
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Signac
NULL

## Class definitions

setClassUnion(name = "AnyMatrix", c("matrix", "dgCMatrix"))

#' The Fragment class
#'
#' The Fragment class is designed to hold information needed for working with
#' fragment files.
#'
#' @slot path Path to the fragment file on disk.
#' See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file, and the name of each element is the
#' corresponding cell barcode as stored in the ChromatinAssay object.
#'
#' @name Fragment-class
#' @rdname Fragment-class
#' @exportClass Fragment
#' @concept fragments
Fragment <- setClass(
  Class = "Fragment",
  slots = list(
    path = "character",
    hash = "character",
    cells = "ANY"
  )
)

#' The Motif class
#'
#' The Motif class is designed to store DNA sequence motif information,
#' including motif PWMs or PFMs, motif positions, and metadata.
#'
#' @slot data A sparse, binary, feature x motif matrix. Columns
#' correspond to motif IDs, rows correspond to genomic features
#' (peaks or bins). Entries in the matrix should be 1 if the
#' genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A named list of position weight matrices
#' @slot motif.names A list containing the name of each motif
#' @slot positions A \code{\link[GenomicRanges]{GRangesList}} object containing
#' exact positions of each motif.
#' @slot meta.data A dataframe for storage of additional
#' information related to each motif. This could include the
#' names of proteins that bind the motif.
#'
#' @name Motif-class
#' @rdname Motif-class
#' @exportClass Motif
#' @concept motifs
Motif <- setClass(
  Class = "Motif",
  slots = list(
    data = "dgCMatrix",
    pwm = "list",
    motif.names = "list",
    positions = "ANY",
    meta.data = "data.frame"
  )
)

#' The ChromatinAssay class
#'
#' The ChromatinAssay object is an extended \code{\link[SeuratObject]{Assay}}
#' for the storage and analysis of single-cell chromatin data.
#'
#' @slot ranges A \code{\link[GenomicRanges]{GRanges}} object describing the
#' genomic location of features in the object
#' @slot motifs A \code{\link{Motif}} object
#' @slot fragments A list of \code{\link{Fragment}} objects.
#' @slot seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome sequence used.
#' @slot annotation A  \code{\link[GenomicRanges]{GRanges}} object containing
#' genomic annotations
#' @slot bias A vector containing Tn5 integration bias information
#' (frequency of Tn5 integration at different kmers)
#' @slot positionEnrichment A named list of matrices containing positional
#' enrichment scores for Tn5 integration (for example, enrichment at the TSS)
#' @slot links A \code{\link[GenomicRanges]{GRanges}} object describing linked
#' genomic positions, such as co-accessible sites or enhancer-gene regulatory
#' relationships. This should be a \code{GRanges} object, where the start and
#' end coordinates are the two linked genomic positions, and must contain a
#' "score" metadata column.
#'
#' @name ChromatinAssay-class
#' @rdname ChromatinAssay-class
#' @exportClass ChromatinAssay
#' @concept assay
ChromatinAssay <- setClass(
  Class = "ChromatinAssay",
  contains = "Assay",
  slots = list(
    "ranges" = "GRanges",
    "motifs" = "ANY",
    "fragments" = "list",
    "seqinfo" = "ANY",
    "annotation" = "ANY",
    "bias" = "ANY",
    "positionEnrichment" = "list",
    "links" = "GRanges"
  )
)

#' Create ChromatinAssay object
#'
#' Create a \code{\link{ChromatinAssay}} object from a count matrix or
#' normalized data matrix. The expected format of the input matrix is features x
#' cells. A set of genomic ranges must be supplied along with the matrix, with
#' the length of the ranges equal to the number of rows in the matrix. If a set
#' of genomic ranges are not supplied, they will be extracted from the
#' row names of the matrix.
#'
#' @param counts Unnormalized data (raw counts)
#' @param data Normalized data; if provided, do not pass counts
#' @param min.cells Include features detected in at least this many cells.
#' Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a lower cutoff.
#' @param max.cells Include features detected in less than this many cells.
#' Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a higher cutoff.
#' This can be useful for chromatin assays where certain artefactual loci
#' accumulate reads in all cells. A percentage cutoff can also be set using
#' 'q' followed by the percentage of cells, for example 'q90' will discard
#' features detected in 90 percent of cells.
#' If NULL (default), do not apply any maximum value.
#' @param min.features Include cells where at least this many features are
#' detected.
#' @param ranges A set of \code{\link[GenomicRanges]{GRanges}} corresponding to
#' the rows of the input matrix
#' @param motifs A Motif object (not required)
#' @param fragments Path to a tabix-indexed fragments file for the data
#' contained in the input matrix. If multiple fragment files are required,
#' you can add additional \code{\link{Fragment}} object to the assay after it is
#' created using the \code{\link{CreateFragmentObject}} and
#' \code{\link{Fragments}} functions. Alternatively, a list of
#' \code{\link{Fragment}} objects can be provided.
#' @param genome A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotation A set of \code{\link[GenomicRanges]{GRanges}} containing
#' annotations for the genome used
#' @param bias A Tn5 integration bias matrix
#' @param positionEnrichment A named list of matrices containing positional
#' signal enrichment information for each cell. Should be a cell x position
#' matrix, centered on an element of interest (for example, TSS sites).
#' @param sep Separators to use for strings encoding genomic coordinates.
#' First element is used to separate the chromosome from the coordinates,
#' second element is used to separate the start from end coordinate. Only
#' used if \code{ranges} is NULL.
#' @param validate.fragments Check that cells in the assay are present in the
#' fragment file.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{CreateFragmentObject}}
#'
#' @importFrom Seurat CreateAssayObject
#' @importFrom Matrix rowSums colSums
#' @importFrom GenomicRanges isDisjoint
#' @concept assay
#'
#' @export
CreateChromatinAssay <- function(
  counts,
  data,
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  ranges = NULL,
  motifs = NULL,
  fragments = NULL,
  genome = NULL,
  annotation = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c("-", "-"),
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
) {
  if (missing(x = counts) && missing(x = data)) {
    stop("Must provide either 'counts' or 'data'")
  } else if (!missing(x = counts) && !missing(x = data)) {
    stop("Either 'counts' or 'data' must be missing; both cannot be provided")
  } else if (!missing(x = counts)) {
    data.use <- counts
  } else {
    data.use <- data
  }
  if (!is.null(x = ranges)) {
    if (length(x = ranges) != nrow(x = data.use)) {
      stop("Length of supplied genomic ranges does not match number
           of rows in matrix")
    }
  } else {
    ranges <- StringToGRanges(regions = rownames(x = data.use), sep = sep)
  }
  if (!isDisjoint(x = ranges)) {
    warning("Overlapping ranges supplied. Ranges should be non-overlapping.")
  }
  if (!is.null(x = annotation) & !inherits(x = annotation, what = "GRanges")) {
    stop("Annotation must be a GRanges object.")
  }
  # remove low-count cells
  ncount.cell <- colSums(x = data.use > 0)
  data.use <- data.use[, ncount.cell > min.features]

  if (ncol(x = data.use) == 0) {
    stop("No cells retained due to minimum feature cutoff supplied")
  }

  ncell.feature <- rowSums(x = data.use > 0)
  if (!is.null(x = max.cells)) {
    if (is(object = max.cells, class2 = "character")) {
      percent.cutoff <- as.numeric(
        x = gsub(pattern = "q", replacement = "", x = max.cells)
      )
      max.cells <- (percent.cutoff / 100) * ncol(x = data.use)
    }
  } else {
    max.cells <- ncol(x = data.use)
  }
  features.keep <- (ncell.feature >= min.cells) & (ncell.feature <= max.cells)
  if (sum(features.keep) == 0) {
    stop("No features retained due to minimum cell cutoff supplied")
  }
  data.use <- data.use[features.keep, ]
  ranges <- ranges[features.keep, ]
  # re-assign row names of matrix so that it's a known granges transformation
  new.rownames <- GRangesToString(grange = ranges, sep = c("-", "-"))
  rownames(x = data.use) <- new.rownames
  if (!missing(x = counts)) {
    seurat.assay <- CreateAssayObject(
      counts = data.use,
      data = data,
      min.cells = -1,
      min.features = -1 # min cell/feature filtering already done
    )
  } else {
    seurat.assay <- CreateAssayObject(
      counts = counts,
      data = data.use,
      min.cells = min.cells,
      min.features = min.features
    )
  }
  if (inherits(x = fragments, what = "list")) {
    # check each object in the list is a fragment object
    # fragment list usually supplied when doing object merge,
    # so don't validate cells here, we can assume that was done in
    # individual object creation
    obj.class <- sapply(
      X = fragments, FUN = function(x) inherits(x = x, what = "Fragment")
    )
    if (!all(obj.class)) {
      stop("All objects in fragments list must be Fragment-class objects")
    }
    frags <- lapply(
      X = fragments,
      FUN = AssignFragCellnames,
      cellnames = colnames(x = seurat.assay)
    )
   } else if (inherits(x = fragments, what = "Fragment")) {
    # single Fragment object supplied
    frags <- AssignFragCellnames(
      fragments = fragments, cellnames = colnames(x = seurat.assay)
    )
  } else {
    # path to fragment file supplied, create fragment object
    frags <- list()
    if (!is.null(x = fragments)) {
      if (nchar(x = fragments) > 0) {
        cells <- colnames(x = seurat.assay)
        names(x = cells) <- cells
        frags[[1]] <- CreateFragmentObject(
          path = fragments,
          cells = cells,
          validate.fragments = validate.fragments,
          verbose = verbose,
          ...
        )
      }
    }
  }

  chrom.assay <- as.ChromatinAssay(
    x = seurat.assay,
    ranges = ranges,
    seqinfo = genome,
    motifs = motifs,
    fragments = frags,
    annotation = annotation,
    bias = bias,
    positionEnrichment = positionEnrichment
  )
  return(chrom.assay)
}

#' @param ranges A GRanges object
#' @param seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotation Genomic annotation
#' @param motifs A \code{\link{Motif}} object
#' @param fragments A list of \code{\link{Fragment}} objects
#' @param bias Tn5 integration bias matrix
#' @param positionEnrichment A named list of position enrichment matrices.
#' @param sep Characters used to separate the chromosome, start, and end
#' coordinates in the row names of the data matrix
#'
#' @rdname as.ChromatinAssay
#' @export
#' @method as.ChromatinAssay Assay
#' @concept assay
#'
as.ChromatinAssay.Assay <- function(
  x,
  ranges = NULL,
  seqinfo = NULL,
  annotation = NULL,
  motifs = NULL,
  fragments = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c("-", "-"),
  ...
) {
  new.assay <- as(object = x, Class = "ChromatinAssay")
  ranges <- SetIfNull(
    x = ranges,
    y = StringToGRanges(regions = rownames(x = x), sep = sep)
  )
  new.assay <- SetAssayData(
    object = new.assay,
    slot = "ranges",
    new.data = ranges
  )
  if (!is.null(x = fragments)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = "fragments",
      new.data = fragments
    )
  }
  if (!is.null(x = seqinfo)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = "seqinfo",
      new.data = seqinfo
    )
  }
  if (!is.null(x = annotation)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = "annotation",
      new.data = annotation
    )
  }
  if (!is.null(x = motifs)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = "motifs",
      new.data = motifs
    )
  }
  if (!is.null(x = bias)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = "bias",
      new.data = bias
    )
  }
  if (!is.null(x = positionEnrichment)) {
    new.assay <- SetAssayData(
      object = new.assay,
      slot = "positionEnrichment",
      new.data = positionEnrichment
    )
  }
  return(new.assay)
}

setAs(
  from = "Assay",
  to = "ChromatinAssay",
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
        "Class" = "ChromatinAssay"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)

## Functions

#' Create motif object
#'
#' Create a \code{\link{Motif-class}} object.
#'
#' @param data A motif x region matrix
#' @param pwm A named list of position weight matrices or position frequency
#' matrices matching the motif names in \code{data}.
#' Can be of class PFMatrixList.
#' @param motif.names A named list of motif names. List element names
#' must match the names given in \code{pwm}. If NULL, use the names from the
#' list of position weight or position frequency matrices. This can be used to
#' set a alternative common name for the motif. If a PFMatrixList is passed to
#' \code{pwm}, it will pull the motif name from the PFMatrixList.
#' @param positions A \code{\link[GenomicRanges]{GRangesList}} object containing
#' exact positions of each motif.
#' @param meta.data A data.frame containing metadata
#' @export
#' @return Returns a \code{\link{Motif}} object
#' @concept motifs
#' @examples
#' motif.matrix <- matrix(
#'   data = sample(c(0,1),
#'     size = 100,
#'     replace = TRUE),
#'   ncol = 5
#' )
#' motif <- CreateMotifObject(data = motif.matrix)
CreateMotifObject <- function(
  data = NULL,
  pwm = NULL,
  motif.names = NULL,
  positions = NULL,
  meta.data = NULL
) {
  data <- SetIfNull(x = data, y = new(Class = "dgCMatrix"))
  meta.data <- SetIfNull(x = meta.data, y = data.frame())
  if (
    !(inherits(x = data, what = "matrix") |
      inherits(x = data, what = "dgCMatrix"))
    ) {
    stop("Data must be matrix or sparse matrix class. Supplied ",
         class(x = data))
  }
  if (inherits(x = data, what = "matrix")) {
    data <- as(Class = "dgCMatrix", object = data)
  }
  if ((nrow(x = data) > 0) & (length(x = pwm) > 0)) {
    if (!all(names(x = pwm) == colnames(x = data))) {
      stop("Motif names in data matrix and PWM list are inconsistent")
    }
  }
  if ((nrow(x = data) > 0) & (nrow(x = meta.data) > 0)) {
    if (!all(rownames(x = meta.data) == rownames(x = data))) {
      stop("Motif names in data matrix and metadata are inconsistent")
    }
  }
  if (inherits(x = pwm, what = "list")) {
    if (is.null(names(x = pwm))) {
      stop("PWM must be a named list")
    }
  }
  if (!is.null(x = motif.names)) {
    if (length(x = motif.names) != length(x = pwm)) {
      stop("Number of motif names supplied does not match the number of motifs")
    }
  }
  if (
    inherits(x = pwm, what = "PFMatrixList") |
    inherits(x = pwm, what = "PWMatrixList")
    ) {
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
    Class = "Motif",
    data = data,
    pwm = pwm,
    motif.names = motif.names,
    positions = positions,
    meta.data = meta.data
  )
  return(motif.obj)
}

#' @importFrom Seurat GetAssayData
#' @method GetAssayData ChromatinAssay
#' @export
#' @concept assay
GetAssayData.ChromatinAssay <- function(
  object,
  slot = "data",
  assay = NULL,
  ...
) {
  if (!(slot %in% slotNames(x = object))) {
    stop(
      "slot must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  return(slot(object = object, name = slot))
}

#' Get Fragment object data
#'
#' Extract data from a \code{\link{Fragment-class}} object
#'
#' @param object A \code{\link{Fragment}} object
#' @param slot Information to pull from object (path, hash, cells, prefix, suffix)
#' @export
#' @concept assay
GetFragmentData <- function(object, slot = "path") {
  return(slot(object = object, name = slot))
}

#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
#' @concept motifs
#' @examples
#' motif.obj <- Seurat::GetAssayData(
#'   object = atac_small[['peaks']], slot = "motifs"
#' )
#' GetMotifData(object = motif.obj)
GetMotifData.Motif <- function(object, slot = "data", ...) {
  return(slot(object = object, name = slot))
}

#' @importFrom Seurat GetAssayData
#' @rdname GetMotifData
#' @concept motifs
#' @method GetMotifData ChromatinAssay
#' @export
GetMotifData.ChromatinAssay <- function(object, slot = "data", ...) {
  motif.obj <- GetAssayData(object = object, slot = "motifs")
  if (is.null(x = motif.obj)) {
    stop("Motif object not present in assay")
  } else {
    return(GetMotifData(object = motif.obj, slot = slot, ...))
  }
}

#' @param assay Which assay to use. Default is the current active assay
#' @rdname GetMotifData
#' @method GetMotifData Seurat
#' @concept motifs
#' @importFrom Seurat DefaultAssay GetAssay
#' @export
#' @examples
#' GetMotifData(object = atac_small)
GetMotifData.Seurat <- function(object, assay = NULL, slot = "data", ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  return(GetMotifData(
    object = GetAssay(object = object, assay = assay),
    slot = slot,
    ...
  ))
}

#' @importFrom Seurat RenameCells
#' @importFrom Seurat GetAssayData
#' @concept assay
#' @method RenameCells ChromatinAssay
#' @export
RenameCells.ChromatinAssay <- function(object, new.names = NULL, ...) {
  names(x = new.names) <- colnames(x = object)
  for (i in seq_along(along.with = Fragments(object = object))) {
    slot(object = object, name = "fragments")[[i]] <- RenameCells(
      object = slot(object = object, name = "fragments")[[i]],
      new.names = new.names
    )
  }
  pos.enrich <- GetAssayData(object = object, slot = "positionEnrichment")
  for (i in seq_along(along.with = pos.enrich)) {
    mat <- pos.enrich[[i]]
    mat <- mat[colnames(x = object), ]
    rownames(x = mat) <- new.names[rownames(x = mat)]
    pos.enrich[[i]] <- mat
  }
  slot(object = object, name = "positionEnrichment") <- pos.enrich

  # TODO need to convert to standard assay, rename cells, convert back
  # this would account for possibility of SCT-normalized data in a ChrAssay
  names(x = new.names) <- NULL
  for (data.slot in c("counts", "data", "scale.data")) {
    old.data <- GetAssayData(object = object, slot = data.slot)
    if (ncol(x = old.data) <= 1) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }
  return(object)
}

#' @importFrom Seurat RenameCells
#' @concept fragments
#' @method RenameCells Fragment
#' @export
RenameCells.Fragment <- function(object, new.names, ...) {
  cells <- GetFragmentData(object = object, slot = "cells")
  if (is.null(x = cells)) {
    stop("Cannot rename cells in Fragment object ",
         "with no cell information stored")
  }
  cells <- cells[names(x = new.names)]
  names(x = cells) <- new.names[names(x = cells)]
  slot(object = object, name = "cells") <- cells
  return(object)
}

#' @importFrom Seurat SetAssayData
#' @importFrom GenomeInfoDb genome Seqinfo
#' @method SetAssayData ChromatinAssay
#' @concept assay
#' @export
SetAssayData.ChromatinAssay <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop(
      "slot must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  if (slot %in% c("counts", "data", "scale.data")) {
    if (!(is(object = new.data, class2 = "AnyMatrix"))) {
      stop("Data must be a matrix or sparseMatrix")
    }
    if (ncol(x = object) != ncol(x = new.data)) {
      stop("Number of columns in the provided matrix does not match
           the number of cells in the object")
    }
    if (slot %in% c("counts", "data")) {
      if (nrow(x = object) != nrow(x = new.data)) {
        stop("Number of rows in provided matrix does not match
           the number of rows in the object")
      }
    } else {
      # scale data
      if (nrow(x = object) < nrow(x = new.data)) {
        stop("Number of rows in provided matrix is greater than
             the number of rows in the object")
      }
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "seqinfo") {
    if (inherits(x = new.data, what = "Seqinfo")) {
      slot(object = object, name = slot) <- new.data
    } else if (is(object = new.data, class2 = "character")) {
      slot(object = object, name = slot) <- Seqinfo(genome = new.data)
    } else if(is.null(x = new.data)) {
      slot(object = object, name = slot) <- NULL
    } else {
      stop("Unknown object supplied. Choose a Seqinfo object or the name
           of a UCSC genome")
    }
  } else if (slot == "fragments") {
    if (inherits(x = new.data, what = "list")) {
      # check that it's a list containing fragment class objects
      for (i in seq_along(new.data)) {
        if (!inherits(x = new.data[[i]], what = "Fragment")) {
          stop("New data is not a Fragment object")
        }
      }
    } else if (inherits(x = new.data, what = "Fragment")) {
      # single fragment object
      new.data <- list(new.data)
    }
    frag.list <- GetAssayData(object = object, slot = "fragments")
    if (length(x = frag.list) != 0) {
      warning("Overwriting existing fragment objects")
    }
    slot(object = object, name = "fragments") <- new.data
  } else if (slot == "annotation") {
    if (!is(object = new.data, class2 = "GRanges")) {
      stop("Must provide a GRanges object")
    }
    current.genome <- unique(x = genome(x = object))
    annotation.genome <- unique(x = genome(x = new.data))
    if (!is.null(x = current.genome)) {
      if (!is.na(x = annotation.genome) &
          (current.genome != annotation.genome)) {
        stop("Annotation genome does not match genome of the object")
        }
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "bias") {
    if (!is(object = new.data, class2 = "vector")) {
      stop("Bias must be provided as a vector")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "positionEnrichment") {
    if (inherits(x = new.data, what = "list")) {
      # list of position enrichment matrices being added
      if (length(x = new.data) == 0) {
        # if list is empty, assign and overwrite slot
        slot(object = object, name = slot) <- new.data
      } else if (is.null(x = names(x = new.data))) {
        stop("If supplying a list of position enrichment matrices,
             each element must be named")
      } else {
        current.data <- GetAssayData(object = object, slot = slot)
        if (length(x = current.data) != 0) {
          warning("Overwriting current list of position enrichement matrices")
        }
        for (i in seq_along(along.with = new.data)) {
          if (!is(object = new.data[[i]], class2 = "AnyMatrix")) {
            stop(
              "Position enrichment must be provided as a matrix or sparseMatrix"
              )
          }
        }
        slot(object = object, name = slot) <- new.data
      }
    } else if (!is(object = new.data, class2 = "AnyMatrix")) {
      stop("Position enrichment must be provided as a matrix or sparseMatrix")
    } else {
      # single new matrix being added, needs a key
      args <- list(...)
      if (!("key" %in% names(x = args))) {
        stop("Must supply a key when adding positionEnrichment data")
      } else {
        key <- args$key
      }
      current.pos <- slot(object = object, name = slot)
      current.pos[[key]] <- new.data
      slot(object = object, name = slot) <- current.pos
    }
  } else if (slot == "ranges") {
    if (!is(object = new.data, class2 = "GRanges")) {
      stop("Must provide a GRanges object")
    } else if (length(x = new.data) != nrow(x = object)) {
      stop("Number of ranges provided is not equal to the number
           of features in the assay")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "motifs") {
    if (!inherits(x = new.data, what = "Motif")) {
      stop("Must provide a Motif class object")
    }
    # TODO allow mismatching row names, but check that the genomic ranges
    # are equivalent. Requires adding a granges slot to the motif class
    if (!all(rownames(x = object) == rownames(x = new.data))) {
      keep.features <- intersect(x = rownames(x = new.data),
                                 y = rownames(x = object))
      if (length(x = keep.features) == 0) {
        stop("No features in common between the ChromatinAssay
             and Motif objects")
      }
      else {
        warning("Features do not match in ChromatinAssay and Motif object.
                Subsetting the Motif object.")
        new.data <- new.data[keep.features, ]
      }
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "links") {
    slot(object = object, name = slot) <- new.data
  }
  return(object)
}

#' @rdname SetMotifData
#' @method SetMotifData Motif
#' @export
#' @concept motifs
#' @examples
#' motif.obj <- Seurat::GetAssayData(
#'   object = atac_small[['peaks']], slot = "motifs"
#' )
#' SetMotifData(object = motif.obj, slot = 'data', new.data = matrix())
SetMotifData.Motif <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop("slot must be one of ",
         paste(slotNames(x = object), collapse = ", "),
         call. = FALSE)
  }
  if (slot == "data") {
    if (inherits(x = new.data, what = "matrix")) {
      new.data <- as(Class = "dgCMatrix", object = new.data)
    }
  }
  # TODO check that new data is compatible with existing slots
  # rownames of data must match rownames of meta.data and names of pwm, if not
  # empty
  slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param new.data motif matrix to add. Should be matrix or sparse matrix class
#' @param slot Name of slot to use
#' @importFrom Seurat GetAssayData SetAssayData
#' @rdname SetMotifData
#' @export
#' @concept motifs
#' @examples
#' SetMotifData(
#'   object = atac_small[['peaks']], slot = 'data', new.data = matrix()
#' )
#' @method SetMotifData ChromatinAssay
SetMotifData.ChromatinAssay <- function(object, slot, new.data, ...) {
  if (slot == "data") {
    if (
      !(inherits(x = new.data, what = "matrix") |
        inherits(x = new.data, what = "dgCMatrix"))
      ) {
      stop("Data must be matrix or sparse matrix class. Supplied ",
           class(x = new.data))
    }
    if (!all(rownames(x = object) == rownames(x = new.data))) {
      stop("Features do not match existing assay data.
           Column names in motif matrix should match row names in assay data")
    }
    if (inherits(x = new.data, what = "matrix")) {
      new.data <- as(Class = "dgCMatrix", object = new.data)
    }
  }
  motif.obj <- GetAssayData(object = object, slot = "motifs")
  if (is.null(x = motif.obj)) {
    stop("Motif object not present in assay")
  } else {
    motif.obj <- SetMotifData(
      object = motif.obj, slot = slot, new.data = new.data
    )
    object <- SetAssayData(
      object = object, slot = "motifs", new.data = motif.obj
    )
    return(object)
  }
}

#' @param assay Name of assay whose data should be set
#' @rdname SetMotifData
#' @importFrom Seurat DefaultAssay
#' @export
#' @method SetMotifData Seurat
#' @concept motifs
#' @examples
#' motif.matrix <- GetMotifData(object = atac_small)
#' SetMotifData(
#' object = atac_small, assay = 'peaks', slot = 'data', new.data = motif.matrix
#' )
SetMotifData.Seurat <- function(object, assay = NULL, ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  object[[assay]] <- SetMotifData(object = object[[assay]], ...)
  return(object)
}

#' Subset a Motif object
#'
#' Returns a subset of a \code{\link{Motif-class}} object.
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
#' @concept motifs
#' @examples
#' motif.obj <- Seurat::GetAssayData(
#'   object = atac_small[['peaks']], slot = "motifs"
#' )
#' subset(x = motif.obj, features = head(rownames(motif.obj), 10))
subset.Motif <- function(x, features = NULL, motifs = NULL, ...) {
  features <- SetIfNull(x = features, y = rownames(x = x))
  motifs <- SetIfNull(x = motifs, y = colnames(x = x))
  new.data <- GetMotifData(object = x, slot = "data")[features, motifs]
  new.pwm <- GetMotifData(object = x, slot = "pwm")[motifs]
  new.names <- GetMotifData(object = x, slot = "motif.names")[motifs]
  new.meta <- GetMotifData(object = x, slot = "meta.data")[motifs, ]
  new.positions <- GetMotifData(object = x, slot = "positions")
  if (!is.null(x = new.positions)) {
    new.positions <- new.positions[motifs]
  }
  new.motif <- new(
    Class = "Motif",
    data = new.data,
    pwm = new.pwm,
    motif.names = new.names,
    meta.data = new.meta,
    positions = new.positions
  )
  return(new.motif)
}

#' @export
#' @importClassesFrom Seurat Assay
#' @concept assay
#' @method subset ChromatinAssay
subset.ChromatinAssay <- function(
  x,
  features = NULL,
  cells = NULL,
  ...
) {
  # subset elements in the standard assay
  standardassay <- as(object = x, Class = "Assay")
  standardassay <- subset(x = standardassay, features = features, cells = cells)

  # recompute meta features
  standardassay <- FindTopFeatures(
    object = standardassay,
    min.cutoff = NA,
    verbose = FALSE
  )

  # subset genomic ranges
  ranges.keep <- granges(x = x)
  if (!is.null(x = features)) {
    idx.keep <- rownames(x = x) %in% features
    ranges.keep <- ranges.keep[idx.keep]
  }

  # subset motifs
  motifs <- Motifs(object = x)
  if (!is.null(x = motifs)) {
    motifs <- subset(x = motifs, features = features)
  }

  # subset cells in positionEnrichment matrices
  cells <- SetIfNull(x = cells, y = colnames(x = x))
  posmat <- GetAssayData(object = x, slot = "positionEnrichment")
  for (i in seq_along(along.with = posmat)) {
    posmat[[i]] <- posmat[[i]][cells, ]
  }

  # subset cells in Fragments objects
  frags <- Fragments(object = x)
  for (i in seq_along(along.with = frags)) {
    frag.cells <- GetFragmentData(object = frags[[i]], slot = "cells")
    # there can be cells in the assay that are not in the fragment object
    keep <- names(x = frag.cells) %in% cells
    slot(object = frags[[i]], name = "cells") <- frag.cells[keep]
  }

  # convert standard assay to ChromatinAssay
  chromassay <- as.ChromatinAssay(
    x = standardassay,
    ranges = ranges.keep,
    seqinfo = seqinfo(x = x),
    annotation = Annotation(object = x),
    motifs = motifs,
    fragments = frags,
    bias = GetAssayData(object = x, slot = "bias"),
    positionEnrichment = posmat
  )
  return(chromassay)
}

#' @export
#' @concept assay
#' @method merge ChromatinAssay
#' @importFrom GenomicRanges union findOverlaps
#' @importFrom Seurat RowMergeSparseMatrices
#' @importFrom S4Vectors subjectHits queryHits mcols
#' @importMethodsFrom GenomeInfoDb merge
merge.ChromatinAssay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  ...
) {
  # need to do all operations over a list of assays
  assays <- c(x, y)

  # if any are standard Assay class, coerce all to Assay and run merge
  isChromatin <- sapply(
    X = assays, FUN = function(x) inherits(x = x, what = "ChromatinAssay")
  )
  if (!all(isChromatin)) {
    # check that the non-chromatinassays have >1 feature
    nfeature <- sapply(X = assays, FUN = nrow)
    if (all(nfeature > 1)) {
      # genuine assays, coerce to standard assay and run merge.Assay
      warning(
        "Some assays are not ChromatinAssay class, ",
        "coercing ChromatinAssays to standard Assay"
      )
      assays <- sapply(
        X = assays, FUN = function(x) as(object = x, Class = "Assay")
      )
      new.assay <- merge(
        x = assays[[1]], y = assays[[2:length(x = assays)]], ...
      )
      return(new.assay)
    } else {
      # Find which assays are placeholder
      placeholders <- nfeature == 1 & !isChromatin
      # Set feature name as first peak in first real assay
      peak.use <- rownames(x = assays[isChromatin][[1]])[1]
      converted <- sapply(
        X = assays[placeholders], FUN = function(x) {
          rownames(x = x@counts) <- peak.use
          rownames(x = x@data) <- peak.use
          return(x)
        }
      )
      # Covert placeholder assays to ChromatinAssay
      converted <- sapply(
        X = converted, FUN = function(x) as.ChromatinAssay(x = x)
      )
      # Replace original assays
      assays[placeholders] <- converted
      # Continue with merge function
    }
  }

  # rename cells in each assay
  # merge.Seurat already does this, so should only happen here when merging
  # assay objects outside of a Seurat object
  if (is.null(x = add.cell.ids)) {
    # check if any cell names clash, if so add a prefix
    cellnames.all <- sapply(X = assays, FUN = colnames)
    cellnames.all <- Reduce(f = c, x = cellnames.all)
    cellname.freq <- table(cellnames.all)
    if (max(cellname.freq) > 1) {
      message(
        "Cell names not unique, ",
        "adding prefix to enforce unique cell names"
      )
      add.cell.ids <- seq_along(along.with = assays)
    }
  }
  if (!is.null(x = add.cell.ids)) {
    for (i in seq_along(along.with = assays)) {
      assays[[i]] <- RenameCells(
        object = assays[[i]],
        new.names = paste(add.cell.ids[i], colnames(x = assays[[i]]), sep = "_")
      )
    }
  }

  # check genomes are all the same
  genomes <- unlist(
    x = lapply(X = assays, FUN = function(x) unique(x = genome(x = x)))
  )
  if (length(x = unique(x = genomes)) > 1) {
    warning("Genomes do not match, not merging ChromatinAssays")
    return(NULL)
  }

  # merge seqinfo
  all.seqinfo <- lapply(X = assays, FUN = function(x) seqinfo(x = x))
  seqinfo.present <- !sapply(X = all.seqinfo, FUN = is.null)
  if (any(seqinfo.present)) {
    # need at least one non-NULL seqinfo, otherwise just set it as NULL
    all.seqinfo <- all.seqinfo[seqinfo.present]
    if (length(x = all.seqinfo) > 1) {
      seqinfo.use <- all.seqinfo[[1]]
      # iteratively merge seqinfo objects
      for (x in 2:length(x = all.seqinfo)) {
        seqinfo.use <- merge(x = seqinfo.use, y = all.seqinfo[[x]])
      }
    } else {
      seqinfo.use <- all.seqinfo[[1]]
    }
  } else {
    seqinfo.use <- NULL
  }

  # merge annotations
  all.annot <- lapply(X = assays, FUN = function(x) Annotation(object = x))
  annot.present <- !sapply(X = all.annot, FUN = is.null)
  if (any(annot.present)) {
    all.annot <- all.annot[annot.present]
    annot.use <- all.annot[[1]]
    if (length(x = all.annot) > 1) {
      for (x in 2:length(x = all.annot)) {
        if (!identical(x = annot.use, y = all.annot[[x]])) {
          warning("Annotations do not match, keeping annotation from the
            first object only")
        }
      }
    }
  } else {
    annot.use <- NULL
  }

  # merge fragments
  all.frag <- lapply(X = assays, FUN = function(x) Fragments(object = x))
  all.frag <- Reduce(f = c, x = all.frag)
  valid.frags <- sapply(X = all.frag, FUN = ValidateHash, verbose = FALSE)
  if (!all(valid.frags)) {
    warning("Some fragment files are not valid or not indexed.
            Removing invalid files from merged ChromatinAssay")
    all.frag <- all.frag[valid.frags]
  }

  # check that all features are equal
  all.features <- lapply(X = assays, FUN = rownames)
  all.features <- table(do.call(what = c, args = all.features))
  all.identical <- all(all.features == length(x = assays))

  # find whether the unique ranges are all disjoint
  all.nonoverlapping <- NonOverlapping(x = assays, all.features = all.features)

  if (all.identical | all.nonoverlapping) {
    # no non-identical but overlapping features present
    merged.counts <- list()
    merged.data <- list()
    if (all.identical) {
      feat.use <- rownames(x = assays[[1]])
      for (i in seq_along(along.with = assays)) {
        # check that counts are present
        # can be removed by DietSeurat
        assay.counts <- GetAssayData(object = assays[[i]], slot = "counts")
        if (nrow(x = assay.counts) > 0) {
          merged.counts[[i]] <- assay.counts[feat.use, ]
        } else {
          merged.counts[[i]] <- assay.counts
        }
        merged.data[[i]] <- GetAssayData(
          object = assays[[i]], slot = "data"
        )[feat.use, ]
      }
      # exact same features, can just run cbind
      # can also merge data and scaledata
      merged.counts <- do.call(what = cbind, args = merged.counts)
      merged.data <- do.call(what = cbind, args = merged.data)
      reduced.ranges <- granges(x = assays[[1]])
    } else {
      # disjoint
      all.counts <- list()
      all.data <- list()
      for (i in seq_along(along.with = assays)) {
        all.counts[[i]] <- GetAssayData(object = assays[[i]], slot = "counts")
        all.data[[i]] <- GetAssayData(object = assays[[i]], slot = "data")
      }
      count_nonzero <- lapply(X = all.counts, FUN = ncol)
      data_nonzero <- lapply(X = all.data, FUN = ncol)
      if (all(count_nonzero > 0)) {
        merged.counts <- RowMergeSparseMatrices(
          mat1 = all.counts[[1]],
          mat2 = all.counts[[2:length(x = all.counts)]]
        )
        reduced.ranges <- StringToGRanges(regions = rownames(x = merged.counts))
      } else {
        merged.counts <- matrix(nrow = 0, ncol = 0)
        reduced.ranges <- NULL
      }
      if (all(data_nonzero > 0)) {
        merged.data <- RowMergeSparseMatrices(
          mat1 = all.data[[1]],
          mat2 = all.data[[2:length(x = all.data)]]
        )
        reduced.ranges <- SetIfNull(
          x = reduced.ranges,
          y = StringToGRanges(regions = rownames(x = merged.data))
        )
      } else {
        merged.data <- matrix(nrow = 0, ncol = 0)
      }
      if (is.null(x = reduced.ranges)) {
        stop("No counts or data in the assay")
      }
    }

    # create new ChromatinAssay object
    # bias, motifs, positionEnrichment, metafeatures not kept
    # scaledata only kept if features exactly identical
    if (nrow(x = merged.counts) > 0) {
      new.assay <- CreateChromatinAssay(
        counts = merged.counts,
        min.cells = -1,
        min.features = -1,
        max.cells = NULL,
        ranges = reduced.ranges,
        motifs = NULL,
        fragments = all.frag,
        genome = seqinfo.use,
        annotation = annot.use,
        bias = NULL,
        validate.fragments = FALSE
      )
      new.assay <- SetAssayData(
        object = new.assay, slot = "data", new.data = merged.data
      )
    } else {
      new.assay <- CreateChromatinAssay(
        data = merged.data,
        min.cells = -1,
        min.features = -1,
        max.cells = NULL,
        ranges = reduced.ranges,
        motifs = NULL,
        fragments = all.frag,
        genome = seqinfo.use,
        annotation = annot.use,
        bias = NULL,
        validate.fragments = FALSE
      )
    }
  } else {
    # first create a merged set of granges, preserving the assay of origin
    granges.all <- sapply(X = assays, FUN = granges)
    for (i in seq_along(along.with = granges.all)) {
      granges.all[[i]]$dataset <- i
    }
    granges.all <- Reduce(f = c, x = granges.all)

    # create reduced ranges, recording the indices of the merged ranges
    reduced.ranges <- reduce(x = granges.all, with.revmap = TRUE)

    # get the new rownames for the count matrix
    new.rownames <- GRangesToString(grange = reduced.ranges)

    # function to look up original
    tomerge <- GetRowsToMerge(
      assay.list = assays,
      all.ranges = granges.all,
      reduced.ranges = reduced.ranges
    )

    # if the grange is the same, merge matrix rows
    merged.counts <- MergeOverlappingRows(
      mergeinfo = tomerge,
      assay.list = assays,
      slot = "counts",
      verbose = TRUE
    )

    merged.data <- MergeOverlappingRows(
      mergeinfo = tomerge,
      assay.list = assays,
      slot = "data",
      verbose = TRUE
    )

    if (nrow(x = merged.counts[[1]]) > 0) {
      merged.counts <- MergeMatrixParts(
        mat.list = merged.counts,
        new.rownames = new.rownames
      )
      merged.data <- MergeMatrixParts(
        mat.list = merged.data,
        new.rownames = new.rownames
      )
      new.assay <- CreateChromatinAssay(
        counts = merged.counts,
        min.cells = -1,
        min.features = -1,
        max.cells = NULL,
        ranges = reduced.ranges,
        motifs = NULL,
        fragments = all.frag,
        genome = seqinfo.use,
        annotation = annot.use,
        bias = NULL,
        validate.fragments = FALSE
      )
      new.assay <- SetAssayData(
        object = new.assay, slot = "data", new.data = merged.data
      )
    } else {
      merged.data <- MergeMatrixParts(
        mat.list = merged.data,
        new.rownames = new.rownames
      )
      # create new ChromatinAssay object
      # bias, motifs, positionEnrichment, metafeatures not kept
      # need to keep data otherwise integration doesn't work
      new.assay <- CreateChromatinAssay(
        data = merged.data,
        min.cells = 0,
        min.features = 0,
        max.cells = NULL,
        ranges = reduced.ranges,
        motifs = NULL,
        fragments = all.frag,
        genome = seqinfo.use,
        annotation = annot.use,
        bias = NULL,
        validate.fragments = FALSE
      )
    }
  }
  return(new.assay)
}

#' @param i Which columns to retain
#' @param j Which rows to retain
#'
#' @rdname subset.Motif
#' @export
#' @concept motifs
#' @method [ Motif
#' @examples
#' motif.obj <- Seurat::GetAssayData(
#'   object = atac_small, assay = 'peaks', slot = 'motifs'
#' )
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
  f = "show",
  signature = "Motif",
  definition = function(object) {
    cat(
      "A Motif object containing",
      ncol(x = slot(object = object, name = "data")),
      "motifs in",
      nrow(x = slot(object = object, name = "data")),
      "regions\n"
    )
  }
)

setMethod(
  f = "show",
  signature = "Fragment",
  definition = function(object) {
    cat(
      "A Fragment object for",
      length(x = slot(object = object, name = "cells")),
      "cells\n"
    )
  }
)

setMethod(
  f = "show",
  signature = "ChromatinAssay",
  definition = function(object) {
    cat(
      "ChromatinAssay data with",
      nrow(x = object),
      "features for",
      ncol(x = object),
      "cells\n"
    )
    cat(
      "Variable features:",
      length(x = VariableFeatures(object = object)),
      "\n"
    )
    cat(
      "Genome:",
      unique(x = genome(x = object)),
      "\n"
    )
    cat(
      "Annotation present:",
      ifelse(
        test = is.null(x = Annotation(object = object)), yes = FALSE, no = TRUE
      ),
      "\n"
    )
    cat(
      "Motifs present:",
      ifelse(
        test = is.null(x = Motifs(object = object)),
        yes = FALSE,
        no = TRUE
      ),
      "\n"
    )
    cat(
      "Fragment files:",
      length(x = Fragments(object = object)),
      "\n"
    )
  }
)

#' @rdname Annotation
#' @method Annotation ChromatinAssay
#' @export
#' @concept assay
#' @examples
#' \donttest{
#' Annotation(atac_small[["peaks"]])
#' }
Annotation.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = "annotation"))
}

#' @param object A Seurat object or ChromatinAssay object
#' @importFrom Seurat DefaultAssay
#' @rdname Annotation
#' @method Annotation Seurat
#' @export
#' @concept assay
#' @examples
#' \donttest{
#' Annotation(atac_small)
#' }
Annotation.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Annotation(object = object[[assay]]))
}

#' @rdname Fragments
#' @method Fragments ChromatinAssay
#' @export
#' @concept assay
#' @concept fragments
#' @examples
#' Fragments(atac_small[["peaks"]])
Fragments.ChromatinAssay <- function(object, ...) {
  return(slot(object, name = "fragments"))
}

#' @param object A Seurat object or ChromatinAssay object
#' @importFrom Seurat DefaultAssay
#' @rdname Fragments
#' @method Fragments Seurat
#' @export
#' @concept assay
#' @concept fragments
#' @examples
#' Fragments(atac_small)
Fragments.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Fragments(object = object[[assay]]))
}

#' @rdname Motifs
#' @method Motifs ChromatinAssay
#' @export
#' @concept assay
#' @concept motifs
#' @examples
#' Motifs(atac_small[["peaks"]])
Motifs.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = "motifs"))
}

#' @param object A Seurat object
#' @rdname Motifs
#' @importFrom Seurat DefaultAssay
#' @method Motifs Seurat
#' @export
#' @concept assay
#' @concept motifs
#' @examples
#' Motifs(atac_small)
Motifs.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Motifs(object = object[[assay]]))
}

#' @rdname Links
#' @method Links ChromatinAssay
#' @export
#' @concept assay
#' @concept links
#' @examples
#' Links(atac_small[["peaks"]])
Links.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = "links"))
}

#' @param object A Seurat object
#' @rdname Links
#' @method Links Seurat
#' @importFrom Seurat DefaultAssay
#' @export
#' @concept links
#' @concept assay
#' @examples
#' Links(atac_small)
Links.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Links(object = object[[assay]]))
}

#' @method dimnames Motif
#' @concept motifs
#' @export
dimnames.Motif <- function(x) {
  return(dimnames(x = GetMotifData(object = x)))
}

#' @method dim Motif
#' @concept motifs
#' @export
dim.Motif <- function(x) {
  return(dim(x = GetMotifData(object = x)))
}

#' @export
#' @rdname Motifs
#' @method Motifs<- ChromatinAssay
#' @concept assay
#' @concept motifs
#' @examples
#' motifs <- Motifs(atac_small)
#' Motifs(atac_small[["peaks"]]) <- motifs
"Motifs<-.ChromatinAssay" <- function(object, ..., value) {
  object <- SetAssayData(object = object, slot = "motifs", new.data = value)
  return(object)
}

#' @export
#' @rdname Motifs
#' @method Motifs<- Seurat
#' @importFrom Seurat DefaultAssay
#' @concept assay
#' @concept motifs
#' @examples
#' motifs <- Motifs(atac_small)
#' Motifs(atac_small) <- motifs
"Motifs<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  Motifs(object = object[[assay]]) <- value
  return(object)
}

#' @export
#' @rdname Links
#' @method Links<- ChromatinAssay
#' @concept assay
#' @concept links
#' @examples
#' links <- Links(atac_small)
#' Links(atac_small[["peaks"]]) <- links
"Links<-.ChromatinAssay" <- function(object, ..., value) {
  object <- SetAssayData(object = object, slot = "links", new.data = value)
  return(object)
}

#' @export
#' @rdname Links
#' @method Links<- Seurat
#' @concept assay
#' @concept links
#' @examples
#' links <- Links(atac_small)
#' Links(atac_small) <- links
"Links<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  Links(object[[assay]]) <- value
  return(object)
}

#' @export
#' @rdname Annotation
#' @concept assay
#' @method Annotation<- ChromatinAssay
#' @examples
#' genes <- Annotation(atac_small)
#' Annotation(atac_small[["peaks"]]) <- genes
"Annotation<-.ChromatinAssay" <- function(object, ..., value) {
  object <- SetAssayData(object = object, slot = "annotation", new.data = value)
  return(object)
}

#' @export
#' @importFrom Seurat DefaultAssay
#' @method Annotation<- Seurat
#' @concept assay
#' @rdname Annotation
#' @examples
#' genes <- Annotation(atac_small)
#' Annotation(atac_small) <- genes
"Annotation<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  Annotation(object = object[[assay]]) <- value
  return(object)
}

#' @export
#' @method Fragments<- ChromatinAssay
#' @rdname Fragments
#' @importFrom Seurat SetAssayData
#' @concept assay
#' @concept fragments
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small[["bins"]]) <- fragments
"Fragments<-.ChromatinAssay" <- function(object, ..., value) {
  if (is.null(x = value)) {
    slot(object = object, name = "fragments") <- list()
    return(object)
  }
  if (inherits(x = value, what = "list")) {
    for (i in seq_along(along.with = value)) {
      object <- AddFragments(object = object, fragments = value[[i]])
    }
  } else if (is.null(x = value)) {
    object <- SetAssayData(
      object = object,
      slot = "fragments",
      new.data = list()
    )
  } else {
    object <- AddFragments(object = object, fragments = value)
  }
  return(object)
}

#' @export
#' @method Fragments<- Seurat
#' @rdname Fragments
#' @concept assay
#' @concept fragments
#' @importFrom Seurat DefaultAssay
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
"Fragments<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  Fragments(object = object[[assay]]) <- value
  return(object)
}

# Add a single Fragment object to a ChromatinAssay
# @param object A \code{\link{ChromatinAssay}} object
# @param fragments A \code{\link{Fragment}} object
AddFragments <- function(object, fragments) {
  # validate hash
  if (!ValidateHash(object = fragments, verbose = FALSE)) {
    stop("Invalid Fragment object")
  }
  # if cells is NULL, set to all cells in the assay
  # ValidateCells is run in the Cells<- method
  # only allowed if there is no fragment object currently set
  if (is.null(x = Cells(x = fragments))) {
    if (length(x = Fragments(object = object)) != 0) {
      stop("Fragment objects already present in the assay.
           To assign more fragment objects, you must provide a list
           of cells that are contained in each fragment object.")
    } else {
      # each element is the cell name as it appears in the fragment file
      # each element name is the cell name as it appears in the assay
      # here they are assumed to be the same
      cells <- colnames(x = object)
      names(x = cells) <- cells
      Cells(x = fragments) <- cells
    }
  } else {
    # subset cells in the fragment file to those in the assay
    # Cells method returns the names as they appear in the assay
    keep.cells <- Cells(x = fragments) %in% colnames(x = object)
    if (!all(keep.cells)) {
      if (sum(keep.cells) == 0) {
        stop(
          "None of the cells in the fragment object are present in the assay"
        )
      } else {
        # subset the fragment cells, don't need to validate cells again
        # need to make sure to retain the original barcode
        # not the version of the cel name that's stored in the assay
        cell.barcodes <- GetFragmentData(object = fragments, slot = "cells")
        slot(object = fragments, name = "cells") <- cell.barcodes[keep.cells]
      }
    }
    # check that cells not found in any existing fragment objects
    current.frags <- GetAssayData(object = object, slot = "fragments")
    for (i in seq_along(along.with = current.frags)) {
      if (any(Cells(x = fragments) %in% Cells(x = current.frags[[i]]))) {
        stop("Cells already present in a fragment object")
      }
    }
  }
  # append fragments to list
  current.frags <- GetAssayData(object = object, slot = "fragments")
  current.frags[[length(x = current.frags) + 1]] <- fragments
  slot(object = object, name = "fragments") <- current.frags
  return(object)
}

#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod is slot slot<- new as
#' slotNames
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

## Class definitions

setClassUnion(name = "AnyMatrix", c("matrix", "dgCMatrix"))

#' The Fragment class
#'
#' @slot path Path to the fragment file on disk.
#' See \url{https://support.10xgenomics.com/single-cell-atac/software/
#' pipelines/latest/output/fragments}
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file, and the name of each element is the
#' corresponding cell barcode as stored in the ChromatinAssay object.
#'
#' @name Fragment-class
#' @rdname Fragment-class
#' @exportClass Fragment
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
Motif <- setClass(
  Class = "Motif",
  slots = list(
    data = "dgCMatrix",
    pwm = "list",
    motif.names = "list",
    meta.data = "data.frame"
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
#' @slot fragments A list of \code{\link{Fragment}} objects.
#' @slot seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome sequence used.
#' @slot annotation A  \code{\link[GenomicRanges]{GRanges}} object containing
#' genomic annotations
#' @slot bias A matrix containing Tn5 integration bias information
#' (frequency of Tn5 integration at different kmers)
#' @slot positionEnrichment A named list of matrices containing positional
#' enrichment scores for Tn5 integration (for example, enrichment at the TSS)
#'
#' @name ChromatinAssay-class
#' @rdname ChromatinAssay-class
#' @exportClass ChromatinAssay
#'
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
    "positionEnrichment" = "list"
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
#' \code{\link{Fragments}} functions.
#' @param genome A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotation A set of \code{\link[GenomicRanges]{GRanges}} containing
#' annotations for the genome used
#' @param bias A Tn5 integration bias matrix
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
#' @importFrom Matrix rowSums
#'
#' @export
CreateChromatinAssayObject <- function(
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
  if (!is.null(x = annotation) & !inherits(x = annotation, what = "GRanges")) {
    stop("Annotation must be a GRanges object.")
  }
  ncell.feature <- rowSums(data.use > 0)
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
  features.keep <- (ncell.feature > min.cells) & (ncell.feature < max.cells)
  if (!missing(x = counts)) {
    counts <- counts[features.keep, ]
  } else {
    data <- data[features.keep, ]
  }
  ranges <- ranges[features.keep, ]
  seurat.assay <- CreateAssayObject(
    counts = counts,
    data = data,
    min.cells = 0,
    min.features = 0
  )
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
  chrom.assay <- as.ChromatinAssay(
    x = seurat.assay,
    ranges = ranges,
    seqinfo = genome,
    motifs = motifs,
    fragments = frags,
    annotation = annotation,
    bias = bias,
    positionEnrichment = list()
  )
  return(chrom.assay)
}

#' @rdname as.Assay
#' @method as.Assay ChromatinAssay
as.Assay.ChromatinAssay <- function(x) {
  # TODO
  # remove the ChromatinAssay-specific slots and recreate as a standard Assay
  return(x)
}

#' @param ranges A GRanges object
#' @param genome A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotation Genomic annotation
#' @param motifs A \code{\link{Motif}} object
#' @param fragments A list of \code{\link{Fragment}} objects
#' @param bias Tn5 integration bias matrix
#' @param sep Charaters used to separate the chromosome, start, and end
#' coordinates in the row names of the data matrix
#'
#' @rdname as.ChromatinAssay
#' @export
#' @method as.ChromatinAssay Assay
#'
as.ChromatinAssay.Assay <- function(
  x,
  ranges = NULL,
  seqinfo = NULL,
  annotation = NULL,
  motifs = NULL,
  fragments = NULL,
  bias = NULL,
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

#' Create Seurat object with a chromatin assay
#'
#' This is a wrapper for \code{\link[Seurat]{CreateSeuratObject}} to create
#' an object containing a \code{\link{ChromatinAssay}} rather than a
#' standard Seurat \code{\link[Seurat]{Assay}}.
#'
#' @param counts Unnormalized data (raw counts). Should be a sparse matrix.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param genome Name of genome used
#' @param ranges A \code{\link[GenomicRanges]{GRanges}} object containing the
#' genomic position of each row of the counts matrix
#' @param project Sets the project name for the object
#' @param min.cells Include features detected in at least this many cells.
#' Will subset the counts matrix as well. To reintroduce excluded features,
#' create a new object with a lower cutoff.
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
#' @param names.delim For the initial identity class for each cell, choose this
#' delimiter from the cell's column name. E.g. If your cells are named as
#' BARCODE-CLUSTER-CELLTYPE, set this to "-" to separate the cell name into its
#' component parts for picking the relevant field.
#' @param names.field For the initial identity class for each cell, choose this
#' field from the cell's name. E.g. If your cells are named as
#' BARCODE_CLUSTER_CELLTYPE in the input matrix, set names.field to 3 to set the
#' initial identities to CELLTYPE.
#' @param meta.data Additional cell-level metadata to add to the Seurat object.
#' Should be a data frame where the rows are cell names and the columns are
#' additional metadata fields.
#' @param fragments A path to a fragment file on disk. A \code{\link{Fragment}}
#' object will be created from the file path and stored in the object. Only
#' one path should be supplied, and it is assumed that all cells in the input
#' matrix should be present in the fragment file. If you need to create multiple
#' fragment files, you can add them after creating the Seurat object using
#' the \code{\link{CreateFragmentObject}} and \code{\link{Fragments}} functions.
#' @param annotation A \code{\link[GenomicRanges]{GRanges}} object containing
#' genomic annotations for the genome used
#' @param motifs A \code{\link{Motif}} object
#' @param sep Charaters used to separate the chromosome, start, and end
#' coordinates in the row names of the data matrix
#' @param validate.fragments Check that cells in the assay are present in the
#' fragment file.
#' @param verbose Display messages
#' @param ... Additional arguments passed to
#' \code{\link{CreateChromatinAssayObject}}
#'
#' @importFrom Seurat Key<-
#' @export
CreateSignacObject <- function(
  counts,
  assay = "ATAC",
  project = "SignacProject",
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  meta.data = NULL,
  names.delim = "_",
  names.field = 1,
  ranges = NULL,
  fragments = NULL,
  annotation = NULL,
  genome = NULL,
  motifs = NULL,
  sep = c("-", "-"),
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
) {
  ranges <- SetIfNull(
    x = ranges,
    y = StringToGRanges(regions = rownames(x = counts), sep = sep)
  )
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata.
           Please ensure that rownames of metadata match
           column names of data matrix")
    }
    if (
      length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = counts))
      )
    ) {
      warning("Some cells in meta.data not present in provided counts matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data),
                                       y = colnames(x = counts)), ]
    }
    if (class(x = meta.data) == "data.frame") {
      new.meta.data <- data.frame(row.names = colnames(x = counts))
      for (ii in seq_len(length.out = ncol(x = meta.data))) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <-
          meta.data[, ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
  assay.data <- CreateChromatinAssayObject(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    max.cells = max.cells,
    ranges = ranges,
    fragments = fragments,
    annotation = annotation,
    genome = genome,
    motifs = motifs,
    validate.fragments = validate.fragments,
    verbose = verbose,
    ...
  )
  Key(object = assay.data) <- paste0(tolower(x = assay), "_")
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay]]))
  idents <- factor(x = unlist(
    x = lapply(
      X = colnames(x = assay.data),
      FUN = ExtractField,
      field = names.field,
      delim = names.delim)
    )
  )
  if (any(is.na(x = idents))) {
    warning("Input parameters result in NA values for initial cell identities.
            Setting all initial idents to the project name")
  }
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 ||
      ident.levels == length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = assay.data))
  }
  names(x = idents) <- colnames(x = assay.data)
  object <- new(
    Class = "Seurat",
    assays = assay.list,
    meta.data = init.meta.data,
    active.assay = assay,
    active.ident = idents,
    project.name = project,
    version = packageVersion(pkg = "Seurat")
  )
  object[["orig.ident"]] <- idents
  n.calc <- CalcN(object = assay.data)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = "_")
    object[[names(x = n.calc)]] <- n.calc
  }
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

## Functions

#' CreateMotifObject
#'
#' Create an object of class \code{Motif}
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
#' @param meta.data A data.frame containing metadata
#' @export
#' @return Returns a \code{\link{Motif}} object
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
    meta.data = meta.data
  )
  return(motif.obj)
}

#' @importFrom Seurat GetAssayData
#' @method GetAssayData ChromatinAssay
#' @export
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
#' @param object A \code{\link{Fragment}} object
#' @param slot Information to pull from object (path, hash, cells, prefix, suffix)
#' @export
GetFragmentData <- function(object, slot = "path") {
  return(slot(object = object, name = slot))
}

#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
#' @examples
#' motif.obj <- GetMotifObject(object = atac_small[['peaks']])
#' GetMotifData(object = motif.obj)
GetMotifData.Motif <- function(object, slot = "data", ...) {
  return(slot(object = object, name = slot))
}

#' @importFrom Seurat GetAssayData
#' @rdname GetMotifData
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

#' @rdname RenameCells
#' @importFrom Seurat RenameCells
#' @export
#' @method RenameCells ChromatinAssay
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
  object <- Seurat:::RenameCells.Assay(
    object = object, new.names = new.names, ...
  )
  return(object)
}

#' @rdname RenameCells
#' @export
#' @method RenameCells Fragment
RenameCells.Fragment <- function(object, new.names, ...) {
  cells <- GetFragmentData(object = object, slot = "cells")
  # TODO subset cells in Fragment object when subsetting object
  cells <- cells[names(x = new.names)]
  names(x = cells) <- new.names[names(x = cells)]
  slot(object = object, name = "cells") <- cells
  return(object)
}

#' @importFrom Seurat SetAssayData
#' @importFrom GenomeInfoDb genome Seqinfo
#' @method SetAssayData ChromatinAssay
#' @export
SetAssayData.ChromatinAssay <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop(
      "slot must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  if (slot %in% c("counts", "data")) {
    if (!(is(object = new.data, class2 = "AnyMatrix"))) {
      stop("Data must be a matrix or sparseMatrix")
    }
    if (nrow(x = object) != nrow(x = new.data)) {
      stop("Number of rows in provided matrix does not match
           the number of rows in the object")
    }
    if (ncol(x = object) != ncol(x = new.data)) {
      stop("Number of columns in the provided matrix does not match
           the number of cells in the object")
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
    # check that it's a list containing fragment class objects
    for (i in seq_along(along.with = new.data)) {
      if (!inherits(x = new.data[[i]], what = "Fragment")) {
        stop("New data is not a Fragment object")
      }
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
    current.genome <- genome(x = object)
    annotation.genome <- unique(x = genome(x = new.data))
    if (!is.null(x = current.genome) &
        !is.na(x = annotation.genome) &
        (current.genome != annotation.genome)) {
      stop("Annotation genome does not match genome of the object")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "bias") {
    if (!is(object = new.data, class2 = "AnyMatrix")) {
      stop("Bias must be provided as a matrix or sparseMatrix")
    }
    slot(object = object, name = slot) <- new.data
  } else if (slot == "positionEnrichment") {
    if (!is(object = new.data, class2 = "AnyMatrix")) {
      stop("Position enrichment must be provided as a matrix or sparseMatrix")
    }
    args <- list(...)
    if (!("key" %in% names(x = args))) {
      stop("Must supply a key when adding positionEnrichment data")
    } else {
      key <- args$key
    }
    current.pos <- slot(object = object, name = slot)
    current.pos[[key]] <- new.data
    slot(object = object, name = slot) <- current.pos
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
  }
  return(object)
}

#' @rdname SetMotifData
#' @method SetMotifData Motif
#' @export
#' @examples
#' motif.obj <- GetMotifObject(object = atac_small)
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
#' @examples
#' SetMotifData(
#' object = atac_small[['peaks']], slot = 'data', new.data = matrix()
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
  new.data <- GetMotifData(object = x, slot = "data")[features, motifs]
  new.pwm <- GetMotifData(object = x, slot = "pwm")[motifs]
  new.names <- GetMotifData(object = x, slot = "motif.names")[motifs]
  new.meta <- GetMotifData(object = x, slot = "meta.data")[motifs, ]
  new.motif <- new(
    Class = "Motif",
    data = new.data,
    pwm = new.pwm,
    motif.names = new.names,
    meta.data = new.meta
  )
  return(new.motif)
}

#' @export
#' @method subset ChromatinAssay
subset.ChromatinAssay <- function(
  x,
  features = NULL,
  cells = NULL,
  ...
) {
  # TODO
  # need to coerce to standard Assay class, then subset the other parts
  # and re-build the ChromatinAssay
  standardassay <- as.Assay(x = x) # TODO internal conversion function

  # subet genomic ranges
  ranges.keep <- granges(x = x)
  if (!is.null(x = features)) {
    idx.keep <- which(colnames(x = object) == features)
    ranges.keep <- ranges.keep[idx.keep]
  }
  # need to subsect matrix first, otherwise will give errors
  # when dimension doesn't match the matrix dimension
  object <- SetAssayData(
    object = object,
    slot = "ranges",
    new.data = ranges.keep
  )
  # subset motifs
  motifs <- Motif(object = x)
  object <- SetAssayData(
    object = object,
    slot = "motifs",
    new.data = subset(x = motifs, features = features)
  )
  # TODO subset cells in positionEnrichment matrices
  # TODO subset cells in Fragments objects
  return(x)
}

#' @rdname merge.ChromatinAssay
#' @export
#' @method merge ChromatinAssay
#' @importFrom GenomicRanges union findOverlaps
#' @importFrom Seurat RowMergeSparseMatrices
#' @importFrom S4Vectors subjectHits queryHits
merge.ChromatinAssay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL, # TODO add cell IDs, esp for fragment file paths
  merge.data = TRUE,
  ...
) {
  # check that genome is the same if both not NULL
  genome.1 <- genome(x = x)
  genome.2 <- genome(x = y)
  if ((!is.null(x = genome.1) &
       !is.null(x = genome.2)) &
      !(identical(x = genome.1, y = genome.2))) {
    warning("Genomes do not match, not merging ChromatinAssays")
    return(NULL)
  } else {
    genome.use <- SetIfNull(x = genome.1, y = genome.2)
  }
  # check the annotations slot
  annot.1 <- Annotation(object = x)
  annot.2 <- Annotation(object = y)
  if ((!is.null(x = annot.1) &
       !is.null(x = annot.2)) &
      !(identical(x = annot.1, y = annot.2))) {
    warning("Annotations do not match, keeping annotation from the
            first object")
    annot <- annot.1
  } else {
    annot <- SetIfNull(x = annot.1, y = annot.2)
  }
  # check fragments
  # TODO update with Fragments objects
  frag.1 <- GetAssayData(object = x, slot = "fragments")
  frag.2 <- GetAssayData(object = y, slot = "fragments")
  merged.frag <- c(frag.1, frag.2)
  valid.frags <- sapply(X = merged.frag, FUN = ValidFragments)
  if (!all(valid.frags)) {
    warning("Some fragment files are not valid or not indexed.
            Removing invalid files from merged ChromatinAssay")
    merged.frag <- merged.frag[valid.frags]
  }

  # find fraction overlap in granges in each direction
  overlaps <- findOverlaps(query = granges(x = x), subject = granges(x = y))
  percent.overlap.x <- length(
    x = unique(x = queryHits(x = overlaps))
  ) / length(granges(x = x)) * 100
  percent.overlap.y <- length(
    x = unique(x = subjectHits(x = overlaps))
  ) / length(granges(x = y)) * 100
  if (max(percent.overlap.x, percent.overlap.y) < 30) {
    warning("Few overlapping ranges between the assays to be merged:\n",
            "\t", round(x = percent.overlap.x, digits = 1), "% for x\n",
            "\t", round(x = percent.overlap.y, digits = 1), "% for y")
  }

  # merge matrix rows that intersect the same range
  # condense individual matrices first
  # should not typically have overlapping features in a single object
  condensed.a <- MergeIntersectingRows(
    mat.a = GetAssayData(object = x, slot = "counts"),
    ranges.a = granges(x = x),
    verbose = FALSE
  )
  condensed.b <- MergeIntersectingRows(
    mat.a = GetAssayData(object = y, slot = "counts"),
    ranges.a = granges(x = y),
    verbose = FALSE
  )
  # TODO this can be sped up significantly
  condensed <- MergeIntersectingRows(
    mat.a = condensed.a[[1]],
    mat.b = condensed.b[[1]],
    ranges.a = condensed.a[[2]],
    ranges.b = condensed.b[[2]],
    verbose = TRUE
  )
  # returns list: matrix A, ranges A, matrix B, ranges B

  # rename matrix rows in B that intersect A
  condensed[[3]] <- RenameIntersectingRows(
    mat.a = condensed[[1]],
    mat.b = condensed[[3]],
    ranges.a = condensed[[2]],
    ranges.b = condensed[[4]],
    verbose = TRUE
  )
  # should now have 1-1 intersection of ranges and matching row names

  # merge matrices
  # RowMergeSparseMatrices only exported in Seurat release Dec-2019 (3.1.2)
  merged.counts <- RowMergeSparseMatrices(
    mat1 = condensed.matrices[[1]],
    mat2 = condensed.matrices[[2]]
  )

  # perform same operation to ranges
  # TODO check this works
  grange.union <- union(x = condensed[[2]], y = condensed[[4]])

  # create new ChromatinAssay object
  # bias, motifs, positionEnrichment, metafeatures, scaledata and data not kept
  new.assay <- CreateChromatinAssayObject(
    counts = merged.counts,
    data = NULL,
    min.cells = 0,
    min.features = 0,
    max.cells = NULL,
    ranges = grange.union,
    motifs = NULL,
    fragments = merged.frag,
    genome = genome.use,
    annotation = annot,
    bias = NULL
  )
  return(new.assay)
}

#' @param i Which columns to retain
#' @param j Which rows to retain
#'
#' @rdname subset.Motif
#' @export
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
Annotation.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = "annotation"))
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

#' @rdname Fragments
#' @method Fragments ChromatinAssay
#' @export
Fragments.ChromatinAssay <- function(object, ...) {
  return(slot(object, name = "fragments"))
}

#' @param object A Seurat object or ChromatinAssay object
#' @importFrom Seurat DefaultAssay
#' @rdname Fragments
#' @method Fragments Seurat
#' @export
Fragments.Seurat <- function(object, ...) {
  assay <- DefaultAssay(object = object)
  return(Fragments(object = object[[assay]]))
}

#' @rdname Motifs
#' @method Motifs ChromatinAssay
#' @export
Motifs.ChromatinAssay <- function(object, ...) {
  return(slot(object = object, name = "motifs"))
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

#' @export
#' @method Annotation<- ChromatinAssay
"Annotation<-.ChromatinAssay" <- function(object, ..., value) {
  object <- SetAssayData(object = object, slot = "annotation", new.data = value)
  return(object)
}

#' @export
#' @importFrom Seurat DefaultAssay
#' @method Annotation<- Seurat
"Annotation<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  Annotation(object = object[[assay]]) <- value
  return(object)
}

#' @export
#' @method Fragments<- ChromatinAssay
"Fragments<-.ChromatinAssay" <- function(object, ..., value) {
  if (is.null(x = value)) {
    slot(object = object, name = "fragments") <- list()
    return(object)
  }
  if (inherits(x = value, what = "list")) {
    for (i in seq_along(along.with = value)) {
      object <- AddFragments(object = object, fragments = value[[i]])
    }
  } else {
    object <- AddFragments(object = object, fragments = value)
  }
  return(object)
}

#' @export
#' @importFrom Seurat DefaultAssay
#' @method Fragments<- Seurat
"Fragments<-.Seurat" <- function(object, ..., value) {
  assay <- DefaultAssay(object = object)
  Fragments(object = object[[assay]]) <- value
  return(object)
}

#' Add a single Fragment object to a ChromatinAssay
#' @param object A \code{\link{ChromatinAssay}} object
#' @param fragments A \code{\link{Fragment}} object
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

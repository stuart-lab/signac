#' @include generics.R
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass setClassUnion setMethod is slot slot<- new as
#' slotNames
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Signac
NULL

#' ChromatinAssay5 and GRangesAssay object constructors
#'
#' The [ChromatinAssay5-class()] and [GRangesAssay-class()]
#' classes extend [SeuratObject::Assay5()] objects for single-cell
#' chromatin data, with or without genomic ranges for each feature in the assay.
#' Use `CreateGRangesAssay` to construct a `GRangesAssay` object or
#' `CreateChromatinAssay5` to construct a `ChromatinAssay5` object.
#'
#' @rdname CreateGRangesAssay
#'
#' @param ranges A set of [GenomicRanges::GRanges()] corresponding to
#' the rows of the input matrix. If `NULL`, genomic ranges will be extracted
#' from the rownames of the input data matrix, and those rownames must be
#' formatted as `chromosome:start-end`.
#' @param ... Additional arguments passed to [CreateChromatinAssay5()]
#' and [SeuratObject::CreateAssay5Object()].
#'
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom Matrix rowSums colSums
#' @importFrom GenomicRanges isDisjoint
#' @importFrom S4Vectors mcols
#' @concept assay
#' @seealso [CreateChromatinAssay5()]
#' @return Returns a Seurat object
#'
#' @export
CreateGRangesAssay <- function(
  counts,
  data,
  ranges = NULL,
  ...
) {
  chrom.assay <- CreateChromatinAssay5(counts = counts, data = data, ...)
  if (!missing(x = counts)) {
    if (any(grepl("_", rownames(x = counts)))) {
      warning("Input matrix contains underscores in feature names. These are not
            allowed by Seurat; replacing underscores with '.'")
      rownames(x = counts) <- gsub(
        pattern = "_",
        replacement = ".",
        x = rownames(x = counts)
      )
    }
    features.keep <- rownames(x = counts) %in% rownames(x = chrom.assay)
  } else {
    if (any(grepl("_", rownames(x = data)))) {
      warning("Input matrix contains underscores in feature names. These are not
            allowed by Seurat; replacing underscores with '.'")
      rownames(x = data) <- gsub(
        pattern = "_",
        replacement = ".",
        x = rownames(x = data)
      )
    }
    features.keep <- rownames(x = data) %in% rownames(x = chrom.assay)
  }
  # subset ranges if there are features removed
  if (!is.null(x = ranges)) {
    ranges <- ranges[features.keep, ]
  }
  granges.assay <- as.GRangesAssay(x = chrom.assay, ranges = ranges, ...)
  return(granges.assay)
}

#' @param counts A two-dimensional matrix
#' @param data Optional pre-normalized count matrix
#' @param fragments Path to a tabix-indexed fragments file for the data
#' contained in the input matrix. If multiple fragment files are required,
#' you can add additional [Fragment()] object to the assay after it is
#' created using the [CreateFragmentObject()] and
#' [Fragments()] functions. Alternatively, a list of
#' [Fragment()] objects can be provided.
#' @param annotation A set of [GenomicRanges::GRanges()] containing
#' annotations for the genome used. It must have the following columns:
#' \itemize{
#'   \item{tx_id or transcript_id: Transcript ID}
#'   \item{gene_name: Gene name}
#'   \item{gene_id: Gene ID}
#'   \item{gene_biotype: Gene biotype (e.g. "protein_coding", "lincRNA")}
#'   \item{type: Annotation type (e.g. "exon", "gap")}
#' }
#' @param bias A Tn5 integration bias matrix
#' @param motifs A Motif object
#' @param links A named list of [InteractionSet::GInteractions()]
#' objects
#' @param region.aggregation A named list of [RegionAggregation-class]
#' objects.
#' @param validate.fragments Check that cells in the assay are present in the
#' fragment file.
#' @param verbose Display messages
#'
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom Matrix rowSums colSums
#' @importFrom S4Vectors mcols
#' @concept assay
#' @rdname CreateGRangesAssay
#' @aliases CreateChromatinAssay5
#'
#' @export
CreateChromatinAssay5 <- function(
  counts,
  data,
  fragments = NULL,
  annotation = NULL,
  motifs = NULL,
  links = NULL,
  bias = NULL,
  region.aggregation = NULL,
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
) {
  # Annotations
  if (!is.null(x = annotation) && !inherits(x = annotation, what = "GRanges")) {
    stop("Annotation must be a GRanges object.")
  }
  if (!is.null(x = annotation)) {
    if (!any(
      c("tx_id", "transcript_id") %in% colnames(x = mcols(x = annotation))
    )
    ) {
      stop(
        "Annotation must have transcript id ",
        "stored in `tx_id` or `transcript_id`."
      )
    }
    if (any(
      !c("gene_name", "gene_id", "gene_biotype", "type") %in%
        colnames(x = mcols(x = annotation))
    )) {
      stop(
        "Annotation must have `gene_name`, `gene_id`, ",
        "`gene_biotype` and `type`."
      )
    }
  }

  # CreateAssay5Object throws error if counts or data is missing
  # rather than NULL
  if (missing(x = counts)) {
    counts <- NULL
  }
  if (missing(x = data)) {
    data <- NULL
  }
  has_underscore <- any(grepl("_", rownames(x = counts)))
  if (has_underscore) {
    warning("Input matrix contains underscores in feature names. These are not
            allowed by Seurat; replacing underscores with '.'")
    rownames(x = counts) <- gsub("_", ".", rownames(x = counts))
  }

  seurat.assay <- CreateAssay5Object(counts = counts, data = data, ...)

  # Fragments
  if (inherits(x = fragments, what = "list")) {
    # check each object in the list is a fragment object
    # fragment list usually supplied when doing object merge,
    # so don't validate cells here, we can assume that was done in
    # individual object creation
    obj.class <- sapply(
      X = fragments, FUN = function(x) inherits(x = x, what = "Fragment2")
    )
    if (!all(obj.class)) {
      stop("All objects in fragments list must be Fragment-class objects")
    }
    frags <- lapply(
      X = fragments,
      FUN = AssignFragCellnames,
      cellnames = colnames(x = seurat.assay)
    )
    # subset to cells in the assay
    frags <- lapply(
      X = fragments,
      FUN = subset,
      cells = colnames(x = seurat.assay)
    )
  } else if (inherits(x = fragments, what = "Fragment2")) {
    # single Fragment object supplied
    frags <- AssignFragCellnames(
      fragments = fragments, cellnames = colnames(x = seurat.assay)
    )
    # subset to cells in the assay
    frags <- subset(x = frags, cells = colnames(x = seurat.assay))
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
          verbose = verbose
        )
      }
    }
  }

  # Motifs
  if (!is.null(x = motifs)) {
    if (!inherits(x = motifs, what = "Motif")) {
      stop("Provided motif object is not a Motif-class object")
    }
  }

  chrom.assay <- as.ChromatinAssay5(
    x = seurat.assay,
    fragments = frags,
    annotation = annotation,
    motifs = motifs,
    links = links,
    bias = bias,
    region.aggregation = region.aggregation
  )
  return(chrom.assay)
}


#' @rdname as.GRangesAssay
#' @method as.GRangesAssay ChromatinAssay
#' @importFrom InteractionSet GInteractions
#' @importFrom GenomicRanges start end
#' @importFrom IRanges IRanges
#' @importFrom Seqinfo seqnames
#' @export
#' @concept assay
as.GRangesAssay.ChromatinAssay <- function(x, ...) {
  # extract information
  frags <- x@fragments
  links <- x@links

  # update fragment objects
  if (length(x = frags) > 0) {
    for (i in seq_along(along.with = frags)) {
      frags[[i]] <- as.Fragment2(x = frags[[i]])
    }
  }

  if (length(x = x@positionEnrichment) > 0) {
    warning(
      "Cannot convert old positionEnrichment information, dropping ",
      length(x = x@positionEnrichment), " positionEnrichment matrices"
    )
  }

  # update links to list of GInteractions objects
  if (length(x = links) > 0) {
    start.gr <- GRanges(
      seqnames = seqnames(x = links),
      ranges = IRanges(start = start(x = links))
    )
    end.gr <- GRanges(
      seqnames = seqnames(x = links),
      ranges = IRanges(start = end(x = links))
    )
    gi <- GInteractions(start.gr, end.gr)
    mcols(x = gi) <- mcols(x = links)
    gi <- list("links" = gi)
  } else {
    gi <- NULL
  }

  # construct new assay object
  newobj <- as(object = x, Class = "Assay5")

  # update rownames
  rownames(x = newobj) <- as.character(x = x@ranges)

  newobj <- as.GRangesAssay(
    x = newobj,
    Class = "GRangesAssay",
    ranges = x@ranges,
    annotation = x@annotation,
    fragments = frags,
    bias = x@bias,
    motifs = x@motifs,
    links = gi,
    region.aggregation = NULL,
  )

  return(newobj)
}

#' @rdname as.Fragment2
#' @method as.Fragment2 Fragment
#' @export
#' @concept fragments
as.Fragment2.Fragment <- function(x, ...) {
  # convert from old Fragment class to new Fragment2 class

  # extract information from old object
  file.path <- x@path
  file.index <- paste0(file.path, ".tbi")
  hash <- x@hash
  cells <- x@cells

  # construct new object
  x <- new(
    Class = "Fragment2",
    file.path = file.path,
    file.index = file.index,
    hash = hash,
    cells = cells,
    seqlevels = NULL
  )
  return(x)
}

#' @param annotation A set of `GenomicRanges::GRanges()` containing gene
#' annotations for the genome used. It must have the following columns:
#' - `tx_id` or `transcript_id`: Transcript ID
#' - `gene_name`: Gene name
#' - `gene_id`: Gene ID
#' - `gene_biotype`: Gene biotype (e.g. "protein_coding", "lincRNA")
#' - `type`: Annotation type (e.g. "exon", "gap")
#' @param fragments Path to a tabix-indexed fragments file for the data
#' contained in the input matrix.
#' @param bias A Tn5 integration bias matrix vector
#' @param region.aggregation A named list of `RegionAggregation-class` objects.
#' @param motifs A `Motif-class` object
#' @param links A named list of `InteractionSet::GInteractions()` objects
#' @rdname as.GRangesAssay
#' @method as.GRangesAssay Assay5
#' @export
#' @concept assay
as.GRangesAssay.Assay5 <- function(
  x,
  annotation = NULL,
  fragments = NULL,
  bias = NULL,
  region.aggregation = NULL,
  ranges = NULL,
  motifs = NULL,
  links = NULL,
  ...
) {
  x <- as.ChromatinAssay5(
    x = x,
    annotation = annotation,
    fragments = fragments,
    bias = bias,
    motifs = motifs,
    links = links,
    region.aggregation = region.aggregation
  )
  x <- as.GRangesAssay(x = x, ranges = ranges)
  return(x)
}

#' @param ranges A GRanges object
#'
#' @rdname as.GRangesAssay
#' @export
#' @method as.GRangesAssay ChromatinAssay5
#' @concept assay
#'
as.GRangesAssay.ChromatinAssay5 <- function(
  x,
  ranges = NULL,
  ...
) {
  if (!is.null(x = ranges)) {
    if (length(x = ranges) != nrow(x = x)) {
      stop("Length of supplied genomic ranges does not match number
           of rows in matrix")
    }
  } else {
    ranges <- GRanges(rownames(x = x))
  }
  if (!isDisjoint(x = ranges)) {
    warning("Overlapping ranges supplied. Ranges should be non-overlapping.")
  }

  # re-assign row names of matrix so that it's a valid granges transformation
  new.rownames <- as.character(x = ranges)
  rownames(x = x) <- new.rownames

  new.assay <- as(object = x, Class = "GRangesAssay")
  new.assay <- SetAssayData(
    object = new.assay,
    layer = "ranges",
    new.data = ranges
  )
  return(new.assay)
}

#' @param annotation Genomic annotation. It must have the following columns:
#' \itemize{
#'   \item{tx_id or transcript_id: Transcript ID}
#'   \item{gene_name: Gene name}
#'   \item{gene_id: Gene ID}
#'   \item{gene_biotype: Gene biotype (e.g. "protein_coding", "lincRNA")}
#'   \item{type: Annotation type (e.g. "exon", "gap")}
#' }
#' @param fragments A list of [Fragment2-class] objects
#' @param bias Tn5 integration bias matrix
#' @param motifs A [Motif-class] object
#' @param links Genomic links TODO
#' @param region.aggregation A named list of [RegionAggregation-class]
#'
#' @rdname as.ChromatinAssay5
#' @export
#' @method as.ChromatinAssay5 Assay5
#' @concept assay
#'
as.ChromatinAssay5.Assay5 <- function(
  x,
  annotation = NULL,
  fragments = NULL,
  bias = NULL,
  motifs = NULL,
  links = NULL,
  region.aggregation = NULL,
  ...
) {
  new.assay <- as(object = x, Class = "ChromatinAssay5")
  if (!is.null(x = fragments)) {
    new.assay <- SetAssayData(
      object = new.assay,
      layer = "fragments",
      new.data = fragments
    )
  }
  if (!is.null(x = annotation)) {
    new.assay <- SetAssayData(
      object = new.assay,
      layer = "annotation",
      new.data = annotation
    )
  }
  if (!is.null(x = bias)) {
    new.assay <- SetAssayData(
      object = new.assay,
      layer = "bias",
      new.data = bias
    )
  }
  if (!is.null(x = region.aggregation)) {
    new.assay <- SetAssayData(
      object = new.assay,
      layer = "region.aggregation",
      new.data = region.aggregation
    )
  }
  if (!is.null(x = motifs)) {
    new.assay <- SetAssayData(
      object = new.assay,
      layer = "motifs",
      new.data = motifs
    )
  }
  if (!is.null(x = links)) {
    new.assay <- SetAssayData(
      object = new.assay,
      layer = "links",
      new.data = links
    )
  }
  return(new.assay)
}

#' Create a RegionAggregation object
#'
#' Create a [RegionAggregation-class()] object to store the insertions over
#' group of regions
#'
#' @param mat A cell-by-position matrix
#' @param regions A [GenomicRanges::granges()] object containing the
#' regions aggregated across.
#' @param upstream Integer denoting number of bases upstream of the
#' centered position that are stored in the matrix
#' @param downstream Integer denoting number of bases downstream of the
#' centered position that are stored in the matrix
#' @param name A name for the set of regions aggregated
#' @param cells A vector of cells where each element is the cell barcode
#' included in the region aggregation matrix. The order of cells in this
#' vector should correspond to the order of cells in the matrix. If `NULL`, cell
#' names will be extracted from the rownames of the matrix.
#' @param expected A vector containing expected number of Tn5 insertions per
#' position. If `NULL`, the expected value will be set uniformly to 1 for each
#' position.
#' @param verbose Display messages.
#'
#' @importFrom rlang is_integerish
#' @export
#' @return Returns a [RegionAggregation-class] object
#' @concept footprinting
CreateRegionAggregationObject <- function(
  mat,
  regions,
  upstream,
  downstream,
  name,
  cells = NULL,
  expected = NULL,
  verbose = TRUE
) {
  if (!inherits(x = mat, what = c("matrix", "CsparseMatrix"))) {
    stop(
      "data must be a matrix or sparse matrix. Supplied ",
      class(x = mat)
    )
  }
  if (any(dim(x = mat) == 0)) {
    stop("Provided matrix has zero dimensions")
  }
  if (!inherits(x = regions, what = "GRanges")) {
    stop("regions must be a GRanges object. Supplied ", class(x = regions))
  }
  # make upstream and downstream integer
  if (!is_integerish(x = upstream)) {
    stop("upstream value must be an integer")
  } else {
    upstream <- as.integer(x = upstream)
  }
  if (!is_integerish(x = downstream)) {
    stop("downstream value must be an integer")
  } else {
    downstream <- as.integer(x = downstream)
  }
  if (!is.null(x = expected)) {
    if (!is.numeric(x = expected)) {
      stop("expected insertions must be a numeric vector")
    }
    if (length(x = expected) != dim(x = mat)[2]) {
      stop(
        "expected insertion must match the number of positions in the matrix"
      )
    }
  } else {
    expected <- rep(x = 1, dim(x = mat)[2])
  }
  # if cells not given, create the cells from rownames of the matrix
  if (is.null(x = cells)) {
    cells <- rownames(x = mat)
    if (is.null(x = cells)) {
      stop(
        "cells information not provided, and the provided matrix has no",
        " row names"
      )
    }
  } else {
    cells <- as.character(x = cells)
    if (length(x = cells) != dim(x = mat)[1]) {
      stop(
        "Number of cells: (", length(x = cells),
        ") does not match number of matrix rows (", dim(x = mat)[1], ")"
      )
    }
  }

  # strip matrix dimnames
  dimnames(x = mat) <- NULL

  agg.obj <- new(
    Class = "RegionAggregation",
    matrix = mat,
    regions = regions,
    upstream = upstream,
    downstream = downstream,
    name = name,
    expected = expected,
    cells = cells
  )
  return(agg.obj)
}

#' Create a Fragment object
#'
#' Create a `Fragment` object to store fragment file information.
#' This object stores a 32-bit MD5 hash of the fragment file and the fragment
#' file index so that any changes to the files on-disk can be detected. A check
#' is also performed to ensure that the expected cells are present in the
#' fragment file.
#'
#' @param path A path to the fragment file.
#' @param index Path to the fragment file index. If NULL, the index file is
#' assumed to be in the same directory as the fragment file.
#' @param cells A named character vector containing cell barcodes contained in
#' the fragment file. This does not need to be all cells in the fragment file,
#' but there should be no cells in the vector that are not present in the
#' fragment file. A search of the file will be performed until at least one
#' fragment from each cell is found. If NULL, don't check for expected cells.
#'
#' Each element of the vector should be a cell barcode that appears in the
#' fragment file, and the name of each element should be the corresponding cell
#' name in the object.
#' @param seqlevels A named vector of sequence levels (eg, chromosome name)
#' where each element is the sequence name as it appears in the fragment file,
#' and the name of each element is the corresponding sequence name as stored in
#' the object. If NULL, the sequence names are assumed to be the same in the
#' fragment file and object.
#'
#' @param validate.fragments Check that expected cells are present in the
#' fragment file.
#' @param verbose Display messages
#' @param ... Additional arguments passed to `ValidateCells`
#'
#' @importFrom tools md5sum file_ext
#' @export
#' @concept fragments
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' cells <- colnames(x = atac_small)
#' names(x = cells) <- paste0("test_", cells)
#' frags <- CreateFragmentObject(
#'   path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5
#' )
CreateFragmentObject <- function(
  path,
  index = NULL,
  cells = NULL,
  seqlevels = NULL,
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
) {
  # check that file exists and is indexed
  # don't check if supplying remote file
  is.remote <- isRemote(x = path)
  if (is.remote) {
    validate.fragments <- FALSE
  }
  if (!file.exists(path) && !is.remote) {
    stop("Fragment file does not exist.")
  }
  index.file <- index %||% GetIndexFile(fragment = path, verbose = verbose)
  if (!file.exists(index.file) && !is.remote) {
    stop("Fragment file index does not exist.")
  }
  if (is.remote) {
    con <- gzcon(con = url(description = path))
  } else {
    con <- path
  }
  df <- readLines(con = con, n = 10000)
  for (i in df) {
    if (grepl(pattern = "^#", x = i)) {
      next
    } else {
      ncol_frag <- length(x = strsplit(x = i, split = "\t")[[1]])
      if (!(ncol_frag == 5 || ncol_frag == 6)) {
        # cellranger-atac v2.2 introduces strand column in fragment file
        stop("Incorrect number of columns found in fragment file")
      } else {
        break
      }
    }
  }
  if (!is.null(x = cells)) {
    if (is.null(names(x = cells))) {
      # assume cells are as they appear in the assay
      names(x = cells) <- cells
    }
  }
  if (!is.null(x = seqlevels)) {
    if (is.null(names(x = seqlevels))) {
      # assume seqlevels are as they appear in the assay
      names(x = seqlevels) <- seqlevels
    }
  }
  # compute hash of the file and index
  if (verbose) {
    message("Computing hash")
  }
  if (!is.remote) {
    path <- normalizePath(path = path, mustWork = TRUE)
    index.file <- normalizePath(path = index.file, mustWork = TRUE)
  }
  # will be NA if file remote
  hashes <- md5sum(files = c(path, index.file))
  # create object
  frags <- new(
    Class = "Fragment2",
    file.path = path,
    file.index = index.file,
    hash = unname(obj = hashes),
    cells = cells,
    seqlevels = seqlevels
  )
  # validate cells
  if (!is.null(x = cells) && validate.fragments) {
    if (ValidateCells(object = frags, verbose = verbose, ...)) {
      return(frags)
    } else {
      stop("Not all cells requested could be found in the fragment file.")
    }
  } else {
    return(frags)
  }
}

#' Create motif object
#'
#' Create a [Motif-class()] object.
#'
#' @param data A motif x region matrix
#' @param pwm A named list of position weight matrices or position frequency
#' matrices matching the motif names in `data`.
#' Can be of class PFMatrixList.
#' @param motif.names A named list of motif names. List element names
#' must match the names given in `pwm`. If NULL, use the names from the
#' list of position weight or position frequency matrices. This can be used to
#' set a alternative common name for the motif. If a PFMatrixList is passed to
#' `pwm`, it will pull the motif name from the PFMatrixList.
#' @param positions A [GenomicRanges::GRangesList()] object containing
#' exact positions of each motif.
#' @param meta.data A data.frame containing metadata
#' @export
#' @return Returns a [Motif-class()] object
#' @concept motifs
#' @examples
#' motif.matrix <- matrix(
#'   data = sample(c(0, 1),
#'     size = 100,
#'     replace = TRUE
#'   ),
#'   ncol = 5
#' )
#' rownames(motif.matrix) <- 1:nrow(motif.matrix)
#' motif <- CreateMotifObject(data = motif.matrix)
CreateMotifObject <- function(
  data = NULL,
  pwm = NULL,
  motif.names = NULL,
  positions = NULL,
  meta.data = NULL
) {
  data <- data %||% new(Class = "dgCMatrix")
  meta.data <- meta.data %||% data.frame()
  if (
    !(inherits(x = data, what = "matrix") ||
      inherits(x = data, what = "CsparseMatrix"))
  ) {
    stop(
      "Data must be matrix or sparse matrix class. Supplied ",
      class(x = data)
    )
  }
  if (inherits(x = data, what = "matrix")) {
    data <- as(Class = "CsparseMatrix", object = data)
  }
  if ((nrow(x = data) > 0) && (length(x = pwm) > 0)) {
    if (!all(names(x = pwm) == colnames(x = data))) {
      stop("Motif names in data matrix and PWM list are inconsistent")
    }
  }
  if ((nrow(x = data) > 0) && (nrow(x = meta.data) > 0)) {
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
  if (is.null(x = rownames(x = data))) {
    stop("Row names of data matrix cannot be NULL")
  }
  if (
    inherits(x = pwm, what = "PFMatrixList") ||
      inherits(x = pwm, what = "PWMatrixList")
  ) {
    pwm.converted <- lapply(X = as.list(x = pwm), FUN = PFMatrixToList)
    pwm <- lapply(X = pwm.converted, FUN = "[[", 1)
    motif.names <- lapply(X = pwm.converted, FUN = "[[", 2)
  }
  # ensure names are unique
  motif.id <- names(x = motif.names)
  mn.stash <- as.character(x = motif.names)
  mn.unique <- make.unique(names = as.character(motif.names))
  if (!identical(x = mn.stash, y = mn.unique)) {
    warning("Non-unique motif names supplied, making unique", immediate. = TRUE)
  }
  motif.names <- as.list(x = mn.unique)
  names(x = motif.names) <- motif.id
  pwm <- pwm %||% list()
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

# Update chromatin object
#
# Create a new \code{\link[SeuratObject]} with the cells only in the
# chromatin/expression assays. This is used for V5 objects, such as an
# extended reference object, that can have more cells in the whole object
# than are in the chromatin assay.
#
# @param object A \code{\link[SeuratObject]}
# @param chromatin.assay A list of the name(s) of the chromatin assay
# @param expression.assay The name of the expression assay
# @param features NULL or a list of features. If features is not null,
# the expression assay will be added to the object
# @return Returns a new \code{\link[SeuratObject]} that only contains
# the cells in the chromatin assay
#' @importFrom SeuratObject Assays CreateSeuratObject
UpdateChromatinObject <- function(
  object,
  chromatin.assay,
  expression.assay = NULL,
  features = NULL
) {
  op <- options(Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  lapply(X = chromatin.assay, FUN = function(x) {
    if (!(x %in% Assays(object))) {
      stop("The requested assay is not in the object.")
    }
  })
  # Create new seurat object
  new.object <- CreateSeuratObject(
    counts = object[[chromatin.assay[[1]]]],
    assay = chromatin.assay[[1]],
    meta.data = slot(object, name = "meta.data")
  )
  if (length(chromatin.assay) > 1) {
    for (i in 2:length(chromatin.assay)) {
      if (!identical(
        colnames(new.object[[chromatin.assay[[1]]]]),
        colnames(object[[chromatin.assay[[i]]]])
      )) {
        stop("All chromatin assays must have the same cells.")
      }
      new.object[[chromatin.assay[[i]]]] <- object[[chromatin.assay[[i]]]]
    }
  }
  # Add expression assay if applicable and if Seurat Object v5 is loaded
  if (!is.null(features)) {
    if (!is.null(expression.assay)) {
      if (utils::packageVersion("SeuratObject") >= package_version("4.9.9")) {
        if (!(expression.assay %in% Assays(object))) {
          stop("The requested assay is not in the object.")
        }
        if (!all(
          colnames(new.object) %in% colnames(object[[expression.assay]])
        )
        ) {
          stop("Chromatin and expression assays have different cells.")
        }
        # Convert BP Cells to sparse matrix
        for (i in SeuratObject::Layers(object[[expression.assay]])) {
          layer.data <- SeuratObject::LayerData(
            object = object,
            assay = expression.assay,
            layer = i
          )
          if (inherits(layer.data, what = "IterableMatrix")) {
            warning("Converting IterableMatrix to sparse dgCMatrix",
              call. = FALSE
            )
            SeuratObject::LayerData(
              object = object,
              assay = expression.assay,
              layer = i
            ) <- as(
              object = layer.data,
              Class = "dgCMatrix"
            )
          }
        }
        # Subset expression data if necessary
        if (!suppressWarnings(
          all(colnames(object[[expression.assay]]) == colnames(new.object))
        )
        ) {
          warning(
            "Subsetting expression assay to have ",
            "same cells as chromatin assay.",
            call. = FALSE
          )
          new.object[[expression.assay]] <- subset(
            x = object[[expression.assay]],
            cells = colnames(new.object)
          )
        } else {
          new.object[[expression.assay]] <- object[[expression.assay]]
        }
      } else {
        warning(
          "Cannot access layers if SeuratObject version is not ",
          "5.0.0 or greater. Please update SeuratObject to also ",
          "visualize expression data when your",
          "object has layers with different numbers of cells.",
          call. = FALSE,
          immediate. = TRUE
        )
      }
    }
  }
  return(new.object)
}

#' @importFrom SeuratObject GetAssayData
#' @method GetAssayData GRangesAssay
#' @importFrom lifecycle deprecated is_present
#' @export
#' @concept assay
GetAssayData.GRangesAssay <- function(
  object,
  layer = "data",
  assay = NULL,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    layer <- slot
  }
  if (
    layer %in% c("counts", "data", "scale.data", "meta.data", "misc", "key")
  ) {
    return(NextMethod())
  }
  if (!(layer %in% slotNames(x = object))) {
    stop(
      "layer must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  return(methods::slot(object = object, name = layer))
}


#' @importFrom SeuratObject GetAssayData
#' @method GetAssayData ChromatinAssay5
#' @importFrom lifecycle deprecated is_present
#' @export
#' @concept assay
GetAssayData.ChromatinAssay5 <- function(
  object,
  layer = "data",
  assay = NULL,
  slot = deprecated(),
  ...
) {
  if (is_present(arg = slot)) {
    layer <- slot
  }
  if (
    layer %in% c("counts", "data", "scale.data", "meta.data", "misc", "key")
  ) {
    return(NextMethod())
  }
  if (!(layer %in% slotNames(x = object))) {
    stop(
      "layer must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  return(methods::slot(object = object, name = layer))
}

#' Get Fragment object data
#'
#' Extract data from a [Fragment2-class] object
#'
#' @param object A [Fragment2-class] object
#' @param slot Information to pull from object
#' (file.path, index.path, hash, cells, seqlevels)
#' @export
#' @concept assay
GetFragmentData <- function(object, slot = "file.path") {
  return(methods::slot(object = object, name = slot))
}

#' @param slot Information to pull from object (data, pwm, meta.data)
#' @rdname GetMotifData
#' @method GetMotifData Motif
#' @export
#' @concept motifs
#' @examples
#' motif.obj <- SeuratObject::GetAssayData(
#'   object = atac_small[["peaks"]], slot = "motifs"
#' )
#' GetMotifData(object = motif.obj)
GetMotifData.Motif <- function(object, slot = "data", ...) {
  return(methods::slot(object = object, name = slot))
}

#' @importFrom SeuratObject GetAssayData
#' @rdname GetMotifData
#' @concept motifs
#' @method GetMotifData ChromatinAssay5
#' @export
GetMotifData.ChromatinAssay5 <- function(object, slot = "data", ...) {
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
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @examples
#' GetMotifData(object = atac_small)
GetMotifData.Seurat <- function(object, assay = NULL, slot = "data", ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(GetMotifData(
    object = object[[assay]],
    slot = slot,
    ...
  ))
}

#' @importFrom SeuratObject RenameCells GetAssayData
#' @concept assay
#' @method RenameCells ChromatinAssay5
#' @export
RenameCells.ChromatinAssay5 <- function(object, new.names = NULL, ...) {
  # name of each element is the existing cell name
  # element itself is the corresponding new name
  # new.names vector names will be set in the Seurat method. Need to set here
  # in case RenameCells is called directly on a ChromatinAssay5 object
  if (is.null(x = names(x = new.names))) {
    if (length(x = new.names) != ncol(x = object)) {
      stop("Insufficient names supplied to rename cells")
    }
    names(x = new.names) <- colnames(x = object)
  }

  # there's nothing that needs to be renamed in the GRangesAssay class
  # rename cells in the parental class
  object <- NextMethod()

  # rename cells in fragment objects
  frags <- Fragments(object = object)
  for (i in seq_along(along.with = frags)) {
    frags[[i]] <- RenameCells(object = frags[[i]], new.names = new.names)
  }
  Fragments(object = object) <- NULL
  Fragments(object = object) <- frags

  # rename cells in RegionAggregation objects
  region.aggr <- RegionAggr(object = object)
  for (i in seq_along(along.with = region.aggr)) {
    region.aggr[[i]] <- RenameCells(
      object = region.aggr[[i]],
      new.names = new.names
    )
  }
  RegionAggr(object = object) <- NULL
  RegionAggr(object = object) <- region.aggr

  return(object)
}

#' @importFrom SeuratObject RenameCells
#' @concept fragments
#' @method RenameCells Fragment2
#' @export
RenameCells.Fragment2 <- function(object, new.names, ...) {
  cells <- GetFragmentData(object = object, slot = "cells")
  if (is.null(x = cells)) {
    stop(
      "Cannot rename cells in Fragment object ",
      "with no cell information stored"
    )
  }
  if (is.null(x = names(x = new.names))) {
    if (length(x = new.names) != length(x = cells)) {
      stop("Insufficient names supplied to rename cells")
    }
    names(x = new.names) <- cells
  }
  cells <- cells[names(x = new.names)]
  names(x = cells) <- new.names[names(x = cells)]
  slot(object = object, name = "cells") <- cells
  return(object)
}

#' @importFrom SeuratObject RenameCells
#' @importFrom methods slot "slot<-"
#' @concept footprinting
#' @method RenameCells RegionAggregation
#' @export
RenameCells.RegionAggregation <- function(object, new.names, ...) {
  # name of each element is the existing cell name
  # element itself is the corresponding new name

  cells <- slot(object = object, name = "cells")
  if (is.null(x = cells)) {
    # invalid object
    stop(
      "Cannot rename cells in RegionAggregation object ",
      "with no cell information stored"
    )
  }
  if (is.null(x = names(x = new.names))) {
    if (length(x = new.names) != length(x = cells)) {
      stop("Insufficient names supplied to rename cells")
    }
    names(x = new.names) <- cells
  }
  new_cells <- unname(new.names[cells])
  slot(object = object, name = "cells") <- new_cells
  return(object)
}

#' @importFrom SeuratObject SetAssayData
#' @importFrom Seqinfo genome Seqinfo
#' @importFrom lifecycle deprecated is_present
#' @importFrom S4Vectors mcols
#' @method SetAssayData GRangesAssay
#' @concept assay
#' @export
SetAssayData.GRangesAssay <- function(
  object,
  layer,
  new.data,
  slot = deprecated(),
  ...
) {
  if (layer %in% c(
    "counts", "data", "scale.data", "meta.data", "misc", "key",
    "fragments", "annotation", "bias", "region.aggregation", "motifs", "links"
  )) {
    return(NextMethod())
  }
  if (is_present(arg = slot)) {
    layer <- slot
  }
  if (!(layer %in% slotNames(x = object))) {
    stop(
      "layer must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  if (layer == "ranges") {
    if (!is(object = new.data, class2 = "GRanges")) {
      stop("Must provide a GRanges object")
    } else if (length(x = new.data) != nrow(x = object)) {
      stop("Number of ranges provided is not equal to the number
           of features in the assay")
    }
    methods::slot(object = object, name = layer) <- new.data
  } else {
    stop("SetAssayData not implemented for ", layer)
  }
  return(object)
}

#' @importFrom SeuratObject SetAssayData
#' @importFrom Seqinfo genome Seqinfo
#' @importFrom lifecycle deprecated is_present
#' @importFrom S4Vectors mcols
#' @method SetAssayData ChromatinAssay5
#' @concept assay
#' @export
SetAssayData.ChromatinAssay5 <- function(
  object,
  layer,
  new.data,
  slot = deprecated(),
  ...
) {
  if (
    layer %in% c("counts", "data", "scale.data", "meta.data", "misc", "key")
  ) {
    return(NextMethod())
  }
  if (is_present(arg = slot)) {
    layer <- slot
  }
  if (!(layer %in% slotNames(x = object))) {
    stop(
      "layer must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }

  if (layer == "fragments") {
    # this will overwrite the slot with new data
    # Fragments<- will append new fragment objects to the list
    # Fragments function is the recommended user-facing function
    if (inherits(x = new.data, what = "list")) {
      # check that it's a list containing fragment class objects
      for (i in seq_along(along.with = new.data)) {
        if (!inherits(x = new.data[[i]], what = "Fragment2")) {
          stop("New data is not a Fragment object")
        }
      }
    } else if (inherits(x = new.data, what = "Fragment2")) {
      # single fragment object
      new.data <- list(new.data)
    }
    frag.list <- GetAssayData(object = object, layer = "fragments")
    if (length(x = frag.list) != 0) {
      warning("Overwriting existing fragment objects")
    }

    # resolve any duplicated fragment file paths
    all.path <- lapply(X = new.data, FUN = GetFragmentData, slot = "file.path")
    unique.paths <- unique(all.path)
    if (length(x = unique.paths) != length(x = new.data)) {
      # consolidate duplicated fragment paths
      new.data <- lapply(unique.paths, function(p) {
        idx <- which(all.path == p)
        objs <- new.data[idx]
        # merge cell vectors
        all.cells <- do.call(c, lapply(objs, GetFragmentData, slot = "cells"))
        duplicate.cells <- duplicated(x = all.cells)
        all.cells <- all.cells[!duplicate.cells]
        objs[[1]]@cells <- all.cells
        objs[[1]]
      })
    }
    methods::slot(object = object, name = "fragments") <- new.data
  } else if (layer == "annotation") {
    if (!is(object = new.data, class2 = "GRanges")) {
      stop("Must provide a GRanges object")
    }
    if (!any(
      c("tx_id", "transcript_id") %in% colnames(x = mcols(x = new.data))
    )
    ) {
      stop(
        "Annotation must have transcript id stored ",
        "in `tx_id` or `transcript_id`."
      )
    }
    if (any(
      !c("gene_name", "gene_id", "gene_biotype", "type") %in%
        colnames(x = mcols(x = new.data))
    )) {
      stop(
        "Annotation must have `gene_name`, `gene_id`, ",
        "`gene_biotype` and `type`."
      )
    }
    if (!"tx_id" %in% colnames(x = mcols(x = new.data))) {
      new.data$tx_id <- new.data$transcript_id
    }
    methods::slot(object = object, name = layer) <- new.data
  } else if (layer == "bias") {
    if (!is(object = new.data, class2 = "vector")) {
      stop("Bias must be provided as a vector")
    }
    methods::slot(object = object, name = layer) <- new.data
  } else if (layer == "region.aggregation") {
    # pull overwrite from ... , only interpret inside layer==region.aggregation
    dots <- list(...)
    overwrite <- dots$overwrite %||% FALSE # default FALSE if not provided
    if (is.null(x = new.data) || length(x = new.data) == 0) {
      # overwrite with empty list
      methods::slot(object = object, name = layer) <- list()
      return(object)
    }
    if (inherits(x = new.data, what = "list")) {
      # check if its a list containing RegionAggregation class objects
      for (i in seq_along(new.data)) {
        if (!inherits(x = new.data[[i]], what = "RegionAggregation")) {
          stop("New data is not a RegionAggregation object")
        }
      }
    } else if (inherits(x = new.data, what = "RegionAggregation")) {
      # single RegionAggregation object
      new.data <- list(new.data)
    }
    # get existing RegionAggregation object
    agg.list <- GetAssayData(object = object, layer = layer)
    if (length(x = agg.list) == 0) {
      # nothing exists yet -> assign directly
      new.data <- MergeRegionAggregation(new.data)
      methods::slot(object, "region.aggregation") <- new.data
      return(object)
    }

    for (i in seq_along(new.data)) {
      new.agg <- new.data[[i]]
      new.cells <- new.agg@cells
      merged <- FALSE
      # compare against same-name objects
      same.name.idx <- which(vapply(
        agg.list, function(x) identical(x@name, new.agg@name), logical(1)
      ))

      if (length(same.name.idx) > 0) {
        if (overwrite) {
          # remove every exisiting RegAggr object with the same feature name
          warning(sprintf(
            "Overwriting RegionAggregation for '%s' ",
            new.agg@name
          ), call. = FALSE)
          agg.list <- agg.list[-same.name.idx]
          agg.list <- append(agg.list, new.agg)
          merged <- TRUE
        } else { # skip the new cells that already exists in the old object
          for (j in same.name.idx) {
            old.agg <- agg.list[[j]]
            overlap.cells <- intersect(old.agg@cells, new.agg@cells)
            new.cells <- setdiff(new.agg@cells, old.agg@cells)
            if (length(overlap.cells) > 0) {
              warning(sprintf(
                paste0(
                  "RegionAggregation '%s' already exists for %d cells and will",
                  " not be recomputed. \n",
                  "Set overwrite=TRUE to replace the existing ",
                  "RegionAggregation, or supply a different name to store it ",
                  "separately"
                ),
                new.agg@name,
                length(overlap.cells)
              ), call. = FALSE)
            }
            if (length(new.cells) > 0) {
              compatible <- IsCompatibleRegionAggregation(old.agg, new.agg)

              if (compatible) {
                # concatenate the matrix and the cells vector
                old.agg@matrix <- rbind(old.agg@matrix, new.agg@matrix)
                old.agg@cells <- c(old.agg@cells, new.agg@cells)
                agg.list[[j]] <- old.agg
                merged <- TRUE
                break
              }
            }
          }
        }
      }
      if (!merged) {
        new.sub <- subset(new.agg, cells = new.cells)
        if (!is.null(new.sub)) {
          agg.list <- append(agg.list, list(new.sub))
        }
      }
    }
    methods::slot(object = object, layer) <- agg.list
  } else if (layer == "motifs") {
    if (is.null(x = new.data)) {
      methods::slot(object = object, name = layer) <- NULL
      return(object)
    }
    if (!inherits(x = new.data, what = "Motif")) {
      stop("Must provide a Motif class object")
    }
    if (!is.null(x = GetAssayData(object = object, layer = "motifs"))) {
      warning("Overwriting motif information")
    }
    methods::slot(object = object, name = layer) <- new.data
  } else if (layer == "links") {
    if (is.null(x = new.data)) {
      # empty links list
      methods::slot(object = object, name = layer) <- list()
    } else {
      # ensure object is a list of GInteractions objects
      if (inherits(x = new.data, what = "list")) {
        # check that it's a list containing GInteraction objects
        for (i in seq_along(new.data)) {
          if (!inherits(x = new.data[[i]], what = "GInteractions")) {
            stop("New data is not a GInteractions object")
          }
        }
        # if passing a list, replace existing data
        methods::slot(object = object, name = layer) <- new.data
      } else {
        if (!inherits(x = new.data, what = "GInteractions")) {
          stop("New data is not a GInteractions object")
        } else {
          # append to list of GInteractions using key
          args <- list(...)
          key <- args[["key"]]
          if (is.null(x = key)) {
            stop("Key not supplied")
          }
          existing <- Links(object = object)
          if (key %in% names(x = existing)) {
            warning("Overwriting existing links key: ", key)
          }
          existing[[key]] <- new.data
          methods::slot(object = object, name = layer) <- existing
        }
      }
    }
  } else {
    stop("SetAssayData not implemented for ", layer)
  }
  return(object)
}
#' @rdname SetMotifData
#' @method SetMotifData Motif
#' @export
#' @concept motifs
#' @examples
#' motif.obj <- SeuratObject::GetAssayData(
#'   object = atac_small[["peaks"]], slot = "motifs"
#' )
#' SetMotifData(object = motif.obj, slot = "data", new.data = matrix(1:2))
SetMotifData.Motif <- function(object, slot, new.data, ...) {
  if (!(slot %in% slotNames(x = object))) {
    stop("slot must be one of ",
      paste(slotNames(x = object), collapse = ", "),
      call. = FALSE
    )
  }
  if (slot == "data") {
    if (inherits(x = new.data, what = "matrix")) {
      new.data <- as(Class = "CsparseMatrix", object = new.data)
    }
  }
  # TODO check that new data is compatible with existing slots
  # rownames of data must match rownames of meta.data and names of pwm, if not
  # empty
  methods::slot(object = object, name = slot) <- new.data
  return(object)
}

#' @param new.data motif matrix to add. Should be matrix or sparse matrix class
#' @param slot Name of slot to use
#' @importFrom SeuratObject GetAssayData SetAssayData
#' @rdname SetMotifData
#' @export
#' @concept motifs
#' @examples
#' new.data <- matrix(sample(c(0, 1),
#'   size = nrow(atac_small[["peaks"]]) * 10,
#'   replace = TRUE
#' ), nrow = nrow(atac_small[["peaks"]]))
#' rownames(new.data) <- rownames(atac_small[["peaks"]])
#' SetMotifData(
#'   object = atac_small[["peaks"]], slot = "data", new.data = new.data
#' )
#' @method SetMotifData ChromatinAssay5
SetMotifData.ChromatinAssay5 <- function(object, slot, new.data, ...) {
  if (slot == "data") {
    if (
      !(inherits(x = new.data, what = "matrix") ||
        inherits(x = new.data, what = "CsparseMatrix"))
    ) {
      stop(
        "Data must be matrix or sparse matrix class. Supplied ",
        class(x = new.data)
      )
    }
    if (inherits(x = new.data, what = "matrix")) {
      new.data <- as(Class = "CsparseMatrix", object = new.data)
    }
  }
  motif.obj <- GetAssayData(object = object, layer = "motifs")
  if (is.null(x = motif.obj)) {
    stop("Motif object not present in assay")
  } else {
    motif.obj <- SetMotifData(
      object = motif.obj, slot = slot, new.data = new.data
    )
    object <- SetAssayData(
      object = object, layer = "motifs", new.data = motif.obj
    )
    return(object)
  }
}

#' @param assay Name of assay whose data should be set
#' @rdname SetMotifData
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @method SetMotifData Seurat
#' @concept motifs
#' @examples
#' motif.matrix <- GetMotifData(object = atac_small)
#' SetMotifData(
#'   object = atac_small,
#'   assay = "peaks",
#'   slot = "data",
#'   new.data = motif.matrix
#' )
SetMotifData.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  object[[assay]] <- SetMotifData(object = object[[assay]], ...)
  return(object)
}

#' Subset a Motif object
#'
#' Returns a subset of a [Motif-class()] object.
#'
#' @param x A Motif object
#' @param features Which features to retain
#' @param motifs Which motifs to retain
#' @param ... Arguments passed to other methods
#'
#' @method subset Motif
#'
#' @seealso [base::subset()]
#' @return Returns a subsetted [Motif-class] object
#' @export
#' @concept motifs
#' @examples
#' motif.obj <- SeuratObject::GetAssayData(
#'   object = atac_small[["peaks"]], layer = "motifs"
#' )
#' subset(x = motif.obj, features = head(rownames(motif.obj), 10))
subset.Motif <- function(x, features = NULL, motifs = NULL, ...) {
  features <- features %||% rownames(x = x)
  motifs <- motifs %||% colnames(x = x)
  new.data <- GetMotifData(
    object = x, slot = "data"
  )[features, motifs, drop = FALSE]
  new.pwm <- GetMotifData(object = x, slot = "pwm")[motifs]
  new.names <- GetMotifData(object = x, slot = "motif.names")[motifs]
  new.meta <- GetMotifData(
    object = x, slot = "meta.data"
  )[motifs, , drop = FALSE]
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
#' @importClassesFrom SeuratObject Assay5
#' @concept assay
#' @method subset GRangesAssay
subset.GRangesAssay <- function(
  x,
  features = NULL,
  cells = NULL,
  ...
) {
  if (inherits(x = features, what = "Rle")) {
    features <- as.vector(x = features)
  }
  if (inherits(x = cells, what = "Rle")) {
    cells <- as.vector(x = cells)
  }

  # subset genomic ranges
  ranges.keep <- granges(x = x)
  if (!is.null(x = features)) {
    if (is.logical(x = features)) {
      idx.keep <- features
    } else {
      idx.keep <- rownames(x = x) %in% features
    }
    ranges.keep <- ranges.keep[idx.keep]
  }

  # subset elements in the parent classes assay
  x <- NextMethod()

  # subset granges
  x <- SetAssayData(object = x, layer = "ranges", new.data = ranges.keep)

  return(x)
}

#' @export
#' @importClassesFrom SeuratObject Assay5
#' @concept assay
#' @method subset ChromatinAssay5
subset.ChromatinAssay5 <- function(
  x,
  features = NULL,
  cells = NULL,
  ...
) {
  if (inherits(x = features, what = "Rle")) {
    features <- as.vector(x = features)
  }
  if (inherits(x = cells, what = "Rle")) {
    cells <- as.vector(x = cells)
  }

  if (is.logical(x = features)) {
    if (length(x = features) != nrow(x = x)) {
      stop("Incorrect number of logical values provided to subset features")
    } else {
      features <- Features(x = x)[features]
    }
  }

  if (is.logical(x = cells)) {
    if (length(x = cells) != ncol(x = x)) {
      stop("Incorrect number of logical values provided to subset cells")
    } else {
      cells <- Cells(x = x)[cells]
    }
  }


  # subset elements in the standard assay
  x <- NextMethod()

  # subset cells in RegionAggregation objects
  if (!is.null(x = cells)) {
    ragg <- RegionAggr(object = x)
    if (length(x = ragg) > 0) {
      # list of region aggregation objects
      ragg <- lapply(X = ragg, FUN = subset, cells = cells)
      ragg <- Filter(f = Negate(f = is.null), x = ragg)
    }
    RegionAggr(object = x, overwrite = TRUE) <- ragg
  }

  # subset cells in Fragments objects
  if (!is.null(x = cells)) {
    frags <- Fragments(object = x)
    Fragments(object = x) <- NULL
    for (i in seq_along(along.with = frags)) {
      frags[[i]] <- subset(x = frags[[i]], cells = cells)
    }
    Fragments(object = x) <- frags
  }

  # subset motifs
  motifs <- Motifs(object = x)
  if (!is.null(x = motifs)) {
    motifs <- subset(x = motifs, features = features)
  }
  Motifs(object = x) <- motifs

  return(x)
}

#' Subset a RegionAggregation object
#'
#' @param x A [RegionAggregation-class] object
#' @param cells Vector of cells to retain
#' @param ... Arguments passed to other methods
#'
#' @method subset RegionAggregation
#' @importFrom fastmatch fmatch
#' @importFrom methods slot "slot<-"
#' @rdname RegionAggregation-class
#' @export
#' @concept footprinting
subset.RegionAggregation <- function(
  x,
  cells = NULL,
  ...
) {
  if (is.null(x = cells)) {
    return(x)
  }
  mat <- slot(object = x, name = "matrix")
  agg.cells <- slot(object = x, name = "cells")

  keep <- fmatch(x = agg.cells, table = cells, nomatch = 0L) > 0
  if (!any(keep)) {
    return(NULL)
  }

  slot(object = x, name = "matrix") <- mat[keep, , drop = FALSE]
  slot(object = x, name = "cells") <- agg.cells[keep]

  return(x)
}

#' Subset a Fragment object
#'
#' Returns a subset of a [Fragment()] object.
#'
#' @param x A Fragment object
#' @param cells Vector of cells to retain
#' @param ... Arguments passed to other methods
#'
#' @method subset Fragment2
#'
#' @importFrom fastmatch fmatch
#' @seealso [base::subset()]
#' @return Returns a subsetted [Fragment()] object
#' @export
#' @concept fragments
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' cells <- colnames(x = atac_small)
#' names(x = cells) <- paste0("test_", cells)
#' frags <- CreateFragmentObject(
#'   path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5
#' )
#' subset(frags, head(names(cells)))
subset.Fragment2 <- function(
  x,
  cells = NULL,
  ...
) {
  frag.cells <- GetFragmentData(object = x, slot = "cells")
  keep <- fmatch(x = names(x = frag.cells), table = cells, nomatch = 0L) > 0
  if (!any(keep)) {
    return(NULL)
  }
  slot(object = x, name = "cells") <- frag.cells[keep]
  return(x)
}

#' @export
#' @concept assay
#' @method merge GRangesAssay
#' @importFrom GenomicRanges union findOverlaps
#' @importFrom SeuratObject RowMergeSparseMatrices Key Key<-
#' @importFrom S4Vectors subjectHits queryHits mcols
merge.GRangesAssay <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  ...
) {
  # call NextMethod to get merged ChromatinAssay5 object
  merged <- NextMethod()

  # need to do all operations over a list of assays
  assays <- c(x, y)

  # check that all features are equal
  all.features <- lapply(X = assays, FUN = rownames)
  all.features <- table(do.call(what = c, args = all.features))
  all.identical <- all(all.features == length(x = assays))

  if (!all.identical) {
    warning("Merging objects with different feature sets")
  }

  # if any are standard Assay class, coerce all to Assay and run merge
  isGranges <- sapply(
    X = assays, FUN = function(x) inherits(x = x, what = "GRangesAssay")
  )
  if (!all(isGranges)) {
    warning(
      "Some assays are not GRangesAssay class, ",
      "coercing GRangesAssay to ChromatinAssay5"
    )
    return(merged)
  } else {
    # list of GRangesAssay

    # merge granges across all objects
    granges.all <- GRanges(rownames(x = merged))
    merged <- as.GRangesAssay(x = merged, ranges = granges.all)
    return(merged)
  }
}

# Condense a list of RegionAggregation objects
#
# Takes a list of RegionAggregation objects and merges compatible
# objects that share the same feature name. Compatibility is defined by
# identical upstream/downstream extension, regions, and expected insertion
# values. Objects with different feature names are never merged.
#
# The function performs strict invariant checks and will error if:
# \itemize{
#  \item RegionAggregation objects with the same feature name have different
#        region widths.
#  \item RegionAggregation objects with the same feature name contain
#        overlapping cell barcodes.
# }
#
# @param x A list of RegionAggregation objects
#
# @return A list of RegionAggregation objects, where compatible objects
#  have been merged.
#
# @details
# Compatibility between two RegionAggregation objects is determined by
# \code{IsCompatibleRegionAggregation()}
#
MergeRegionAggregation <- function(
  x = NULL
) {
  stopifnot(is.list(x))
  stopifnot(all(vapply(x, inherits, logical(1), "RegionAggregation")))

  # extract feature names
  feature.names <- vapply(x, FUN = function(x) x@name, FUN.VALUE = character(1))

  # group by feature names
  grouped.agg <- split(x = x, f = feature.names)

  condensed.list <- lapply(X = grouped.agg, FUN = function(aggs) {
    if (length(x = aggs) == 1) {
      # only one agg obj with this feature name
      # nothing to do
      out <- aggs[[1]]
    } else {
      ## single motif width per feature
      w <- unique(
        x = unlist(
          x = lapply(
            X = aggs,
            FUN = function(x) width(x = x@regions)
          )
        )
      )
      if (length(x = w) != 1) {
        stop(
          "RegionAggregation for ",
          aggs[[1]]@name,
          " have different widths ",
          w
        )
      }
      ## feature should not have duplicated cell barcodes
      cells.all <- unlist(
        x = lapply(
          X = aggs,
          FUN = function(x) x@cells
        ),
        use.names = FALSE
      )
      dup.cells <- cells.all[duplicated(x = cells.all)]
      if (length(x = dup.cells) > 0) {
        stop(
          "Duplicated cell barcodes across RegionAggregation objects for ",
          aggs[[1]]@name
        )
      }

      ## opportunistic merging
      merged.obj.list <- list()
      for (i in seq_along(along.with = aggs)) {
        if (length(x = merged.obj.list) == 0) {
          merged.obj.list <- list(aggs[[i]])
        } else {
          merged <- FALSE
          for (w in seq_along(along.with = merged.obj.list)) {
            if (IsCompatibleRegionAggregation(
              x = aggs[[i]], y = merged.obj.list[[w]]
            )
            ) {
              # replace merged.obj.list[w] with the merged
              merged.obj.list[[w]]@matrix <- rbind(
                merged.obj.list[[w]]@matrix, aggs[[i]]@matrix
              )
              merged.obj.list[[w]]@cells <- c(
                merged.obj.list[[w]]@cells,
                aggs[[i]]@cells
              )
              merged <- TRUE
              break # stop looping through merged.obj.list
            }
          }
          if (!merged) {
            merged.obj.list <- append(x = merged.obj.list, values = aggs[[i]])
          }
        }
      }
      out <- merged.obj.list
    }
    out
  })
  # return as a flatten list
  condensed.list <- unname(obj = do.call(what = c, args = condensed.list))
  return(condensed.list)
}

# Check compatibility of two RegionAggregation objects
#
# Determines whether two [RegionAggregation-class] objects can be merged.
# Compatibility id  defined as having identical non-cell-specific slots,
# including:
# \itemize{
#    \item \code{upstream} and \code{downstream} extensions
#    \item \code{expected} insertion vector
#    \item \code{regions} (\code{Granges} object)
# }
#
# @returns `TRUE` if `x` and `y` are compatible for merging `FALSE` otherwise.
IsCompatibleRegionAggregation <- function(
  x = NULL,
  y = NULL
) {
  stopifnot(inherits(x, "RegionAggregation"))
  stopifnot(inherits(y, "RegionAggregation"))

  identical(x@upstream, y@upstream) &&
    identical(x@downstream, y@downstream) &&
    identical(x@expected, y@expected) &&
    identical(x@regions, y@regions)
}

#' @export
#' @concept assay
#' @method merge ChromatinAssay5
merge.ChromatinAssay5 <- function(
  x = NULL,
  y = NULL,
  add.cell.ids = NULL,
  ...
) {
  # call NextMethod to get merged standard assay object
  merged <- NextMethod()

  # need to do all operations over a list of assays
  assays <- c(x, y)

  # if any are standard Assay class, coerce all to Assay and run merge
  isChromatin <- sapply(
    X = assays, FUN = function(x) inherits(x = x, what = "ChromatinAssay5")
  )
  if (!all(isChromatin)) {
    warning(
      "Some assays are not ChromatinAssay5 class, ",
      "coercing ChromatinAssays to standard Assay"
    )
    assays <- sapply(
      X = assays, FUN = function(x) as(object = x, Class = "Assay5")
    )
    new.assay <- merge(
      x = assays[[1]], y = assays[2:length(x = assays)], ...
    )
    return(new.assay)
  } else {
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
          new.names = paste(
            add.cell.ids[i], colnames(x = assays[[i]]),
            sep = "_"
          )
        )
      }
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
            Removing invalid files from merged ChromatinAssay5")
      all.frag <- all.frag[valid.frags]
    }

    # merge region.aggregations
    all.agg <- lapply(X = assays, FUN = function(x) {
      RegionAggr(object = x)
    })
    # flatten list-of-lists
    all.agg <- Reduce(f = c, x = all.agg)

    # create new ChromatinAssay5 object
    # bias, motifs, links, metafeatures not kept
    merged <- as.ChromatinAssay5(
      x = merged,
      annotation = annot.use,
      fragments = all.frag,
      bias = NULL,
      region.aggregation = all.agg
    )
    return(merged)
  }
}

#' @param i Which columns to retain
#' @param j Which rows to retain
#'
#' @rdname subset.Motif
#' @export
#' @concept motifs
#' @method [ Motif
#' @examples
#' motif.obj <- SeuratObject::GetAssayData(
#'   object = atac_small, assay = "peaks", layer = "motifs"
#' )
#' motif.obj[1:10, 1:10]
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
  signature = "Fragment2",
  definition = function(object) {
    sl <- seqlevels(x = object)
    if (is.null(x = sl)) {
      sl.text <- "Unknown seqlevels\n"
    } else {
      sl.text <- paste0(
        "Seqlevels: ",
        paste(head(x = sl), collapse = ", "),
        " ...\n"
      )
    }
    cat(
      "A Fragment v2 object for",
      length(x = slot(object = object, name = "cells")),
      "cells\n",
      sl.text
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
  signature = "GRangesAssay",
  definition = function(object) {
    cat(
      "GRangesAssay data with",
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
      "Annotation present:",
      ifelse(
        test = is.null(x = Annotation(object = object)), yes = FALSE, no = TRUE
      ),
      "\n"
    )
    cat(
      "Fragment files:",
      length(x = Fragments(object = object)),
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
      "Links present:", length(x = Links(object = object)),
      "\n"
    )
    ragg <- RegionAggr(object = object)
    if (length(x = ragg) == 0) {
      cat("Region aggregation matrices: 0\n")
    } else {
      ragstr <- ifelse(length(x = ragg) > 6, "...", "")
      cat(
        length(x = ragg),
        "region aggregation",
        ifelse(length(x = ragg) == 1, "matrix:", "matrices:"),
        RegionAggNames(object = object),
        ragstr, "\n"
      )
    }
  }
)

setMethod(
  f = "show",
  signature = "ChromatinAssay5",
  definition = function(object) {
    cat(
      "ChromatinAssay5 data with",
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
      "Annotation present:",
      ifelse(
        test = is.null(x = Annotation(object = object)), yes = FALSE, no = TRUE
      ),
      "\n"
    )
    cat(
      "Fragment files:",
      length(x = Fragments(object = object)),
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
      "Links present:", length(x = Links(object = object)),
      "\n"
    )
    ragg <- RegionAggr(object = object)
    if (length(x = ragg) == 0) {
      cat("Region aggregation matrices: 0\n")
    } else {
      ragstr <- ifelse(length(x = ragg) > 6, "...", "")
      cat(
        length(x = ragg),
        "region aggregation",
        ifelse(length(x = ragg) == 1, "matrix:", "matrices:"),
        ragstr, "\n"
      )
    }
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
#' @method Annotation ChromatinAssay5
#' @export
#' @concept assay
#' @examples
#' \donttest{
#' Annotation(atac_small[["peaks"]])
#' }
Annotation.ChromatinAssay5 <- function(object, ...) {
  return(slot(object = object, name = "annotation"))
}

#' @param object A Seurat, GRangesAssay, or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @importFrom SeuratObject DefaultAssay
#' @rdname Annotation
#' @method Annotation Seurat
#' @export
#' @concept assay
#' @examples
#' \donttest{
#' Annotation(atac_small)
#' }
Annotation.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(Annotation(object = object[[assay]]))
}

#' @rdname Bias
#' @method Bias ChromatinAssay5
#' @export
#' @concept assay
#' @examples
#' \donttest{
#' Bias(atac_small[["peaks"]])
#' }
Bias.ChromatinAssay5 <- function(object, ...) {
  return(slot(object = object, name = "bias"))
}

#' @param object A Seurat, GRangesAssay, or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @rdname Bias
#' @method Bias Seurat
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept assay
#' @examples
#' \donttest{
#' Bias(atac_small)
#' }
Bias.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(Bias(object = object[[assay]]))
}

#' @rdname Fragments
#' @method Fragments ChromatinAssay5
#' @export
#' @concept assay
#' @concept fragments
#' @examples
#' Fragments(atac_small[["peaks"]])
Fragments.ChromatinAssay5 <- function(object, ...) {
  return(slot(object, name = "fragments"))
}

#' @param object A Seurat, GRangesAssay, or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @importFrom SeuratObject DefaultAssay
#' @rdname Fragments
#' @method Fragments Seurat
#' @export
#' @concept assay
#' @concept fragments
#' @examples
#' Fragments(atac_small)
Fragments.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(Fragments(object = object[[assay]]))
}

#' @rdname Motifs
#' @method Motifs ChromatinAssay5
#' @export
#' @concept assay
#' @concept motifs
#' @examples
#' Motifs(atac_small[["peaks"]])
Motifs.ChromatinAssay5 <- function(object, ...) {
  return(slot(object = object, name = "motifs"))
}

#' @param object A Seurat or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @rdname Motifs
#' @importFrom SeuratObject DefaultAssay
#' @method Motifs Seurat
#' @export
#' @concept assay
#' @concept motifs
#' @examples
#' Motifs(atac_small)
Motifs.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(Motifs(object = object[[assay]]))
}

#' @rdname Links
#' @method Links ChromatinAssay5
#' @export
#' @concept assay
#' @concept links
#' @examples
#' Links(atac_small[["peaks"]])
Links.ChromatinAssay5 <- function(object, ...) {
  return(slot(object = object, name = "links"))
}

#' @param object A Seurat or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @rdname Links
#' @method Links Seurat
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept links
#' @concept assay
#' @examples
#' Links(atac_small)
Links.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(Links(object = object[[assay]]))
}

#' @param features Optional character vector of region aggregation names to
#' return
#' @rdname RegionAggr
#' @method RegionAggr ChromatinAssay5
#' @export
#' @concept footprinting
#' @concept assay
#' @examples
#' RegionAggr(atac_small[["peaks"]])
RegionAggr.ChromatinAssay5 <- function(object, features = NULL, ...) {
  agg.list <- slot(object = object, name = "region.aggregation")

  if (!is.null(features)) {
    agg.list.names <- vapply(agg.list, function(x) x@name, character(1))
    agg.list <- agg.list[agg.list.names %in% features]
  }
  return(agg.list)
}


#' @param object A Seurat or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @rdname RegionAggr
#' @method RegionAggr Seurat
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept footprinting
#' @concept assay
#' @examples
#' RegionAggr(atac_small)
RegionAggr.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(RegionAggr(object = object[[assay]], ...))
}

#' List stored RegionAggregation objects
#'
#' Returns the names of all stored  [RegionAggregation-class] objects
#'
#' @param object A [ChromatinAssay5-class]  object
#' @rdname RegionAggNames
#' @method RegionAggNames ChromatinAssay5
#' @return Character vector of stored result names
#' @export
RegionAggNames.ChromatinAssay5 <- function(object, ...) {
  name.list <- vapply(RegionAggr(object), function(x) x@name, character(1))
  return(name.list)
}

#' @param object A Seurat or ChromatinAssay5 object
#' @param assay Name of assay to use
#' @rdname RegionAggNames
#' @method RegionAggNames Seurat
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept footprinting
#' @concept assay
#' @examples
#' RegionAggr(atac_small)
RegionAggNames.Seurat <- function(object, assay = NULL, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  return(RegionAggNames(object = object[[assay]]))
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
#' @method Motifs<- ChromatinAssay5
#' @concept assay
#' @concept motifs
#' @examples
#' motifs <- Motifs(atac_small)
#' Motifs(atac_small[["peaks"]]) <- motifs
"Motifs<-.ChromatinAssay5" <- function(object, ..., value) {
  object <- SetAssayData(object = object, layer = "motifs", new.data = value)
  return(object)
}

#' @export
#' @rdname Motifs
#' @method Motifs<- Seurat
#' @importFrom SeuratObject DefaultAssay
#' @concept assay
#' @concept motifs
#' @examples
#' motifs <- Motifs(atac_small)
#' Motifs(atac_small) <- motifs
"Motifs<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  Motifs(object = object[[assay]]) <- value
  return(object)
}

#' @export
#' @rdname Links
#' @method Links<- ChromatinAssay5
#' @concept assay
#' @concept links
#' @examples
#' links <- Links(atac_small)
#' Links(atac_small[["peaks"]]) <- links
"Links<-.ChromatinAssay5" <- function(object, key = NULL, ..., value) {
  object <- SetAssayData(
    object = object, layer = "links", new.data = value, key = key
  )
  return(object)
}

#' @param key Key to use when adding a new [InteractionSet::GInteractions]
#' object
#' @export
#' @rdname Links
#' @method Links<- Seurat
#' @concept assay
#' @concept links
#' @examples
#' links <- Links(atac_small)
#' Links(atac_small) <- links
"Links<-.Seurat" <- function(object, assay = NULL, key = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  Links(object[[assay]], key = key) <- value
  return(object)
}

#' @export
#' @rdname RegionAggr
#' @method RegionAggr<- ChromatinAssay5
#' @concept assay
#' @concept links
#' @examples
#' ra <- RegionAggr(atac_small)
#' RegionAggr(atac_small[["peaks"]]) <- ra
"RegionAggr<-.ChromatinAssay5" <- function(object, ..., value) {
  object <- SetAssayData(
    object = object, layer = "region.aggregation", new.data = value, ...
  )
  return(object)
}

#' @export
#' @rdname RegionAggr
#' @method RegionAggr<- Seurat
#' @concept assay
#' @concept links
#' @examples
#' ra <- RegionAggr(atac_small)
#' RegionAggr(atac_small) <- ra
"RegionAggr<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  RegionAggr(object[[assay]], ...) <- value
  return(object)
}

#' @export
#' @rdname Annotation
#' @concept assay
#' @method Annotation<- ChromatinAssay5
#' @examples
#' genes <- Annotation(atac_small)
#' Annotation(atac_small[["peaks"]]) <- genes
"Annotation<-.ChromatinAssay5" <- function(object, ..., value) {
  object <- SetAssayData(
    object = object, layer = "annotation", new.data = value
  )
  return(object)
}

#' @export
#' @importFrom SeuratObject DefaultAssay
#' @method Annotation<- Seurat
#' @concept assay
#' @rdname Annotation
#' @examples
#' genes <- Annotation(atac_small)
#' Annotation(atac_small) <- genes
"Annotation<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  Annotation(object = object[[assay]]) <- value
  return(object)
}

#' @export
#' @method Bias<- ChromatinAssay5
#' @concept assay
#' @rdname Bias
"Bias<-.ChromatinAssay5" <- function(object, ..., value) {
  object <- SetAssayData(object = object, layer = "bias", new.data = value)
  return(object)
}

#' @export
#' @importFrom SeuratObject DefaultAssay
#' @method Bias<- Seurat
#' @concept assay
#' @rdname Bias
#' @examples
#' bases <- c("A", "C", "G", "T")
#' hexamers <- apply(expand.grid(rep(list(bases), 6)), 1, paste0, collapse = "")
#' b <- 1:length(hexamers)
#' names(b) <- hexamers
#'
#' # assign bias in a Seurat object
#' Bias(atac_small) <- b
#'
#' # assign bias in a ChromatinAssay5 or GRangesAssay
#' Bias(atac_small[["peaks"]]) <- b
"Bias<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  Bias(object = object[[assay]]) <- value
  return(object)
}

#' @export
#' @method Fragments<- ChromatinAssay5
#' @rdname Fragments
#' @importFrom SeuratObject SetAssayData
#' @concept assay
#' @concept fragments
"Fragments<-.ChromatinAssay5" <- function(object, ..., value) {
  if (is.null(x = value)) {
    slot(object = object, name = "fragments") <- list()
    return(object)
  }
  # append new fragment objects to existing list
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
#' @method Fragments<- Seurat
#' @rdname Fragments
#' @concept assay
#' @concept fragments
#' @importFrom SeuratObject DefaultAssay
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' # Seurat object method
#' Fragments(atac_small) <- fragments
#'
#' # Remove fragment objects
#' Fragments(atac_small) <- NULL
#'
#' # Assay method
#' Fragments(atac_small[["peaks"]]) <- fragments
"Fragments<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  Fragments(object = object[[assay]]) <- value
  return(object)
}

# Add a single Fragment object to a ChromatinAssay5
# @param object A \code{\link{ChromatinAssay5}} object
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
    current.frags <- GetAssayData(object = object, layer = "fragments")
    for (i in seq_along(along.with = current.frags)) {
      if (any(Cells(x = fragments) %in% Cells(x = current.frags[[i]]))) {
        stop("Cells already present in a fragment object")
      }
    }
  }
  # check file is not already linked to the object
  current.frags <- GetAssayData(object = object, layer = "fragments")
  all.path <- lapply(
    X = current.frags,
    FUN = GetFragmentData,
    slot = "file.path"
  )
  new.path <- GetFragmentData(object = fragments, slot = "file.path")
  if (new.path %in% all.path) {
    stop("Fragment file already present in the object")
  }
  # append fragments to list
  current.frags[[length(x = current.frags) + 1]] <- fragments
  slot(object = object, name = "fragments") <- current.frags
  return(object)
}

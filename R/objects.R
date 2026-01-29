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
#' the rows of the input matrix
#' @param sep Separators to use for strings encoding genomic coordinates.
#' First element is used to separate the chromosome from the coordinates,
#' second element is used to separate the start from end coordinate. Only
#' used if `ranges` is NULL.
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
  sep = c("-", "-"),
  ...
) {
  chrom.assay <- CreateChromatinAssay5(counts = counts, data = data, ...)
  if (!missing(x = counts)) {
    features.keep <- rownames(x = counts) %in% rownames(x = chrom.assay)
  } else {
    features.keep <- rownames(x = data) %in% rownames(x = chrom.assay)
  }
  # subset ranges if there are features removed
  if (!is.null(x = ranges)) {
    ranges <- ranges[features.keep, ]
  }
  granges.assay <- as.GRangesAssay(
    x = chrom.assay,
    ranges = ranges,
    sep = sep,
    ...
  )
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
#' @param region.aggregation A named list of [RegionAggregation]
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
  if (!is.null(x = annotation) && !inherits(x = annotation, what = "GRanges")) {
    stop("Annotation must be a GRanges object.")
  }
  if (!is.null(x = annotation)) {
    if (!any(c("tx_id", "transcript_id") %in% colnames(x = mcols(x = annotation)))) {
      stop("Annotation must have transcript id stored in `tx_id` or `transcript_id`.")
    }
    if (any(!c("gene_name", "gene_id", "gene_biotype", "type") %in% colnames(x = mcols(x = annotation)))) {
      stop("Annotation must have `gene_name`, `gene_id`, `gene_biotype` and `type`.")
    }
  }

  # CreateAssay5Object throws error if counts or data is missing rather than NULL
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

  # TODO this should be move to the SetAssayData method
  if (!is.null(x = motifs)) {
    # pre-computed motif object, make sure features are formatted the same
    # as count matrix and subset features
    if (!inherits(x = motifs, what = "Motif")) {
      stop("Provided motif object is not a Motif-class object")
    }
    if (!is.null(x = GetMotifData(object = motifs, slot = "data"))) {
      if (!(all(rownames(x = motifs) == rownames(x = chrom.assay)))) {
        # rownames don't match
        motif.mat <- GetMotifData(object = motifs, slot = "data")
        # subset
        if (!all(rownames(x = chrom.assay) %in% rownames(x = motif.mat))) {
          warning(
            "Some peak regions missing from supplied motif object. ",
            "Motif information will not be added"
          )
          motifs <- NULL
        }
        motif.mat <- motif.mat[rownames(x = chrom.assay), ]
        motifs <- SetMotifData(
          object = motifs, slot = "data", new.data = motif.mat
        )
      }
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
  sep = c("-", "-"),
  ...
) {
  x <- as.ChromatinAssay5(
    object = x,
    annotation = annotation,
    fragments = fragments,
    bias = bias,
    motifs = motifs,
    links = links,
    region.aggregation = region.aggregation,
  )
  x <- as.GRangesAssay(
    object = x,
    ranges = ranges,
    sep = sep
  )
  return(x)
}

#' @param ranges A GRanges object
#' @param sep Characters used to separate the chromosome, start, and end
#' coordinates in the row names of the data matrix
#'
#' @rdname as.GRangesAssay
#' @export
#' @method as.GRangesAssay ChromatinAssay5
#' @concept assay
#'
as.GRangesAssay.ChromatinAssay5 <- function(
  x,
  ranges = NULL,
  sep = c("-", "-"),
  ...
) {
  if (!is.null(x = ranges)) {
    if (length(x = ranges) != nrow(x = x)) {
      stop("Length of supplied genomic ranges does not match number
           of rows in matrix")
    }
  } else {
    ranges <- StringToGRanges(regions = rownames(x = x), sep = sep)
  }
  if (!isDisjoint(x = ranges)) {
    warning("Overlapping ranges supplied. Ranges should be non-overlapping.")
  }

  # re-assign row names of matrix so that it's a known granges transformation
  new.rownames <- GRangesToString(grange = ranges, sep = c("-", "-"))
  rownames(x = x) <- new.rownames

  new.assay <- as(object = x, Class = "GRangesAssay")
  ranges <- ranges %||% StringToGRanges(regions = rownames(x = x), sep = sep)
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
#' @param fragments A list of [Fragment()] objects
#' @param bias Tn5 integration bias matrix
#' @param motifs A [Motif()] object
#' @param links Genomic links TODO
#' @param region.aggregation A named list of [RegionAggregation]
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

setAs(
  from = "ChromatinAssay5",
  to = "GRangesAssay",
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
        "Class" = "GRangesAssay"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)

setAs(
  from = "Assay5",
  to = "ChromatinAssay5",
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
        "Class" = "ChromatinAssay5"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)

## Functions
#' Create a RegionAggregation object 
#'
#' Create a [RegionAggregation-class()] object to store the insertions over group of regions 
#' 
#' @param matrix A cell-by-position matrix 
#' @param regions, A [GenomicRanges::granges()] object containing the
#' regions aggregated across. 
#' @param upstream Integer denoting number of bases upstream of the 
#' centered position that are stored in the matrix 
#' @param downstream Integer denoting number of bases downstream of the 
#' centered position that are stored in the matrix
#' @param name A name for the set of regions aggregated
#' @param expected A vector containing expected number of Tn5 insertions per position
#' @param cells A named vector of cells where each element is the cell barcode
#' and the name of each element is the corresponding cell barcode as stored 
#' in the ChromatinAssay5 object
#' 
#' @export 
#' @return Returns a [RegionAggregation()] object 
#' @concept regionaggregation 
CreateRegionAggregationObject <- function(
    matrix = NULL,
    regions = NULL,
    upstream = NULL, 
    downstream = NULL, 
    name = NULL, 
    expected = NULL, 
    cells = NULL
){
    if (!inherits(matrix, c("matrix", "CsparseMatrix", "dgCMatrix"))) {
        stop("data must be a matrix or sparse matrix. Supplied ", class(x=matrix))
    }
    if (!inherits(regions, "GRanges")){
        stop("regions must be a GRanges object")
    }
    # in footprinting if ocmpute.expected = False, 
    #  expected.insertions <- rep(1, width(x = dna.sequence)[[1]] - 6)
    if (!is.null(expected)){
        if (!is.numeric(expected)) {
            stop("expected insertions must be a numeric vector")
        }
    }
    # if cells not given, create the cells from rownames of the matrix
    if (is.null(cells)){
        cells <- setNames(rownames(matrix), rownames(matrix))
    } else {
        cells <- as.character(cells)
        if (length(cells) != length(rownames(matrix))){
            stop("Number of cells: (", length(cells), 
                 ") does not match number of matrix rows (",length(rownames(matrix)), ")")
        }
        # if unnamed or partially named -> force names from matrix 
        if (is.null(names(cells)) || any(names(cells) == "")) {
            names(cells) <- rownames(matrix)
        }
    }
        
    agg.obj <- new(
        Class = "RegionAggregation", 
        matrix = matrix, 
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
#' @param seqlevels A named vector of sequence levels (eg, chromosome name) where
#' each element is the sequence name as it appears in the fragment file, and the
#' name of each element is the corresponding sequence name as stored in the
#' object. If NULL, the sequence names are assumed to be the same in the
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
#' frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
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
#' @return Returns a [Motif()] object
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
        if (!all(colnames(new.object) %in% colnames(object[[expression.assay]]))) {
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
        if (!suppressWarnings(all(colnames(object[[expression.assay]]) == colnames(new.object)))) {
          warning("Subsetting expression assay to have same cells as chromatin assay.",
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
        warning("Cannot access layers if SeuratObject version is not 5.0.0 or greater.",
          "Please update SeuratObject to also visualize expression data when your",
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
  if (layer %in% c("counts", "data", "scale.data", "meta.data", "misc", "key")) {
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
  if (layer %in% c("counts", "data", "scale.data", "meta.data", "misc", "key")) {
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
#' Extract data from a [Fragment2-class()] object
#'
#' @param object A [Fragment2()] object
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
  # there's currently nothing cell-centric that needs to be renamed in the GRangesAssay class

  # rename cells in the parental class
  object <- NextMethod()

  # rename cells in fragment objects
  frags <- Fragments(object = object)
  for (i in seq_along(along.with = frags)) {
    frags[[i]] <- RenameCells(object = frags[[i]], new.names = new.names)
  }
  Fragments(object = object) <- NULL
  Fragments(object = object) <- frags

  region.aggr <- GetAssayData(object = object, layer = "region.aggregation")
  for (i in seq_along(along.with = region.aggr)) {
    # TODO implement RenameCells.RegionAggregation
    region.aggr[[i]] <- RenameCells(
      object = region.aggr[[i]],
      new.names = new.names
    )
  }
  # TODO implement region aggregation assignment method
  slot(object = object, name = "region.aggregation") <- region.aggr
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
  cells <- cells[names(x = new.names)]
  names(x = cells) <- new.names[names(x = cells)]
  slot(object = object, name = "cells") <- cells
  return(object)
}

#' @importFrom SeuratObject RenameCells
#' @concept RegionAggregation
#' @method RenameCells RegionAggregation
#' @export
RenameCells.RegionAggregation <- function(object, new.names, ...) {
  cells <- object@cells 
  # ^ update to cells <- GetRegionAggregation(object = object, slot = "cells")
  if (is.null(x = cells)) {
    stop(
      "Cannot rename cells in RegionAggregation object ",
      "with no cell information stored"
    )
  }
  # subset and rename 
  cells <- cells[names(x = new.names)]
  names(x = cells) <- new.names[names(x = cells)]
  slot(object = object, name = "cells") <- cells
  # also rename matrix rownames 
  old.rows <- rownames(object@matrix)
  # ^ update to get using GetRegionAggregation(object, slot = "matrix")
  rownames(object@matrix) <- new.names(old.rows, names(new.names))

  return(object)
}

#' @importFrom SeuratObject SetAssayData CheckFeaturesNames
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
    "fragments", "annotation", "bias", "region.aggregation"
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
  } else if (layer == "motifs") {
    if (is.null(x = new.data)) {
      methods::slot(object = object, name = layer) <- NULL
      return(object)
    }
    if (!inherits(x = new.data, what = "Motif")) {
      stop("Must provide a Motif class object")
    }
    # Set the feature names compatible with Seurat
    new.data <- SetMotifData(
      object = new.data,
      slot = "data",
      new.data = CheckFeaturesNames(
        data = GetMotifData(object = new.data, slot = "data")
      )
    )

    # TODO allow mismatching row names, but check that the genomic ranges
    # are equivalent. Requires adding a granges slot to the motif class
    if (nrow(x = object) != nrow(x = new.data) ||
      !all(rownames(x = object) == rownames(x = new.data))) {
      keep.features <- intersect(
        x = rownames(x = new.data),
        y = rownames(x = object)
      )
      if (length(x = keep.features) == 0) {
        stop("No features in common between the GRangesAssay
             and Motif objects")
      } else {
        warning("Features do not match in GRangesAssay and Motif object.
                Subsetting/Filling the Motif object.")
        new.data <- new.data[keep.features, ]

        new.data <- SetMotifData(
          object = new.data,
          slot = "data",
          new.data = AddMissing(
            GetMotifData(object = new.data, slot = "data"),
            features = rownames(x = object)
          )
        )
      }
    }
    methods::slot(object = object, name = layer) <- new.data
  } else if (layer == "links") {
    methods::slot(object = object, name = layer) <- new.data
  }
  return(object)
}

#' @importFrom SeuratObject SetAssayData CheckFeaturesNames
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
  if (layer %in% c("counts", "data", "scale.data", "meta.data", "misc", "key")) {
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
    if (inherits(x = new.data, what = "list")) {
      # check that it's a list containing fragment class objects
      for (i in seq_along(new.data)) {
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
    if (!any(c("tx_id", "transcript_id") %in% colnames(x = mcols(x = new.data)))) {
      stop("Annotation must have transcript id stored in `tx_id` or `transcript_id`.")
    }
    if (any(!c("gene_name", "gene_id", "gene_biotype", "type") %in% colnames(x = mcols(x = new.data)))) {
      stop("Annotation must have `gene_name`, `gene_id`, `gene_biotype` and `type`.")
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
    if (inherits(x = new.data, what = "list")){
      # check if its a list containing RegionAggregation class objects
      for (i in seq_along(new.data)){
        if (!inherits(x = new.data[[i]], what = "RegionAggregation")){
          stop("New data is not a RegionAggregation object")
        }
      }
    } else if (inherits(x = new.data, what = "RegionAggregation")){
      # single RegionAggregation object
      new.data <- list(new.data)
    }
    # get existing RegionAggregation object 
    agg.list <- GetAssayData(object = object, layer = layer)
    if (length(x = agg.list) == 0) {
      # nothing exists yet -> assign directly 
      methods::slot(object, "region.aggregation") <- new.data
      return(object)
    } 
    # try to merge compatible RegionAggregation objects 
    for (i in seq_along(new.data)){
      new.agg <- new.data[[i]]
      merged <- FALSE
      # compare against same-name objects 
      same.name.idx <- which(vapply(
        agg.list, function(x) identical(x@name, new.agg@name), logical(1)))
      if (length(same.name.idx) > 0) {
        for (j in same.name.idx){
          old.agg <- agg.list[[j]]
          compatible <- 
            identical(old.agg@upstream, new.agg@upstream) && 
            identical(old.agg@downstream, new.agg@downstream) && 
            identical(old.agg@expected, new.agg@expected) && 
            identical(old.agg@regions, new.agg@regions)
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
      if (!merged) { 
        agg.list <- c(agg.list, list(new.agg))
      }
    }
    methods::slot(object = object, layer) <- agg.list
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
    if (!all(rownames(x = object) == rownames(x = new.data))) {
      stop("Features do not match existing assay data.
           Column names in motif matrix should match row names in assay data")
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
#'   object = atac_small, assay = "peaks", slot = "data", new.data = motif.matrix
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
#' @return Returns a subsetted [Motif()] object
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
  new.data <- GetMotifData(object = x, slot = "data")[features, motifs, drop = FALSE]
  new.pwm <- GetMotifData(object = x, slot = "pwm")[motifs]
  new.names <- GetMotifData(object = x, slot = "motif.names")[motifs]
  new.meta <- GetMotifData(object = x, slot = "meta.data")[motifs, , drop = FALSE]
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

  # recompute meta features
  x <- FindTopFeatures(
    object = x,
    min.cutoff = NA,
    verbose = FALSE
  )
  # subset cells in region aggregation matrices
  # cells <- cells %||% colnames(x = x)
  posmat <- GetAssayData(object = x, layer = "region.aggregation")

  # TODO update for RegionAggregation class
  if (length(posmat) > 0) { 
    posmat <- lapply(posmat, subset, cells = cells)
    posmat <- Filter(Negate(is.null), posmat)
  }
  x <- SetAssayData(object = x, layer = "region.aggregation", new.data = posmat)

  # for (i in seq_along(along.with = posmat)) {
    # TODO need to make the formatting for positionEnrichment slot better defined
    # currently the RegionMatrix and Footprint functions write differently
    # formatted information
    # regionmatrix is group x position
    # footprint is cell x position
    # if (inherits(x = posmat[[i]], what = "list")) {
      # from RegionMatrix
      # group x position matrix
      # do not subset as we don't have per-cell information here
      # next
    # } else {
      # from Footprint
      # cell x position matrix
      # with expected and motif position rows
      # added_rows <- c("expected", "motif")
      # added_rows <- added_rows[added_rows %in% rownames(x = posmat[[i]])]
      # posmat[[i]] <- posmat[[i]][c(cells, added_rows), ]
    # }
  # }

  # subset cells in Fragments objects
  frags <- Fragments(object = x)
  Fragments(object = x) <- NULL
  for (i in seq_along(along.with = frags)) {
    frags[[i]] <- subset(x = frags[[i]], cells = cells)
  }
  Fragments(object = x) <- frags
  
  # subset motifs
  motifs <- Motifs(object = x)
  if (!is.null(x = motifs)) {
    motifs <- subset(x = motifs, features = features)
  }
  Motifs(object = x) <- motifs

  return(x)
}

#' Subset a single RegionAggregation object
#'
#' @method subset RegionAggregation 
#' @importFrom methods slot "slot<-"
#' @export
#' @concept RegionAggregation 
subset.RegionAggregation <- function(
    x,
    cells = NULL,
    ...
){
  if (is.null(x = cells)) {
    return(x)
  }
  mat <- slot(object = x, name = "matrix")
  agg.cells <- slot(object = x, name = "cells")
  
  keep <- names(x = agg.cells) %in% cells
  if (!any(keep)){
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
#' frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
#' subset(frags, head(names(cells)))
subset.Fragment2 <- function(
  x,
  cells = NULL,
  ...
) {
  frag.cells <- GetFragmentData(object = x, slot = "cells")
  keep <- fmatch(x = names(x = frag.cells), table = cells, nomatch = 0L) > 0
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
    granges.all <- StringToGRanges(regions = rownames(x = merged))
    merged <- as.GRangesAssay(x = merged, ranges = granges.all)
    return(merged)
  }
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
      x = assays[[1]], y = assays[[2:length(x = assays)]], ...
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
          new.names = paste(add.cell.ids[i], colnames(x = assays[[i]]), sep = "_")
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
      GetAssayData(object = x, layer = "region.aggregation")
    })
    # flatten list-of-lists 
    all.agg <- Reduce(f=c, x= all.agg)
    

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
    cat(
      "Region aggregation matrices:",
      length(x = GetAssayData(
        object = object,
        layer = "region.aggregation"
      )),
      "\n"
    )
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
    cat(
      "Region aggregation matrices:",
      length(x = GetAssayData(
        object = object,
        layer = "region.aggregation"
      )),
      "\n"
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

#' @rdname RegionAggr
#' @method RegionAggr ChromatinAssay5
#' @export
#' @concept footprinting
#' @concept assay
#' @examples
#' RegionAggr(atac_small[["peaks"]])
RegionAggr.ChromatinAssay5 <- function(object, ...) {
  return(slot(object = object, name = "region.aggregation"))
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
  return(RegionAggr(object = object[[assay]]))
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
"Links<-.ChromatinAssay5" <- function(object, ..., value) {
  object <- SetAssayData(object = object, layer = "links", new.data = value)
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
"Links<-.Seurat" <- function(object, assay = NULL, ..., value) {
  assay <- assay %||% DefaultAssay(object = object)
  Links(object[[assay]]) <- value
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
  object <- SetAssayData(object = object, layer = "region.aggregation", new.data = value)
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
  RegionAggr(object[[assay]]) <- value
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
  object <- SetAssayData(object = object, layer = "annotation", new.data = value)
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
#' bases <- c("A","C","G","T")
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
  if (inherits(x = value, what = "list")) {
    for (i in seq_along(along.with = value)) {
      object <- AddFragments(object = object, fragments = value[[i]])
    }
  } else if (is.null(x = value)) {
    object <- SetAssayData(
      object = object,
      layer = "fragments",
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
  # append fragments to list
  current.frags <- GetAssayData(object = object, layer = "fragments")
  current.frags[[length(x = current.frags) + 1]] <- fragments
  slot(object = object, name = "fragments") <- current.frags
  return(object)
}

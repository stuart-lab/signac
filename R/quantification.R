#' @include generics.R
NULL

#' @param genome genome A vector of chromosome sizes for the genome. This is
#' used to construct the genome bin coordinates. The can be obtained by calling
#' seqlengths on a BSgenome-class object.
#' @param assay Name of assay to use
#' @param new.assay.name Name of new assay to create containing aggregated
#' genome tiles
#' @param min_counts Minimum number of counts for a tile to be retained prior to
#' aggregation
#' @param binsize Size of the genome bins (tiles) in base pairs
#' @param verbose Display messages
#'
#' @rdname AggregateTiles
#' @importFrom Seurat DefaultAssay
#' @export
#' @method AggregateTiles Seurat
#' @concept quantification
#' @return When running on a Seurat object, returns the Seurat object with a new
#' \code{\link{ChromatinAssay}} added.
AggregateTiles.Seurat <- function(
  object,
  genome,
  assay = NULL,
  new.assay.name = "tiles",
  min_counts = 5,
  binsize = 5000,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  object[[new.assay.name]] <- AggregateTiles(
    object = object[[assay]],
    genome = genome,
    min_counts = min_counts,
    binsize = binsize,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @rdname AggregateTiles
#' @export
#' @method AggregateTiles ChromatinAssay
#' @concept quantification
#' @return When running on a \code{\link{ChromatinAssay}}, returns a new
#' \code{ChromatinAssay} containing the aggregated genome tiles.
AggregateTiles.ChromatinAssay <- function(
  object,
  genome,
  min_counts = 5,
  binsize = 5000,
  verbose = TRUE,
  ...
) {
  frags <- Fragments(object = object)
  bins <- AggregateTiles(
    object = frags,
    genome = genome,
    cells = colnames(x = object),
    min_counts = min_counts,
    binsize = binsize,
    verbose = verbose,
    ...
  )
  if (verbose) {
    message("Constructing assay")
  }
  assay.obj <- CreateChromatinAssay(
    counts = bins,
    fragments = frags
  )
  return(assay.obj)
}

#' @param cells Cells to include
#' @rdname AggregateTiles
#' @importFrom Matrix rowSums
#' @export
#' @method AggregateTiles default
#' @concept quantification
#' @return When running on a fragment file, returns a sparse region x cell
#' matrix.
AggregateTiles.default <- function(
  object,
  genome,
  cells = NULL,
  min_counts = 5,
  binsize = 5000,
  verbose = TRUE,
  ...
) {
  # quantify genome bins
  bins <- GenomeBinMatrix(
    fragments = object,
    genome = genome,
    cells = cells,
    binsize = binsize,
    verbose = verbose
  )

  # filter out low coverage bins
  keep.rows <- rowSums(x = bins) > min_counts
  if (sum(x = keep.rows) == 0) {
    stop("No bins found with over ", min_counts, " cells")
  }
  bins <- bins[keep.rows, ]

  # join adjacent bins
  if (verbose) {
    message("Combining adjacent tiles")
  }
  aggregate.tiles <- CombineTiles(bins = bins)
  return(aggregate.tiles)
}

#' Genome bin matrix
#'
#' Construct a bin x cell matrix from a fragments file.
#'
#' This function bins the genome and calls \code{\link{FeatureMatrix}} to
#' construct a bin x cell matrix.
#'
#' @param fragments Path to tabix-indexed fragments file or a list of
#' \code{\link{Fragment}} objects
#' @param genome A vector of chromosome sizes for the genome. This is used to
#' construct the genome bin coordinates. The can be obtained by calling
#' \code{\link[GenomeInfoDb]{seqlengths}} on a
#' \code{\link[BSgenome]{BSgenome-class}} object.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param binsize Size of the genome bins to use
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param verbose Display messages
#'
#' @importFrom GenomicRanges tileGenome
#' @export
#' @concept quantification
#' @return Returns a sparse matrix
#' @examples
#' \donttest{
#' genome <- 780007
#' names(genome) <- 'chr1'
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(fpath)
#' GenomeBinMatrix(
#'   fragments = fragments,
#'   genome = genome,
#'   binsize = 1000
#' )
#' }
GenomeBinMatrix <- function(
  fragments,
  genome,
  cells = NULL,
  binsize = 5000,
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
) {
  tiles <- tileGenome(
    seqlengths = genome,
    tilewidth = binsize,
    cut.last.tile.in.chrom = TRUE
  )
  binmat <- FeatureMatrix(
    fragments = fragments,
    features = tiles,
    cells = cells,
    process_n = process_n,
    sep = sep,
    verbose = verbose
  )
  return(binmat)
}

#' Feature Matrix
#'
#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments A list of \code{\link{Fragment}} objects.
#' @param features A GRanges object containing a set of genomic intervals.
#' These will form the rows of the matrix, with each entry recording the number
#' of unique reads falling in the genomic region for each cell.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param verbose Display messages
#'
#' @export
#' @importFrom Seurat RowMergeSparseMatrices
#' @concept quantification
#' @return Returns a sparse matrix
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(fpath)
#' FeatureMatrix(
#'   fragments = fragments,
#'   features = granges(atac_small)
#' )
FeatureMatrix <- function(
  fragments,
  features,
  cells = NULL,
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
) {
  if (!inherits(x = fragments, what = "list")) {
    if (inherits(x = fragments, what = "Fragment")) {
      fragments <- list(fragments)
    } else {
      stop("fragments should be a list of Fragment objects")
    }
  }
  # if cells is not NULL, iterate over all fragment objects
  # and find which objects contain cells that are requested
  if (!is.null(x = cells)) {
    obj.use <- c()
    for (i in seq_along(along.with = fragments)) {
      if (any(cells %in% Cells(x = fragments[[i]]))) {
        obj.use <- c(obj.use, i)
      }
    }
  } else {
    obj.use <- seq_along(along.with = fragments)
  }
  # create a matrix from each fragment file
  mat.list <- sapply(
    X = obj.use,
    FUN = function(x) {
      SingleFeatureMatrix(
        fragment = fragments[[x]],
        features = features,
        cells = cells,
        sep = sep,
        verbose = verbose,
        process_n = process_n
      )
    })
  # merge all the matrices
  if (length(x = mat.list) == 1) {
    return(mat.list[[1]])
  } else {
    featmat <- Reduce(f = RowMergeSparseMatrices, x = mat.list)
    return(featmat)
  }
}

#### Not Exported ####

# matrix multiplication method for summing matrix rows
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors elementNROWS
#' @importFrom Matrix crossprod sparseMatrix
#' @importMethodsFrom Matrix t
CombineTiles <- function(bins) {
  ranges <- StringToGRanges(regions = rownames(x = bins))
  reduced.tiles <- reduce(x = ranges, with.revmap = TRUE)
  rmap <- reduced.tiles$revmap

  # construct matrix
  collapse_matrix <- sparseMatrix(
    i = unlist(x = rmap),
    j = rep(x = seq_along(rmap), times = elementNROWS(x = rmap)),
    x = 1
  )

  # sum bin matrix rows via matrix multiplication
  collapsed <- crossprod(x = bins, y = collapse_matrix)
  collapsed <- t(x = collapsed)
  rownames(x = collapsed) <- GRangesToString(grange = reduced.tiles)

  return(collapsed)
}

# Run FeatureMatrix on a single Fragment object
# @inheritParams FeatureMatrix
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom Matrix sparseMatrix
#' @importMethodsFrom GenomicRanges intersect
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom fastmatch fmatch
SingleFeatureMatrix <- function(
  fragment,
  features,
  cells = NULL,
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
) {
  fragment.path <- GetFragmentData(object = fragment, slot = "path")
  if (!is.null(cells)) {
    # only look for cells that are in the fragment file
    frag.cells <- GetFragmentData(object = fragment, slot = "cells")
    # first subset frag.cells
    cell.idx <- fmatch(
      x = names(x = frag.cells),
      table = cells,
      nomatch = 0L
    ) > 0
    cells <- frag.cells[cell.idx]
  }
  tbx <- TabixFile(file = fragment.path)
  features <- keepSeqlevels(
    x = features,
    value = intersect(
      x = seqnames(x = features),
      y = seqnamesTabix(file = tbx)
    ),
    pruning.mode = "coarse"
  )
  if (length(x = features) == 0) {
    stop("No matching chromosomes found in fragment file.")
  }

  feature.list <- ChunkGRanges(
    granges = features,
    nchunk = ceiling(x = length(x = features) / process_n)
  )
  if (verbose) {
    message("Extracting reads overlapping genomic regions")
  }
  if (nbrOfWorkers() > 1) {
    matrix.parts <- future_lapply(
      X = feature.list,
      FUN = PartialMatrix,
      tabix = tbx,
      cells = cells,
      sep = sep,
      future.globals = list(),
      future.scheduling = FALSE
    )
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    matrix.parts <- mylapply(
      X = feature.list,
      FUN = PartialMatrix,
      tabix = tbx,
      cells = cells,
      sep = sep
    )
  }
  # remove any that are NULL (no fragments for any cells in the region)
  null.parts <- sapply(X = matrix.parts, FUN = is.null)
  matrix.parts <- matrix.parts[!null.parts]
  if (is.null(x = cells)) {
    all.cells <- unique(
      x = unlist(x = lapply(X = matrix.parts, FUN = colnames))
    )
    matrix.parts <- lapply(
      X = matrix.parts,
      FUN = AddMissingCells,
      cells = all.cells
    )
  }
  featmat <- do.call(what = rbind, args = matrix.parts)
  if (!is.null(x = cells)) {
    # cells supplied, rename with cell name from object rather than file
    cell.convert <- names(x = cells)
    names(x = cell.convert) <- cells
    colnames(x = featmat) <- unname(obj = cell.convert[colnames(x = featmat)])
  }
  # reorder features
  feat.str <- GRangesToString(grange = features, sep = sep)
  featmat <- featmat[feat.str, ]
  return(featmat)
}

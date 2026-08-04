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
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @method AggregateTiles Seurat
#' @concept quantification
#' @return When running on a Seurat object, returns the Seurat object with a new
#' [ChromatinAssay5-class] assay added.
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
  assay <- assay %||% DefaultAssay(object = object)
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
#' @method AggregateTiles ChromatinAssay5
#' @concept quantification
#' @return When running on a [ChromatinAssay5-class], returns a new
#' `ChromatinAssay5` containing the aggregated genome tiles.
AggregateTiles.ChromatinAssay5 <- function(
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
  assay.obj <- CreateGRangesAssay(counts = bins, fragments = frags)
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
    verbose = verbose,
    ...
  )

  # filter out low coverage bins
  keep.rows <- rowSums(x = bins) >= min_counts
  if (sum(x = keep.rows) == 0) {
    stop("No bins found with at least ", min_counts, " counts")
  }
  bins <- bins[keep.rows, ]

  # join adjacent bins
  if (verbose) {
    message("Combining adjacent tiles")
  }
  aggregate.tiles <- CombineTiles(bins = bins)
  return(aggregate.tiles)
}

#' Create gene activity matrix
#'
#' Compute counts per cell in gene body and promoter region.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. If NULL, use the default assay
#' @param features Genes to include. If NULL, use all protein-coding genes in
#' the annotations stored in the object
#' @param extend.upstream Number of bases to extend upstream of the TSS
#' @param extend.downstream Number of bases to extend downstream of the TTS
#' @param biotypes Gene biotypes to include. If NULL, use all biotypes in the
#' gene annotation.
#' @param max.width Maximum allowed gene width for a gene to be quantified.
#' Setting this parameter can avoid quantifying extremely long transcripts that
#' can add a relatively long amount of time. If NULL, do not filter genes based
#' on width.
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param fragtk Use `fragtk` for fast and memory-efficient data
#' quantification. Can be TRUE/FALSE or a character vector. If TRUE,
#' `fragtk` will be used and attempt to find the `fragtk` executable
#' in the path. If FALSE, use the R implementation to produce the data matrix.
#' If a character vector is provided, this should be the path to the
#' `fragtk` executable and `fragtk` will be used. See
#' <https://crates.io/crates/fragtk> for fragtk documentation.
#' @param pic Use Paired Insertion Counting. If TRUE (default), each fragment
#' contributes at most 1 count per gene region. If FALSE, each insertion site
#' is counted separately (fragments with both ends in a region contribute 2
#' counts).
#' @param gene.id Record gene IDs in output matrix rather than gene name.
#' @param bpcells If `TRUE`, return a `BPCells::IterableMatrix` backed by an
#' on-disk BPCells directory instead of an in-memory sparse matrix. See
#' [FeatureMatrix()].
#' @param bpcells.dir Optional path to persist the BPCells directory beyond
#' the R session. See [FeatureMatrix()].
#' @param verbose Display messages
#'
#' @return Returns a sparse matrix, or a `BPCells::IterableMatrix` when
#' `bpcells = TRUE`.
#'
#' @concept utilities
#' @export
#' @importFrom SeuratObject DefaultAssay
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' GeneActivity(atac_small, fragtk = FALSE)
GeneActivity <- function(
  object,
  assay = NULL,
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = "protein_coding",
  max.width = 500000,
  process_n = 2000,
  fragtk = TRUE,
  pic = TRUE,
  gene.id = FALSE,
  bpcells = FALSE,
  bpcells.dir = NULL,
  verbose = TRUE
) {
  if (!is.null(x = features)) {
    if (length(x = features) == 0) {
      stop("Empty list of features provided")
    }
  }
  ValidateBPCellsArgs(bpcells = bpcells, bpcells.dir = bpcells.dir)
  # collapse to longest protein coding transcript
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("The requested assay is not a ChromatinAssay5.")
  }
  annotation <- Annotation(object = object[[assay]])
  if (is.null(x = annotation) || length(x = annotation) == 0) {
    stop("No gene annotations present in object")
  }
  # replace NA names with gene ID
  annotation$gene_name <- ifelse(
    test = is.na(x = annotation$gene_name) | (annotation$gene_name == ""),
    yes = annotation$gene_id,
    no = annotation$gene_name
  )
  if (verbose) {
    message("Extracting gene coordinates")
  }
  transcripts <- GetGeneRanges(ranges = annotation)
  if (gene.id) {
    transcripts$gene_name <- transcripts$gene_id
  }
  if (!is.null(x = biotypes)) {
    transcripts <- transcripts[transcripts$gene_biotype %in% biotypes]
    if (length(x = transcripts) == 0) {
      stop("No genes remaining after filtering for requested biotypes")
    }
  }

  # filter genes if provided
  if (!is.null(x = features)) {
    transcripts <- transcripts[transcripts$gene_name %in% features]
    if (length(x = transcripts) == 0) {
      stop("None of the requested genes were found in the gene annotation")
    }
  }
  if (!is.null(x = max.width)) {
    transcript.keep <- which(x = width(x = transcripts) < max.width)
    transcripts <- transcripts[transcript.keep]
    if (length(x = transcripts) == 0) {
      stop("No genes remaining after filtering for max.width")
    }
  }

  # extend to include promoters
  transcripts <- Extend(
    x = transcripts,
    upstream = extend.upstream,
    downstream = extend.downstream
  )

  # quantify. When bpcells is requested, route FeatureMatrix output through an
  # intermediate session-lifetime directory so that gene-name row renaming and
  # empty-row dropping can happen before the final matrix is persisted to the
  # user-supplied `bpcells.dir`.
  intermediate.dir <- if (isTRUE(x = bpcells)) {
    tempfile(pattern = "signac_geneactivity_")
  } else {
    NULL
  }
  counts <- FeatureMatrix(
    object = object[[assay]],
    features = transcripts,
    process_n = process_n,
    fragtk = fragtk,
    pic = pic,
    bpcells = bpcells,
    bpcells.dir = intermediate.dir,
    verbose = verbose
  )
  # set row names
  gene.key <- transcripts$gene_name
  names(x = gene.key) <- as.character(x = transcripts)
  rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
  counts <- counts[rownames(x = counts) != "", ]

  if (isTRUE(x = bpcells)) {
    counts <- AsBPCells(
      mat = counts, bpcells = TRUE, bpcells.dir = bpcells.dir
    )
    unlink(x = intermediate.dir, recursive = TRUE)
  }
  return(counts)
}

#' Genome bin matrix
#'
#' Construct a bin x cell matrix from a fragments file.
#'
#' This function bins the genome and calls [FeatureMatrix()] to
#' construct a bin x cell matrix.
#'
#' @param fragments Path to tabix-indexed fragments file or a list of
#' [Fragment()] objects
#' @param genome A vector of chromosome sizes for the genome. This is used to
#' construct the genome bin coordinates. The can be obtained by calling
#' `seqlengths` on a
#' [BSgenome::BSgenome-class()] object.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param binsize Size of the genome bins to use
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param bpcells If `TRUE`, return a `BPCells::IterableMatrix` backed by an
#' on-disk BPCells directory instead of an in-memory sparse matrix. See
#' [FeatureMatrix()].
#' @param bpcells.dir Optional path to persist the BPCells directory beyond
#' the R session. See [FeatureMatrix()].
#' @param verbose Display messages.
#' @param ... Arguments passed to [FeatureMatrix()].
#'
#' @importFrom GenomicRanges tileGenome
#' @export
#' @concept quantification
#' @return Returns a sparse matrix, or a `BPCells::IterableMatrix` when
#' `bpcells = TRUE`.
#' @examples
#' \donttest{
#' genome <- 780007
#' names(genome) <- "chr1"
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' fragments <- CreateFragmentObject(fpath, cells = colnames(atac_small))
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
  bpcells = FALSE,
  bpcells.dir = NULL,
  verbose = TRUE,
  ...
) {
  tiles <- tileGenome(
    seqlengths = genome,
    tilewidth = binsize,
    cut.last.tile.in.chrom = TRUE
  )
  if (inherits(x = fragments, what = "list")) {
    binmat <- FeatureMatrixList(
      frags = fragments,
      features = tiles,
      cells = cells,
      process_n = process_n,
      bpcells = bpcells,
      bpcells.dir = bpcells.dir,
      verbose = verbose,
      ...
    )
  } else {
    binmat <- FeatureMatrix(
      object = fragments,
      features = tiles,
      cells = cells,
      process_n = process_n,
      bpcells = bpcells,
      bpcells.dir = bpcells.dir,
      verbose = verbose,
      ...
    )
  }
  return(binmat)
}

#' @param features A [GenomicRanges::GRanges()] object containing the genomic
#' intervals to quantify. These form the rows of the returned matrix. A
#' character vector of UCSC-style coordinates (e.g. `"chr1:1-1000"`) is also
#' accepted and will be coerced to a `GRanges`.
#' @param assay Name of assay to use. If NULL, the default assay is used. The
#' assay must be a [ChromatinAssay5-class] with fragment information attached
#' via [Fragments()].
#' @param cells Character vector of cell barcodes to include as columns in the
#' output. If NULL, all cells present in the fragment file(s) are included;
#' this requires the R backend, since `fragtk` needs the set of cells up front
#' (see `fragtk`). When called on a [SeuratObject::Seurat] or
#' [ChromatinAssay5-class] object, the assay's `colnames` are used as the
#' default. Cell names are matched against the object-level names recorded in
#' the [Fragment2-class] `cells` slot; file-level barcodes are resolved
#' internally.
#' @param fragtk Use the `fragtk` backend for quantification. `TRUE` (default)
#' looks up `fragtk` on `PATH`; `FALSE` uses the R implementation; a character
#' string is interpreted as an explicit path to the `fragtk` executable. The R
#' backend is typically faster for small numbers of features (e.g. a handful
#' of peaks); `fragtk` is faster for genome-scale quantifications. See
#' <https://crates.io/crates/fragtk>.
#' @param pic Use Paired Insertion Counting. If `TRUE` (default), each
#' fragment contributes at most 1 count per feature, regardless of whether
#' one or both insertion sites fall inside the feature. If `FALSE`, each
#' insertion site (fragment start and end) is counted independently, so a
#' fragment with both ends inside a feature contributes 2 counts.
#' @param group Group features and sum counts across grouped rows. If `FALSE`
#' (default), no grouping is performed and each feature becomes a row in the
#' output. If `TRUE`, group by the first metadata column of `features` (from
#' `mcols(features)`). If a character string, group by the named metadata
#' column of `features`. Useful for collapsing exon- or transcript-level
#' intervals to gene-level counts.
#' @param keep_all_features By default, features on chromosomes not present in
#' the fragment file are dropped with a warning. Set `keep_all_features =
#' TRUE` to keep every feature in `features`; rows for absent chromosomes are
#' filled with zeros. Only honored by the R backend (`fragtk = FALSE`).
#' @param file.index Path to the tabix index (`.tbi` or `.csi`) for the
#' fragment file. If `NULL`, the index is located automatically from the
#' fragment file path (preferring `.tbi` when both are present). Only used by
#' the R backend (`fragtk = FALSE`).
#' @param frag.cells A named character vector mapping object-level cell names
#' (the `names()`) to file-level barcodes (the values). Used to translate
#' between the cell names used in a Seurat object and the barcodes written in
#' the fragment file. Typically `NULL` when calling with a raw fragment file
#' path; populated automatically when dispatched from a [Fragment2-class]
#' object via the `cells` slot.
#' @param seqlevels A named character vector specifying a seqname conversion
#' used to rename `features` to match the fragment file before quantification.
#' The names of the vector are the seqnames as they appear in `features` and
#' the values are the corresponding seqnames in the fragment file (e.g.
#' `c(chr1 = "1", chr2 = "2")` to convert UCSC-style names to Ensembl). Row
#' names of the returned matrix are mapped back to the original `features`
#' seqnames after quantification, so the output always matches the input
#' coordinate system regardless of what is stored on disk. Typically `NULL`
#' for raw path input; set automatically from the [Fragment2-class]
#' `seqlevels` slot.
#' @param process_n Number of regions to load into memory at a time, per
#' worker. Larger values can be faster but use more memory. Only affects the R
#' backend.
#' @param bpcells Logical. If `TRUE`, write the count matrix to disk in
#' BPCells format at `bpcells.dir` and return a `BPCells::IterableMatrix`
#' instead of an in-memory sparse matrix. Requires the `BPCells` package.
#' Default `FALSE`.
#' @param bpcells.dir Character. Path to a directory where the BPCells output
#' will be written. Required when `bpcells = TRUE`. The directory must not
#' already exist as a non-empty directory. Ignored when `bpcells = FALSE`.
#' @param verbose Display progress messages.
#'
#' @rdname FeatureMatrix
#' @export
#' @method FeatureMatrix Seurat
#' @concept quantification
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' FeatureMatrix(atac_small, features = granges(atac_small), fragtk = FALSE)
FeatureMatrix.Seurat <- function(
  object,
  features,
  assay = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  FeatureMatrix(object = object[[assay]], features = features, ...)
}

#' @rdname FeatureMatrix
#' @export
#' @method FeatureMatrix ChromatinAssay5
#' @concept quantification
FeatureMatrix.ChromatinAssay5 <- function(
  object,
  features,
  cells = NULL,
  ...
) {
  frags <- Fragments(object = object)
  if (length(x = frags) == 0) {
    stop("No fragment information found for requested assay")
  }
  cells <- cells %||% colnames(x = object)
  FeatureMatrixList(
    frags = frags,
    features = features,
    cells = cells,
    ...
  )
}

#' @rdname FeatureMatrix
#' @export
#' @method FeatureMatrix Fragment2
#' @concept quantification
FeatureMatrix.Fragment2 <- function(
  object,
  features,
  cells = NULL,
  ...
) {
  FeatureMatrix(
    object = GetFragmentData(object = object, slot = "file.path"),
    features = features,
    cells = cells,
    file.index = GetFragmentData(object = object, slot = "file.index"),
    frag.cells = GetFragmentData(object = object, slot = "cells"),
    seqlevels = GetFragmentData(object = object, slot = "seqlevels"),
    ...
  )
}

#' @rdname FeatureMatrix
#' @export
#' @method FeatureMatrix default
#' @importFrom SeuratObject RowMergeSparseMatrices
#' @importFrom GenomeInfoDb renameSeqlevels
#' @importFrom fastmatch fmatch
#' @concept quantification
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' FeatureMatrix(
#'   object = fpath,
#'   features = granges(atac_small),
#'   fragtk = FALSE
#' )
FeatureMatrix.default <- function(
  object,
  features,
  cells = NULL,
  file.index = NULL,
  frag.cells = NULL,
  seqlevels = NULL,
  fragtk = TRUE,
  pic = TRUE,
  group = FALSE,
  keep_all_features = FALSE,
  process_n = 2000,
  bpcells = FALSE,
  bpcells.dir = NULL,
  verbose = TRUE,
  ...
) {
  if (!is.character(x = object) || length(x = object) != 1) {
    stop("object should be a single fragment file path")
  }
  if (!inherits(x = features, what = "GRanges")) {
    if (inherits(x = features, what = "character")) {
      features <- GRanges(features)
    } else {
      stop("features should be a GRanges object")
    }
  }
  ValidateBPCellsArgs(bpcells = bpcells, bpcells.dir = bpcells.dir)

  # resolve file-level cells and object<->file barcode mapping
  if (!is.null(x = frag.cells)) {
    if (!is.null(x = cells)) {
      keep <- fmatch(
        x = names(x = frag.cells), table = cells, nomatch = 0L
      ) > 0
      cells.use <- frag.cells[keep]
    } else {
      cells.use <- frag.cells
    }
  } else if (!is.null(x = cells)) {
    cells.use <- cells
    names(x = cells.use) <- cells
  } else {
    cells.use <- NULL
  }

  # rename features to match the fragment file's seqlevels (pre)
  if (!is.null(x = seqlevels)) {
    feat.use <- suppressWarnings(
      expr = renameSeqlevels(x = features, value = seqlevels)
    )
  } else {
    feat.use <- features
  }

  fragtk.path <- NULL
  if (is.character(x = fragtk)) {
    fragtk.path <- fragtk
    fragtk <- TRUE
  }
  grouped <- (is.logical(x = group) && group) || is.character(x = group)

  if (fragtk && is.null(x = cells.use)) {
    # `fragtk matrix` requires a --cells file, so the set of cells has to be
    # known up front
    stop(
      "`fragtk = TRUE` requires a set of cells. Supply `cells`, or use ",
      "`fragtk = FALSE` to quantify every cell found in the fragment file."
    )
  }
  # fragtk writes its BPCells output to an intermediate directory, since the
  # row and column names are only finalized below. It is removed once the
  # matrix has been re-persisted to the user-supplied bpcells.dir.
  intermediate.dir <- NULL
  if (fragtk) {
    intermediate.dir <- if (isTRUE(x = bpcells)) {
      tempfile(pattern = "signac_fragtk_")
    } else {
      NULL
    }
    mat <- RunFragtk(
      fragments = object,
      features = feat.use,
      cells = unname(obj = cells.use),
      group = group,
      pic = pic,
      fragtk.path = fragtk.path,
      seqlevels = NULL,
      bpcells = bpcells,
      bpcells.dir = intermediate.dir,
      verbose = verbose,
      cleanup = TRUE
    )
    if (!grouped) {
      rownames(x = mat) <- as.character(x = feat.use)
    }
  } else {
    file.index <- file.index %||% GetIndexFile(
      fragment = object, verbose = verbose
    )
    mat <- SingleFeatureMatrix(
      path = object,
      file.index = file.index,
      features = feat.use,
      cells = cells.use,
      pic = pic,
      group = group,
      keep_all_features = keep_all_features,
      process_n = process_n,
      verbose = verbose
    )
  }

  # remap column names back to object-level cell names
  if (!is.null(x = cells.use) && !is.null(x = names(x = cells.use))) {
    to.obj <- names(x = cells.use)
    names(x = to.obj) <- unname(obj = cells.use)
    colnames(x = mat) <- unname(obj = to.obj[colnames(x = mat)])
  }

  # map row names back to original seqlevels (non-grouped only)
  if (!is.null(x = seqlevels) && !grouped) {
    sl <- names(x = seqlevels)
    names(x = sl) <- seqlevels
    orig <- suppressWarnings(
      expr = renameSeqlevels(x = feat.use, value = sl)
    )
    rownames(x = mat) <- as.character(x = orig)
  }

  # persist final BPCells matrix to disk with final dimnames
  mat <- AsBPCells(
    mat = mat, bpcells = bpcells, bpcells.dir = bpcells.dir
  )
  if (!is.null(x = intermediate.dir)) {
    unlink(x = intermediate.dir, recursive = TRUE)
  }

  return(mat)
}

# Iterate FeatureMatrix over a list of Fragment2 objects and merge
# @param frags A list of Fragment2 objects
# @param features A GRanges object
# @param cells Object-level cells to include (or NULL)
# @param ... Additional arguments passed to FeatureMatrix.Fragment2
#' @importFrom SeuratObject Cells
FeatureMatrixList <- function(
  frags,
  features,
  cells = NULL,
  group = FALSE,
  keep_all_features = FALSE,
  bpcells = FALSE,
  bpcells.dir = NULL,
  ...
) {
  if (!inherits(x = features, what = "GRanges")) {
    if (inherits(x = features, what = "character")) {
      features <- GRanges(features)
    } else {
      stop("features should be a GRanges object")
    }
  }
  ValidateBPCellsArgs(bpcells = bpcells, bpcells.dir = bpcells.dir)
  # filter frags to those containing any of the requested cells
  if (!is.null(x = cells)) {
    obj.use <- c()
    for (i in seq_along(along.with = frags)) {
      if (is.null(x = Cells(x = frags[[i]]))) {
        obj.use <- c(obj.use, i)
      } else if (any(cells %in% Cells(x = frags[[i]]))) {
        obj.use <- c(obj.use, i)
      }
    }
  } else {
    obj.use <- seq_along(along.with = frags)
  }
  if (length(x = obj.use) == 0) {
    stop(
      "None of the requested cells were found in any of the fragment objects ",
      "for this assay"
    )
  }
  # quantify each fragment in-memory; we BPCells-persist only at the end so
  # the Reduce(`+`) merge below operates on sparse matrices.
  mat.list <- lapply(
    X = frags[obj.use],
    FUN = FeatureMatrix,
    features = features,
    cells = cells,
    group = group,
    keep_all_features = keep_all_features,
    bpcells = FALSE,
    bpcells.dir = NULL,
    ...
  )
  if (length(x = mat.list) == 1) {
    return(AsBPCells(
      mat = mat.list[[1]], bpcells = bpcells, bpcells.dir = bpcells.dir
    ))
  }
  grouped <- (is.logical(x = group) && group) || is.character(x = group)
  all.cells <- unique(
    x = unlist(x = lapply(X = mat.list, FUN = colnames))
  )
  if (grouped) {
    feat.names <- unique(
      x = unlist(x = lapply(X = mat.list, FUN = rownames))
    )
  } else {
    feat.names <- as.character(x = features)
  }
  mat.list <- lapply(
    X = mat.list,
    FUN = AddMissing,
    cells = all.cells,
    features = feat.names
  )
  mat <- Reduce(f = `+`, x = mat.list)
  AsBPCells(mat = mat, bpcells = bpcells, bpcells.dir = bpcells.dir)
}

# Validate the bpcells / bpcells.dir arguments.
#
# Errors if `bpcells` is not a single logical, if `bpcells = TRUE` without a
# `bpcells.dir`, if `bpcells.dir` is not a single character, or if
# `bpcells.dir` points to a non-empty directory. Warns if `bpcells.dir` is
# supplied while `bpcells = FALSE`.
ValidateBPCellsArgs <- function(bpcells, bpcells.dir) {
  if (!is.logical(x = bpcells) || length(x = bpcells) != 1 || is.na(x = bpcells)) {
    stop("`bpcells` must be a single logical value")
  }
  if (isTRUE(x = bpcells) && is.null(x = bpcells.dir)) {
    stop(
      "`bpcells.dir` must be supplied when `bpcells = TRUE`. ",
      "Provide a path to a persistent directory where the BPCells output ",
      "will be written."
    )
  }
  if (!is.null(x = bpcells.dir)) {
    if (!is.character(x = bpcells.dir) || length(x = bpcells.dir) != 1) {
      stop("`bpcells.dir` must be a single character string")
    }
    if (dir.exists(paths = bpcells.dir) &&
        length(x = list.files(path = bpcells.dir)) > 0) {
      stop(
        "`bpcells.dir` (", bpcells.dir, ") already exists and is not empty"
      )
    }
  }
  if (isTRUE(x = bpcells) &&
      !requireNamespace("BPCells", quietly = TRUE)) {
    stop(
      "`bpcells = TRUE` requires the BPCells package. Install from ",
      "https://github.com/bnprks/BPCells"
    )
  }
  if (!isTRUE(x = bpcells) && !is.null(x = bpcells.dir)) {
    warning("`bpcells.dir` is ignored when `bpcells = FALSE`")
  }
  invisible(x = NULL)
}

# Persist a matrix to BPCells format and return an IterableMatrix.
#
# If `bpcells = FALSE`, returns `mat` unchanged. Otherwise writes `mat` to
# `bpcells.dir` via `BPCells::write_matrix_dir` and returns the opened
# IterableMatrix.
AsBPCells <- function(mat, bpcells, bpcells.dir) {
  if (!isTRUE(x = bpcells)) {
    return(mat)
  }
  # BPCells::write_matrix_dir requires dgCMatrix specifically (or an existing
  # IterableMatrix). Other CsparseMatrix subclasses produced upstream are
  # coerced here.
  if (!inherits(x = mat, what = c("IterableMatrix", "dgCMatrix"))) {
    mat <- as(object = mat, Class = "CsparseMatrix")
    if (!inherits(x = mat, what = "dgCMatrix")) {
      mat <- as(object = mat, Class = "dgCMatrix")
    }
  }
  BPCells::write_matrix_dir(mat = mat, dir = bpcells.dir)
  BPCells::open_matrix_dir(dir = bpcells.dir)
}

# Run fragtk matrix
#
# Wrapper function to run `fragtk matrix` and return the output as a sparse
# matrix in R.
#
# See <https://crates.io/crates/fragtk> for fragtk documentation.
#
# @param fragments Path to a fragment file
# @param features A GRanges object containing a set of genomic intervals to
# quantify
# @param cells List of cells to include
# @param group Group genomic ranges according to a grouping variable. If FALSE,
# no grouping variable is used. If TRUE, group by the first metadata column.
# If a character string is provided, this should match a column in the
# provided GRanges object supplied in the `features` parameter.
# @param pic Use paired insertion counting
# @param fragtk.path Path to fragtk executable. If NULL, try to find fragtk
# automatically.
# @param seqlevels Named character vector for seqlevels conversion. If
# provided, rename the seqlevels in `features` to match the fragment file
# before quantification, and reverse-map row names in the output.
# @param outdir Path for output directory
# @param cleanup Remove output files created by fragtk
# @param bpcells If TRUE, import the fragtk MTX output into a BPCells
# IterableMatrix instead of a sparse in-memory matrix.
# @param bpcells.dir Directory to write the BPCells matrix to. This directory
# backs the returned matrix and so is not removed by `cleanup`; the caller is
# responsible for removing it. If NULL, a session temporary directory is used.
# @param verbose Display messages
# @return Returns a CsparseMatrix, or a BPCells IterableMatrix if
# `bpcells = TRUE`.
#
#' @importFrom S4Vectors mcols
#' @importFrom Matrix readMM
#' @importFrom GenomeInfoDb renameSeqlevels
#' @keywords internal
RunFragtk <- function(
  fragments,
  features,
  cells,
  group = FALSE,
  pic = TRUE,
  fragtk.path = NULL,
  seqlevels = NULL,
  outdir = tempdir(),
  cleanup = TRUE,
  bpcells = FALSE,
  bpcells.dir = NULL,
  verbose = TRUE
) {
  # find fragtk
  fragtk.path <- fragtk.path %||% unname(obj = Sys.which(names = "fragtk"))
  if (nchar(x = fragtk.path) == 0) {
    stop(
      "fragtk not found. Please install fragtk:",
      "https://crates.io/crates/fragtk"
    )
  }

  if (!dir.exists(paths = outdir)) {
    stop("Requested output directory does not exist")
  }

  # convert seqlevels to match fragment file if needed
  if (!is.null(x = seqlevels)) {
    features <- suppressWarnings(
      expr = renameSeqlevels(x = features, value = seqlevels)
    )
  }

  # temp files
  bed.path <- tempfile(pattern = "signac_fragtk_bed", tmpdir = outdir)
  cells.path <- tempfile(pattern = "signac_fragtk_cells", tmpdir = outdir)
  out.path <- tempfile(pattern = "signac_fragtk_matrix", tmpdir = outdir)

  additional.args <- ""

  # write cells and regions files
  # convert GRanges 1-based closed to 0-based half-open BED format
  feat <- as.data.frame(x = features)[, 1:3]
  feat[, 2] <- feat[, 2] - 1L
  if (is.logical(x = group)) {
    if (group) {
      # passed group=TRUE, assume group by the first metadata column
      if (verbose) {
        message("Grouping regions by column: ", names(mcols(features))[1])
      }
      feat$group <- mcols(x = features)[[1]]
      additional.args <- paste0(additional.args, " --group")
    }
  } else {
    if (is.character(x = group)) {
      # passed column name
      if (!(group %in% names(x = mcols(x = features)))) {
        stop("Requested grouping column '", group, "' does not exist")
      } else {
        if (verbose) {
          message("Grouping regions by column: ", group)
        }
        feat$group <- mcols(x = features)[[group]]
        additional.args <- paste0(additional.args, " --group")
      }
    }
  }

  if (pic) {
    additional.args <- paste0(additional.args, " --pic")
  }
  write.table(
    x = feat,
    file = bed.path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  writeLines(text = cells, con = cells.path)

  # call fragtk
  cmd <- paste0(
    shQuote(string = fragtk.path),
    " matrix --fragments ",
    shQuote(string = fragments),
    " --bed ",
    shQuote(string = bed.path),
    " --cells ",
    shQuote(string = cells.path),
    " --outdir ",
    shQuote(string = out.path),
    " ",
    additional.args
  )

  exit_code <- system(
    command = cmd,
    wait = TRUE,
    ignore.stderr = !verbose,
    ignore.stdout = !verbose
  )
  if (exit_code != 0) {
    stop("fragtk returned a non-zero exit code (", exit_code, ")")
  }

  # read results
  if (verbose) {
    message("Loading count matrix")
  }
  matrix.file <- paste0(out.path, .Platform$file.sep, "matrix.mtx.gz")
  rownames.file <- paste0(out.path, .Platform$file.sep, "features.tsv.gz")
  colnames.file <- paste0(out.path, .Platform$file.sep, "barcodes.tsv.gz")

  if (isTRUE(x = bpcells)) {
    # stream the fragtk MTX directory into an on-disk BPCells matrix. This
    # directory backs the returned matrix, so the caller owns it and is
    # responsible for removing it once the matrix has been re-persisted.
    bpcells.tmp <- bpcells.dir %||%
      tempfile(pattern = "signac_fragtk_bpcells_")
    imported <- BPCells::import_matrix_market_10x(mtx_dir = out.path)
    BPCells::write_matrix_dir(mat = imported, dir = bpcells.tmp)
    counts <- BPCells::open_matrix_dir(dir = bpcells.tmp)
    # BPCells sets colnames from barcodes.tsv and rownames from features.tsv
    # (first column). Replicate the second-column row name (interval string)
    # so it matches the behavior of the sparse code path.
    rownames(x = counts) <- readLines(con = rownames.file)
    colnames(x = counts) <- readLines(con = colnames.file)
  } else {
    counts <- readMM(file = matrix.file)
    rownames(x = counts) <- readLines(con = rownames.file)
    colnames(x = counts) <- readLines(con = colnames.file)
    counts <- as(object = counts, Class = "CsparseMatrix")
  }

  # reverse-map seqlevels in row names back to original names
  if (!is.null(x = seqlevels)) {
    sl <- names(x = seqlevels)
    names(x = sl) <- seqlevels
    mapped.features <- suppressWarnings(
      expr = renameSeqlevels(x = features, value = sl)
    )
    feat.str <- as.character(x = mapped.features)
    # only remap row names for non-grouped output
    if (is.logical(x = group) && !group) {
      rownames(x = counts) <- feat.str
    }
  }

  # remove temp files
  if (cleanup) {
    files.to.remove <- c(
      cells.path,
      bed.path,
      matrix.file,
      rownames.file,
      colnames.file
    )
    for (i in files.to.remove) {
      if (file.exists(i)) {
        file.remove(i)
      }
    }
    unlink(x = out.path, recursive = TRUE)
  }
  return(counts)
}

#### Not Exported ####

# Collapse matrix rows by group labels
#
# Sum rows of a sparse matrix according to a grouping variable, using
# matrix multiplication (same pattern as CombineTiles).
#
# @param mat A sparse matrix
# @param features A GRanges object with metadata columns
# @param group If TRUE, group by the first metadata column. If a character
# string, group by the named metadata column. If FALSE, return mat unchanged.
# @return Returns a sparse matrix with rows collapsed by group
#' @importFrom S4Vectors mcols
#' @importFrom Matrix sparseMatrix crossprod
#' @importMethodsFrom Matrix t
GroupMatrix <- function(mat, features, group) {
  if (is.logical(x = group) && !group) {
    return(mat)
  }
  if (is.logical(x = group) && group) {
    group_labels <- mcols(x = features)[[1]]
  } else if (is.character(x = group)) {
    if (!(group %in% names(x = mcols(x = features)))) {
      stop("Requested grouping column '", group, "' does not exist")
    }
    group_labels <- mcols(x = features)[[group]]
  } else {
    return(mat)
  }
  unique_groups <- unique(x = group_labels)
  group_idx <- match(x = group_labels, table = unique_groups)
  collapse <- sparseMatrix(
    i = seq_along(along.with = group_labels),
    j = group_idx,
    x = 1,
    dims = c(length(x = group_labels), length(x = unique_groups))
  )
  collapsed <- crossprod(x = mat, y = collapse)
  collapsed <- t(x = collapsed)
  rownames(x = collapsed) <- unique_groups
  return(as(object = collapsed, Class = "CsparseMatrix"))
}

# matrix multiplication method for summing matrix rows
#' @importFrom GenomicRanges reduce
#' @importFrom S4Vectors elementNROWS
#' @importFrom Matrix crossprod sparseMatrix
#' @importMethodsFrom Matrix t
CombineTiles <- function(bins) {
  ranges <- GRanges(rownames(x = bins))
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
  rownames(x = collapsed) <- as.character(x = reduced.tiles)

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
  path,
  features,
  file.index = NULL,
  cells = NULL,
  pic = TRUE,
  group = FALSE,
  keep_all_features = FALSE,
  process_n = 2000,
  verbose = TRUE
) {
  # locate the index automatically (handles both .tbi and .csi) when not given
  file.index <- file.index %||% GetIndexFile(fragment = path, verbose = verbose)
  feat.use <- features
  tbx <- TabixFile(file = path, index = file.index)
  n_feat_start <- length(x = feat.use)
  if (keep_all_features) {
    features_to_get <- as.character(x = feat.use)
  } else {
    features_to_get <- NULL
  }
  feat.use <- keepSeqlevels(
    x = feat.use,
    value = intersect(
      x = seqnames(x = feat.use),
      y = seqnamesTabix(file = tbx)
    ),
    pruning.mode = "coarse"
  )
  if (length(x = feat.use) == 0) {
    stop("No matching chromosomes found in fragment file.")
  }
  n_removed <- n_feat_start - length(x = feat.use)
  if (n_removed > 0 && !keep_all_features) {
    if (n_removed == 1) {
      warning(
        n_removed, " feature is on a seqname not present in ",
        "the fragment file. This will be removed."
      )
    } else {
      warning(
        n_removed, " features are on seqnames not present in ",
        "the fragment file. These will be removed."
      )
    }
  }
  feature.list <- ChunkGRanges(
    granges = feat.use,
    nchunk = ceiling(x = length(x = feat.use) / process_n)
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
      pic = pic,
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
      pic = pic
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
      FUN = AddMissing,
      cells = all.cells,
      features = NULL
    )
  }
  featmat <- do.call(what = rbind, args = matrix.parts)
  # add zero rows for features that were not quantified, and reorder features
  if (keep_all_features) {
    featmat <- AddMissing(
      x = featmat, cells = NULL, features = features_to_get
    )
    feat.str <- features_to_get
  } else {
    feat.str <- as.character(x = feat.use)
  }
  featmat <- featmat[feat.str, , drop = FALSE]
  # apply grouping if requested
  if (!is.logical(x = group) || group) {
    featmat <- GroupMatrix(mat = featmat, features = feat.use, group = group)
  }
  return(featmat)
}

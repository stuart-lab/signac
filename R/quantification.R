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
  assay.obj <- CreateChromatinAssay5(counts = bins, fragments = frags)
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
#' `fragtk` executable and `fragtk` will be used. Note that
#' `fragtk` uses the Paired Insertion Counting method, whereas the R
#' implementation counts insertions. See
#' <https://crates.io/crates/fragtk> for fragtk documentation.
#' @param gene.id Record gene IDs in output matrix rather than gene name.
#' @param verbose Display messages
#'
#' @return Returns a sparse matrix
#'
#' @concept utilities
#' @export
#' @importFrom SeuratObject DefaultAssay
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' GeneActivity(atac_small, fragtk=FALSE)
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
    gene.id = FALSE,
    verbose = TRUE
) {
  if (!is.null(x = features)) {
    if (length(x = features) == 0) {
      stop("Empty list of features provided")
    }
  }
  # collapse to longest protein coding transcript
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("The requested assay is not a ChromatinAssay5.")
  }
  annotation <- Annotation(object = object[[assay]])
  # replace NA names with gene ID
  annotation$gene_name <- ifelse(
    test = is.na(x = annotation$gene_name) | (annotation$gene_name == ""),
    yes = annotation$gene_id,
    no = annotation$gene_name
  )
  if (length(x = annotation) == 0) {
    stop("No gene annotations present in object")
  }
  if (verbose) {
    message("Extracting gene coordinates")
  }
  transcripts <- CollapseToLongestTranscript(ranges = annotation)
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
  
  # quantify
  frags <- Fragments(object = object[[assay]])
  if (length(x = frags) == 0) {
    stop("No fragment information found for requested assay")
  }
  cells <- colnames(x = object[[assay]])
  counts <- FeatureMatrix(
    fragments = frags,
    features = transcripts,
    process_n = process_n,
    cells = cells,
    fragtk = fragtk,
    verbose = verbose
  )
  # set row names
  gene.key <- transcripts$gene_name
  names(x = gene.key) <- as.character(x = transcripts)
  rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
  counts <- counts[rownames(x = counts) != "", ]
  
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
#' @param verbose Display messages.
#' @param ... Arguments passed to [FeatureMatrix()].
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
  verbose = TRUE,
  ...
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
    verbose = verbose,
    ...
  )
  return(binmat)
}

#' Feature Matrix
#'
#' Construct a feature x cell matrix from a genomic fragments file
#'
#' @param fragments A list of [Fragment()] objects. Note that if
#' setting the `cells` parameter, the requested cells should be present in
#' the supplied `Fragment` objects. However, if the cells information in
#' the fragment object is not set (`Cells(fragments)` is `NULL`), then
#' the fragment object will still be searched.
#' @param features A GRanges object containing a set of genomic intervals.
#' These will form the rows of the matrix, with each entry recording the number
#' of unique reads falling in the genomic region for each cell.
#' @param fragtk Use `fragtk` for fast and memory-efficient data
#' quantification. Can be TRUE/FALSE or a character vector. If TRUE,
#' `fragtk` will be used and attempt to find the `fragtk` executable
#' in the path. If FALSE, use the R implementation to produce the data matrix.
#' If a character vector is provided, this should be the path to the
#' `fragtk` executable and `fragtk` will be used. Note that
#' `fragtk` uses the Paired Insertion Counting method, whereas the R
#' implementation counts insertions. See
#' <https://crates.io/crates/fragtk> for fragtk documentation. The R
#' implementation (`fragtk=FALSE`) will be faster for quantifying a small
#' number of genomic regions.
#' @param keep_all_features By default, if a genomic region provided is on a
#' chromosome that is not present in the fragment file,
#' it will not be included in the returned matrix. Set `keep_all_features` to
#' TRUE to force output to include all features in the input ranges. Note that
#' features on chromosomes that are not present in the fragment file will be 
#' filled with zero counts.
#' @param cells Vector of cells to include. If NULL, include all cells found
#' in the fragments file
#' @param process_n Number of regions to load into memory at a time, per thread.
#' Processing more regions at once can be faster but uses more memory.
#' @param verbose Display messages
#'
#' @export
#' @importFrom SeuratObject RowMergeSparseMatrices
#' @importFrom GenomeInfoDb renameSeqlevels
#' @importFrom GenomicRanges GRanges
#' @concept quantification
#' @return Returns a sparse matrix
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(fpath)
#' FeatureMatrix(
#'   fragments = fragments,
#'   features = granges(atac_small),
#'   fragtk = FALSE
#' )
FeatureMatrix <- function(
  fragments,
  features,
  fragtk = TRUE,
  keep_all_features = FALSE,
  cells = NULL,
  process_n = 2000,
  verbose = TRUE
) {
  if (!inherits(x = features, what = "GRanges")) {
    if (inherits(x = features, what = "character")) {
      features <- GRanges(features)
    } else {
      stop("features should be a GRanges object")
    }
  }
  if (!inherits(x = fragments, what = "list")) {
    if (inherits(x = fragments, what = "Fragment2")) {
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
      if (is.null(x = Cells(fragments[[i]]))) {
        # cells information not set for fragment object
        obj.use <- c(obj.use, i)
      } else if (any(cells %in% Cells(x = fragments[[i]]))) {
        obj.use <- c(obj.use, i)
      }
    }
  } else {
    obj.use <- seq_along(along.with = fragments)
  }
  # create a matrix from each fragment file
  fragtk.path <- NULL
  if (is.character(x = fragtk)) {
    # fragtk is the path to executable
    fragtk.path <- fragtk
    fragtk <- TRUE
  }
  if (fragtk) {
    # run fragtk on each fragment file
    # update cell names in output matrix
    mat.list <- sapply(
      X = obj.use,
      FUN = function(x) {
        cell.vec <- GetFragmentData(
          object = fragments[[x]],
          slot = "cells"
        )
        seqlevel.conversion <- GetFragmentData(
          object = fragments[[x]],
          slot = "seqlevels"
        )
        
        if (!is.null(x = seqlevel.conversion)) {
          # replace seqnames
          feat.use <- suppressWarnings(
            expr = renameSeqlevels(x = features, value = seqlevel.conversion)
          )
        } else {
          feat.use <- features
        }
        
        mat <- RunFragtk(
          fragments = GetFragmentData(
            object = fragments[[x]],
            slot = "file.path"
          ),
          features = feat.use,
          cells = cell.vec,
          fragtk.path = fragtk.path,
          verbose = verbose,
          cleanup = TRUE
        )
        colnames(x = mat) <- names(x = cell.vec)
        rownames(x = mat) <- as.character(x = features)
        mat
      }
    )
  } else {
    mat.list <- sapply(
      X = obj.use,
      FUN = function(x) {
        SingleFeatureMatrix(
          fragment = fragments[[x]],
          features = features,
          keep_all_features = keep_all_features,
          cells = cells,
          verbose = verbose,
          process_n = process_n
        )
      })
  }

  # merge all the matrices
  if (length(x = mat.list) == 1) {
    return(mat.list[[1]])
  } else {
    # ensure all cells and features present, same order
    all.cells <- unique(
      x = unlist(x = lapply(X = mat.list, FUN = colnames))
    )
    mat.list <- lapply(
      X = mat.list,
      FUN = AddMissing,
      cells = all.cells,
      features = as.character(x = features)
    )
    featmat <- Reduce(f = `+`, x = mat.list)
    return(featmat)
  }
}

#' Run fragtk matrix
#' 
#' Wrapper function to run `fragtk matrix` and return the output as a sparse
#' matrix in R.
#' 
#' See <https://crates.io/crates/fragtk> for fragtk documentation.
#' 
#' @param fragments A list of Fragment objects or fragment file paths
#' @param features A GRanges object containing a set of genomic intervals to
#' quantify. These genomic ranges will be passed to the `--bed` argument in
#' `fragtk matrix`
#' @param cells List of cells to include
#' @param group Group genomic ranges according to a grouping variable. If NULL,
#' no grouping variable is used. If a character string is provided, this should
#' match a column in the provided GRanges object supplied in the `features`
#' parameter.
#' @param pic Use paired insertion counting
#' @param fragtk.path Path to fragtk executable. If NULL, try to find fragtk
#' automatically.
#' @param outdir Path for output directory
#' @param cleanup Remove output files created by fragtk
#' @param verbose Display messages
#' 
#' @importFrom S4Vectors mcols
#' @importFrom Matrix readMM
#' 
#' @concept quantification
#' 
#' @return Returns a CsparseMatrix
#' @export
RunFragtk <- function(
    fragments,
    features,
    cells,
    group = FALSE,
    pic = TRUE,
    fragtk.path = NULL,
    outdir = tempdir(),
    cleanup = TRUE,
    verbose = TRUE
) {
  # find fragtk
  fragtk.path <- fragtk.path %||% unname(obj = Sys.which(names = "fragtk"))
  if (nchar(x = fragtk.path) == 0) {
    stop("fragtk not found. Please install fragtk:",
         "https://crates.io/crates/fragtk")
  }
  
  if (!dir.exists(paths = outdir)) {
    stop("Requested output directory does not exist")
  }
  
  # temp files
  bed.path <- tempfile(pattern = "signac_fragtk_bed", tmpdir = outdir)
  cells.path <- tempfile(pattern = "signac_fragtk_cells", tmpdir = outdir)
  out.path <- tempfile(pattern = "signac_fragtk_matrix", tmpdir = outdir)
  
  additional.args <- ""
  
  # write cells and regions files
  feat <- as.data.frame(features)[, 1:3]
  if (is.logical(x = group)) {
    if (group) {
      # passed group=TRUE, assume group by the 4th bed column
      if (verbose) {
        message("Grouping regions by column: ", names(mcols(features))[1])
      }
      feat$group <- mcols(features)[[1]]
      additional.args <- paste0(additional.args, " --group")
    }
  } else {
    if (is.character(x = group)) {
      # passed column name
      if (!(group %in% names(mcols(features)))) {
        stop("Requested grouping column '", group, "' does not exist")
      } else {
        if (verbose) {
          message("Grouping regions by column: ", group)
        }
        feat$group <- mcols(features)[[group]]
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
    fragtk.path,
    " matrix --fragments ",
    fragments,
    " --bed ",
    bed.path,
    " --cells ",
    cells.path,
    " --outdir ",
    out.path,
    " ",
    additional.args
  )
  
  system(
    command = cmd,
    wait = TRUE,
    ignore.stderr = !verbose,
    ignore.stdout = !verbose
  )
  
  # read results
  if (verbose) {
    message("Loading count matrix")
  }
  matrix.file <- paste0(out.path, .Platform$file.sep, "matrix.mtx.gz")
  rownames.file <- paste0(out.path, .Platform$file.sep, "features.tsv.gz")
  colnames.file <- paste0(out.path, .Platform$file.sep, "barcodes.tsv.gz")
  
  counts <- readMM(file = matrix.file)
  rownames(counts) <- readLines(rownames.file)
  colnames(counts) <- readLines(colnames.file)
  counts <- as(object = counts, Class = "CsparseMatrix")
  
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

# matrix multiplication method for summing matrix rows
#' @importFrom GenomicRanges reduce GRanges
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
#' @importMethodsFrom GenomicRanges intersect GRanges
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom fastmatch fmatch
SingleFeatureMatrix <- function(
  fragment,
  features,
  keep_all_features = FALSE,
  cells = NULL,
  process_n = 2000,
  verbose = TRUE
) {
  fragment.path <- GetFragmentData(object = fragment, slot = "file.path")
  fragment.index <- GetFragmentData(object = fragment, slot = "file.index")
  frag.cells <- GetFragmentData(object = fragment, slot = "cells")
  seqlevel.conversion <- GetFragmentData(object = fragment, slot = "seqlevels")
  
  # rename seqlevels of features to match seqlevels in the fragment file
  # afterwards, map back to the seqlevels in the input features
  if (!is.null(x = seqlevel.conversion)) {
    # replace seqnames with names in the fragment file
    feat.use <- suppressWarnings(
      expr = renameSeqlevels(x = features, value = seqlevel.conversion)
    )
  } else {
    feat.use <- features
  }
  
  if (!is.null(cells)) {
    # only look for cells that are in the fragment file
    if (is.null(x = frag.cells)) {
      # cells information not set in fragment object
      names(x = cells) <- cells
    } else {
      # first subset frag.cells
      cell.idx <- fmatch(
        x = names(x = frag.cells),
        table = cells,
        nomatch = 0L
      ) > 0
      cells <- frag.cells[cell.idx]
    }
  } else {
    # cells not set, but still need to only return cells that 
    # are in the fragment objects, with cell names converted
    if (!is.null(x = frag.cells)) {
      cells <- frag.cells
    }
  }
  tbx <- TabixFile(file = fragment.path, index = fragment.index)
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
  if (n_removed > 0 && ! keep_all_features) {
    if (n_removed == 1) {
      warning(n_removed, " feature is on a seqname not present in ",
              "the fragment file. This will be removed.")
    } else {
      warning(n_removed, " features are on seqnames not present in ",
              "the fragment file. These will be removed.")
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
      future.globals = list(),
      future.scheduling = FALSE
    )
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    matrix.parts <- mylapply(
      X = feature.list,
      FUN = PartialMatrix,
      tabix = tbx,
      cells = cells
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
      features = features_to_get
    )
  } else if (keep_all_features) {
    matrix.parts <- lapply(
      X = matrix.parts,
      FUN = AddMissing,
      cells = NULL,
      features = features_to_get
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
  if (keep_all_features) {
    feat.str <- features_to_get
  } else {
    feat.str <- as.character(x = feat.use)
  }
  featmat <- featmat[feat.str, , drop=FALSE]
  if (!is.null(x = seqlevel.conversion)) {
    # map back to original seqnames
    sl <- names(x = seqlevel.conversion)
    names(x = sl) <- seqlevel.conversion
    features <- suppressWarnings(
      expr = renameSeqlevels(x = feat.use, value = sl)
    )
    feat.str <- as.character(x = features)
    rownames(x = featmat) <- feat.str
  }
  return(featmat)
}

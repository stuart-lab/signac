#' @include generics.R
#'
NULL

#' Return the first rows of a fragment file
#'
#' Returns the first \code{n} rows of a fragment file. This allows the content
#' of a fragment file to be inspected.
#'
#' @param x a \code{Fragment2} object
#' @param n an integer specifying the number of rows to return from the fragment
#' file
#' @param ... additional arguments passed to \code{\link[utils]{read.table}}
#'
#' @return The first \code{n} rows of a fragment file as a \code{data.frame}
#' with the following columns: chrom, start, end, barcode, readCount.
#'
#' @export
#' @method head Fragment2
#' @concept fragments
#'
head.Fragment2 <- function(x, n = 6L, ...) {
  fpath <- GetFragmentData(object = x, slot = "file.path")
  df <- read.table(file = fpath, nrows = n, ...)
  if (ncol(x = df) == 5) {
    colnames(x = df) <- c("chrom", "start", "end", "barcode", "readCount")
  } else if(ncol(x = df) == 6) {
    colnames(x = df) <- c("chrom", "start", "end", "barcode", "readCount", "strand")
  }
  return(df)
}

#' Return the fragment file header
#' 
#' Returns any comment lines present at the start of the fragment file,
#' before any fragment entries.
#' 
#' @param x a \code{Fragment2} object
#' @return Returns a character vector with the fragment file header lines
#' @export
#' @concept fragments
header <- function(x) {
  fpath <- GetFragmentData(object = x, slot = "file.path")
  con <- file(fpath, "r")
  header_lines <- character()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0 || !startsWith(line, "#")) break
    header_lines <- c(header_lines, line)
  }
  close(con)
  return(header_lines)
}

#' Count fragments
#'
#' Count total fragments per cell barcode present in a fragment file.
#'
#' @param fragments Path to a fragment file. If a list of fragment files is
#' provided, the total fragments for each cell barcode across all files will be
#' returned
#' @param cells Cells to include. If NULL, include all cells
#' @param max_lines Maximum number of lines to read from the fragment file. If
#' NULL, read all lines in the file.
#' @param verbose Display messages
#'
#' @rdname CountFragments
#' @export
#' @concept fragments
#' @return Returns a data.frame with the following columns:
#' \itemize{
#'   \item{CB: the cell barcode}
#'   \item{frequency_count: total number of fragments sequenced for the cell}
#'   \item{mononucleosome: total number of fragments with length between 147 bp and 294 bp}
#'   \item{nucleosome_free: total number of fragments with length <147 bp}
#'   \item{reads_count: total number of reads sequenced for the cell}
#' }
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' counts <- CountFragments(fragments = fpath)
CountFragments <- function(
  fragments,
  cells = NULL,
  max_lines = NULL,
  verbose = TRUE
) {
  if (!inherits(x = fragments, what = "list")) {
    fragments <- list(fragments)
  }
  for (i in seq_along(along.with = fragments)) {
    if (isRemote(x = i)) {
      stop("Remote fragment files not supported")
    }
    fragments[[i]] <- normalizePath(path = fragments[[i]], mustWork = TRUE)
    max_lines <- SetIfNull(x = max_lines, y = 0)
    verbose = as.logical(x = verbose)
    if (!is.null(x = cells)) {
      cells <- unique(x = cells)
    }
    counts <- groupCommand(
      fragments = fragments[[i]],
      some_whitelist_cells = cells,
      max_lines = max_lines,
      verbose = verbose
    )
    rownames(x = counts) <- counts$CB
    counts$CB <- NULL
    if (i == 1) {
      # first file
      allcounts <- counts
    } else {
      # merge
      common <- intersect(
        x = rownames(x = allcounts), y = rownames(x = counts)
      )
      allcounts[common, ] <- allcounts[common, ] + counts[common, ]
      missing_cells <- setdiff(x = rownames(x = counts), y = common)
      allcounts <- rbind(allcounts, counts[missing_cells, ])
    }
  }
  # reformat for backwards compatibility
  allcounts$CB <- rownames(x = allcounts)
  rownames(x = allcounts) <- NULL
  allcounts <- allcounts[, c(5, 1, 2, 3, 4)]
  return(allcounts)
}

#' Filter cells from fragment file
#'
#' Remove all fragments that are not from an allowed set of cell barcodes from
#' the fragment file. This will create a new file on disk that only contains
#' fragments from cells specified in the \code{cells} argument. The output file
#' is block gzip-compressed and indexed, ready for use with Signac functions.
#'
#' @param fragments Path to a fragment file
#' @param cells A vector of cells to keep
#' @param outfile Name for output file
#' @param buffer_length Size of buffer to be read from the fragment file. This
#' must be longer than the longest line in the file.
#' @param verbose Display messages
#'
#' @importFrom Rsamtools bgzip indexTabix
#'
#' @export
#' @concept fragments
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' tmpf <- tempfile(fileext = ".gz")
#' FilterCells(
#'   fragments = fpath,
#'   cells = head(colnames(atac_small)),
#'   outfile = tmpf
#' )
#' file.remove(tmpf)
FilterCells <- function(
  fragments,
  cells,
  outfile = NULL,
  buffer_length = 256L,
  verbose = TRUE
) {
  if (isRemote(x = fragments)) {
    stop("Remote fragment files not supported")
  }
  fragments <- normalizePath(path = fragments, mustWork = TRUE)
  if (is.null(x = outfile)) {
    outfile <- paste0(fragments, ".filtered")
  }
  if (file.exists(outfile)) {
    warning("Output file already exists, file will be overwritten",
            immediate. = TRUE)
  }
  verbose <- as.logical(x = verbose)
  tf <- tempfile(pattern = "filtercells", tmpdir = tempdir(), fileext = "")
  buffer_length <- as.integer(x = buffer_length)
  filtered <- filterCells(
    fragments = fragments,
    keep_cells = cells,
    outfile = tf,
    buffer_length = buffer_length,
    verbose = verbose
  )
  if (filtered == 1) {
    stop("Error: cannot open requested file")
  }
  # bgzip and index output
  if (verbose) {
    message("\nCompressing filtered file")
  }
  bgzip(file = tf, dest = outfile, overwrite = TRUE)
  file.remove(tf)
  if (verbose) {
    message("Indexing fragment file")
  }
  idx <- indexTabix(file = outfile, format = "bed")
}

#' Split fragment file by cell identities
#'
#' Splits a fragment file into separate files for each group of cells. If
#' splitting multiple fragment files containing common cell types, fragments
#' originating from different files will be appended to the same file for one
#' group of cell identities.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param group.by Name of grouping variable to group cells by
#' @param idents List of identities to include
#' @param buffer_length Size of buffer to be read from the fragment file. This
#' must be longer than the longest line in the file.
#' @param outdir Directory to write output files
#' @param file.suffix Suffix to add to all file names (before file extension).
#' If splitting multiple fragment files without the \code{append} option set to
#' TRUE, an additional numeric suffix will be added to each file (eg, .1, .2).
#' @param append If splitting multiple fragment files, append cells from the
#' same group (eg cluster) to the same file. Note that this can cause the output
#' file to be unsorted.
#' @param verbose Display messages
#'
#' @importFrom SeuratObject DefaultAssay
#' @concept fragments
#'
#' @export
SplitFragments <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  outdir = getwd(),
  file.suffix = "",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  frags <- Fragments(object = object[[assay]])
  groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  # replace space with underscore
  cells <- names(x = groups)
  groups <- gsub(pattern = " ", replacement = "_", x = groups)
  groups <- gsub(pattern = .Platform$file.sep, replacement = "_", x = groups)
  names(x = groups) <- cells
  buffer_length <- as.integer(x = buffer_length)
  file.suffix <- as.character(x = file.suffix)
  idents <- as.character(x = unname(obj = groups))
  unique_idents <- unique(x = idents)
  outdir <- normalizePath(path = outdir, mustWork = TRUE)

  # split cells from each fragment file
  # append to existing file when more than one fragment file used
  for (i in seq_along(along.with = frags)) {
    if (!append & (length(x = frags) > 1)) {
      suffix.use <- paste0(file.suffix, ".", i)
    } else {
      suffix.use <- file.suffix
    }
    fragpath <- GetFragmentData(object = frags[[i]], slot = "file.path")
    if (isRemote(x = fragpath)) {
      message("Remote fragment files not supported, skipping fragment file")
    } else {
      # convert cell names
      cellmap <- GetFragmentData(object = frags[[i]], slot = "cells")
      cell.in.frag <- cells %in% names(x = cellmap)
      cells.use <- cellmap[cells[cell.in.frag]]
      idents.use <- as.character(x = unname(obj = groups[names(x = cells.use)]))
      frag.cell.name <- as.character(x = unname(obj = cells.use))

      if (verbose) {
        message("Processing file ", fragpath)
      }
      splitfiles <- splitFragments(
        fragments = fragpath,
        outdir = paste0(outdir, .Platform$file.sep),
        suffix = suffix.use,
        append = append,
        cells = frag.cell.name,
        idents = idents.use,
        unique_idents = unique_idents,
        buffer_length = buffer_length,
        verbose = verbose
      )
      if (verbose) {
        message("\n")
      }
      if (splitfiles == 1) {
        stop("Error: cannot open requested file")
      }
    }
  }
}

#' Validate cells present in fragment file
#'
#' Search for a fragment from each cell that should exist in the fragment file.
#' Will iterate through chunks of the fragment file until at least one fragment
#' from each cell barcode requested is found.
#'
#' @param object A \code{\link{Fragment}} object
#' @param cells A character vector containing cell barcodes to search for.
#' If NULL, use the cells stored in the Fragment object.
#' @param tolerance Fraction of input cells that can be unseen before returning
#' TRUE. For example, \code{tolerance = 0.01} will return TRUE when 99% of cells
#' have observed fragments in the file. This can be useful if there are cells
#' present that have much fewer total counts, and would require extensive
#' searching before a fragment from those cells are found.
#' @param max.lines Maximum number of lines to read in without finding the
#' required number of cells before returning FALSE. Setting this value avoids
#' having to search the whole file if it becomes clear that the expected cells
#' are not present. Setting this value to NULL will enable an exhaustive search
#' of the entire file.
#' @param verbose Display messages
#' @export
#' @concept fragments
ValidateCells <- function(
  object,
  cells = NULL,
  tolerance = 0.5,
  max.lines = 5e7,
  verbose = TRUE
) {
  cell_barcodes <- GetFragmentData(object = object, slot = "cells")
  cells <- SetIfNull(x = cells, y = cell_barcodes)
  if (is.null(x = cells)) {
    warning("No cells stored in object")
    return(TRUE)
  }
  max.lines <- SetIfNull(x = max.lines, y = 0)
  filepath <- GetFragmentData(object = object, slot = "file.path")
  filepath <- normalizePath(path = filepath, mustWork = TRUE)
  is.remote <- isRemote(x = filepath)
  find_n <- as.integer(x = length(x = cells) * (1 - tolerance))
  # if remote, return TRUE
  if (is.remote) {
    return(TRUE)
  }
  valid <- validateCells(
    fragments = filepath,
    cells = cells,
    find_n = find_n,
    max_lines = max.lines,
    verbose = verbose
  )
  return(valid)
}

#' Validate hashes for Fragment object
#'
#' @param object A \code{\link{Fragment}} object
#' @param verbose Display messages
#' @export
#' @concept fragments
#' @importFrom tools md5sum
ValidateHash <- function(object, verbose = TRUE) {
  path <- GetFragmentData(object = object, slot = "file.path")
  index.file <- GetFragmentData(object = object, slot = "file.index")
  is.remote <- isRemote(x = path)
  # if remote, return TRUE
  if (is.remote) {
    return(TRUE)
  }
  if (!all(file.exists(path, index.file))) {
    return(FALSE)
  }
  oldhash <- GetFragmentData(object = object, slot = "hash")
  if (verbose) {
    message("Computing hash")
  }
  newhash <- md5sum(files = c(path, index.file))
  return(all(newhash == oldhash))
}

#' Validate Fragment object
#'
#' Verify that the cells listed in the object exist in the fragment file
#' and that the fragment file or index have not changed since creating the
#' fragment object.
#'
#' @param object A \code{\link{Fragment}} object
#' @param verbose Display messages
#' @param ... Additional parameters passed to \code{\link{ValidateCells}}
#' @export
#' @concept fragments
ValidateFragments <- function(
  object,
  verbose = TRUE,
  ...
) {
  valid.cells <- ValidateCells(object = object, verbose = verbose, ...)
  valid.hash <- ValidateHash(object = object, verbose = verbose, ...)
  return(valid.cells & valid.hash)
}

#' Set and get cell barcode information for a \code{\link{Fragment}} object
#'
#' This returns the names of cells in the object that are contained in the
#' fragment file. These cell barcodes may not match the barcodes present in the
#' fragment file. The \code{\link{Fragment}} object contains an internal mapping
#' of the cell names in the \code{\link{ChromatinAssay}} object to the cell
#' names in the fragment file, so that cell names can be changed in the
#' assay without needing to change the cell names on disk.
#'
#' To access the cell names that are stored in the fragment file itself, use
#' \code{GetFragmentData(object = x, name = "cells")}.
#' @param x A Fragment object
#' @param ... Arguments passed to other methods
#' @rdname Cells
#' @concept fragments
#' @method Cells Fragment2
#' @importFrom SeuratObject Cells
#' @export
Cells.Fragment2 <- function(x, ...) {
  cells <- slot(object = x, name = "cells")
  return(names(x = cells))
}

# Re-export SeuratObject generic
#' @importFrom SeuratObject Cells
#' @export
SeuratObject::Cells

#' @param value A vector of cell names to store in the \code{\link{Fragment}}
#' object
#' @rdname Cells
#' @export
#' @concept fragments
#' @method Cells<- Fragment2
"Cells<-.Fragment2" <- function(x, ..., value) {
  if (is.null(x = names(x = value))) {
    stop("Cells must be a named vector")
  }
  slot(object = x, name = "cells") <- value
  if (!ValidateCells(object = x, verbose = FALSE, ...)) {
    stop("Cells not present in fragment file")
  } else {
    return(x)
  }
}

#' @importFrom Seqinfo seqlevels
#' @exportMethod seqlevels
setMethod(
  f = "seqlevels",
  signature = "Fragment2",
  definition = function(x) {
    return(names(x = slot(object = x, name = "seqlevels")))
  }
)

#' @importFrom Seqinfo seqlevels<-
setReplaceMethod(
  f = "seqlevels",
  signature = "Fragment2",
  definition = function(x, value) {
    if (is.null(x = value)) {
      x@seqlevels <- value
    } else if (!inherits(x = value, what = "character")) {
      stop("New seqlevels must be a character vector")
    } else if (is.null(x = names(x = value))) {
      stop("New seqlevels must be a named vector")
    } else {
      x@seqlevels <- value
    }
    return(x)
  }
)

#' @importFrom GenomeInfoDb seqlevelsStyle
#' @exportMethod seqlevelsStyle
setMethod(
  f = "seqlevelsStyle",
  signature = "Fragment2",
  definition = function(x) {
    return(seqlevelsStyle(seqlevels(x = x)))
  }
)

#' @importFrom GenomeInfoDb seqlevelsStyle<-
setReplaceMethod(
  f = "seqlevelsStyle",
  signature = "Fragment2",
  definition = function(x, value) {
    sl <- seqlevels(x = x)
    if (is.null(x = sl)) {
      stop("seqlevels information not set")
    } else {
      seqlevelsStyle(sl) <- value
      current.vec <- x@seqlevels
      names(x = current.vec) <- sl
      x@seqlevels <- current.vec
    }
    return(x)
  }
)

# map new sequence levels to the old (names correspond to the old levels)
#' @importFrom GenomeInfoDb renameSeqlevels
#' @exportMethod renameSeqlevels
setMethod(
  f = "renameSeqlevels",
  signature = "Fragment2",
  definition = function(x, value) {
    current.vec <- x@seqlevels
    if (is.null(x = current.vec)) {
      stop("seqlevels information not set")
    }
    if (is.null(x = names(x = value))) {
      # unnamed vector
      # needs to be the same length as the current seqlevels
      if (length(x = value) != length(x = current.vec)) {
        stop("Must provide the same number of new seqlevels as existing seqlevels")
      } else {
        names(x = current.vec) <- value
      }
    } else {
      # map names
      old.names <- names(x = current.vec)
      names(x = old.names) <- old.names
      old.names[names(value)] <- value
      names(x = current.vec) <- old.names
    }
    x@seqlevels <- current.vec
    return(x)
  }
)

#' Update the file path for a Fragment object
#'
#' Change the path to a fragment file store in a \code{\link{Fragment}}
#' object. Path must be to the same file that was used to create the fragment
#' object. An MD5 hash will be computed using the new path and compared to the
#' hash stored in the Fragment object to verify that the files are the same.
#'
#' @param object A \code{\link{Fragment}} object.
#' @param new.path Path to the fragment file.
#' @param new.index.path Path to the fragment file index. If NULL, the index is
#' assumed to be in the same directory as the fragment file.
#' @param verbose Display messages.
#'
#' @concept fragments
#' @export
UpdatePath <- function(object, new.path, new.index.path = NULL, verbose = TRUE) {
  new.is.remote <- isRemote(x = new.path)
  if (!new.is.remote) {
    new.path <- normalizePath(path = new.path, mustWork = TRUE)
    index.file <- SetIfNull(
      x = new.index.path,
      y = GetIndexFile(fragment = new.path, verbose = verbose)
    )
    if (!file.exists(new.path)) {
      stop("Fragment file not found")
    } else if (!file.exists(index.file)) {
      stop("Fragment file not indexed")
    }
  }
  old.path <- GetFragmentData(object = object, slot = "file.path")
  old.is.remote <- isRemote(x = old.path)
  if (identical(x = old.path, y = new.path)) {
    return(object)
  }
  if (!old.is.remote & new.is.remote) {
    warning("Replacing local file path with a remote file")
  }
  slot(object = object, name = "file.path") <- new.path
  slot(object = object, name = "file.index") <- index.file
  if (ValidateHash(object = object, verbose = verbose)) {
    return(object)
  } else {
    stop("MD5 sum does not match previously computed sum")
  }
}

# If Fragment object does not contain cell names, assign to given list of names
# Used in CreateChromatinAssay when adding a Fragment object with no cell names
#' @importFrom methods slot slot<-
AssignFragCellnames <- function(fragments, cellnames) {
  if (is.null(x = Cells(x = fragments))) {
    cells <- cellnames
    names(x = cells) <- cells
    slot(object = fragments, name = "cells") <- cells
  }
  return(fragments)
}

# Get index file path
# checks csi and tbi extensions
# tbi used if both present
# @param fragment fragment file path
# @param verbose display messages
# @return returns the index file path
GetIndexFile <- function(fragment, verbose = TRUE) {
  is.remote <- isRemote(x = fragment)
  index.filepaths <- c(paste0(fragment, ".tbi"),
                       paste0(fragment, ".csi"))
  index.file <- index.filepaths[file.exists(index.filepaths)]
  if (length(x = index.file) == 0 & !is.remote) {
    stop("Fragment file is not indexed.")
  } else if(length(x = index.file) == 0) {
    if (verbose) {
      message("Fragment file is on a remote server")
    }
    index.file = paste0(fragment, ".tbi")
  } else if (length(x = index.file) == 2) {
    if (verbose) {
      message("TBI and CSI index both present, using TBI index")
    }
    index.file <- index.file[1]
  }
  return(index.file)
}


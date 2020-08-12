#' @include generics.R
#'
NULL

#' Count fragments
#'
#' Count total fragments per cell barcode present in a fragment file.
#'
#' @param fragments Path to a fragment file
#' @param cells Cells to include. If NULL, include all cells
#' @param max_lines Maximum number of lines to read from the fragment file. If
#' NULL, read all lines in the file.
#' @param verbose Display messages
#'
#' @rdname CountFragments
#' @export
#' @concept fragments
#' @return Returns a data.frame
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' counts <- CountFragments(fragments = fpath)
CountFragments <- function(
  fragments,
  cells = NULL,
  max_lines = NULL,
  verbose = TRUE
) {
  fragments <- normalizePath(path = fragments, mustWork = TRUE)
  max_lines <- SetIfNull(x = max_lines, y = 0)
  verbose = as.logical(x = verbose)
  counts <- groupCommand(
    fragments = fragments,
    some_whitelist_cells = cells,
    max_lines = max_lines,
    verbose = verbose
  )
  return(counts)
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

#' Create a Fragment object
#'
#' Create a \code{Fragment} object to store fragment file information.
#' This object stores a 32-bit MD5 hash of the fragment file and the fragment
#' file index so that any changes to the files on-disk can be detected. A check
#' is also performed to ensure that the expected cells are present in the
#' fragment file.
#'
#' @param path A path to the fragment file. The file should contain a tabix
#' index in the same directory.
#' @param cells A named character vector containing cell barcodes contained in
#' the fragment file. This does not need to be all cells in the fragment file,
#' but there should be no cells in the vector that are not present in the
#' fragment file. A search of the file will be performed until at least one
#' fragment from each cell is found. If NULL, don't check for expected cells.
#'
#' Each element of the vector should be a cell barcode that appears in the
#' fragment file, and the name of each element should be the corresponding cell
#' name in the object.
#' @param validate.fragments Check that expected cells are present in the
#' fragment file.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{ValidateCells}
#'
#' @importFrom tools md5sum file_ext
#' @export
#' @concept fragments
#'
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' cells <- colnames(x = atac_small)
#' names(x = cells) <- paste0("test_", cells)
#' frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
CreateFragmentObject <- function(
  path,
  cells = NULL,
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
) {
  # check that file exists and is indexed
  # don't check if supplying remote file
  is.remote <- grepl(pattern = "^http|^ftp", x = path)
  if (!file.exists(path) & !is.remote) {
    stop("Fragment file does not exist.")
  }
  index.file <- paste0(path, ".tbi")
  if (!file.exists(index.file) & !is.remote) {
    stop("Fragment file is not indexed.")
  }
  if (!is.null(x = cells)) {
    if (is.null(names(x = cells))) {
      # assume cells are as they appear in the assay
      names(x = cells) <- cells
    }
  }
  # compute hash of the file and index
  if (verbose) {
    message("Computing hash")
  }
  if (!is.remote) {
    path <- normalizePath(path = path, mustWork = TRUE)
  }
  # will be NA if file remote
  hashes <- md5sum(files = c(path, index.file))
  # create object
  frags <- new(
    Class = "Fragment",
    path = path,
    hash = unname(obj = hashes),
    cells = cells
  )
  # validate cells
  if (!is.null(x = cells) & validate.fragments) {
    if (ValidateCells(object = frags, verbose = verbose, ...)) {
      return(frags)
    } else {
      stop("Not all cells requested could be found in the fragment file.")
    }
  } else {
    return(frags)
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
#' @importFrom data.table fread
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
  filepath <- GetFragmentData(object = object, slot = "path")
  filepath <- normalizePath(path = filepath, mustWork = TRUE)
  is.remote <- grepl(pattern = "^http|^ftp", x = filepath)
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
  path <- GetFragmentData(object = object, slot = "path")
  index.file <- paste0(path, ".tbi")
  is.remote <- grepl(pattern = "^http|^ftp", x = path)
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
#' @method Cells Fragment
#' @importFrom Seurat Cells
#' @export
Cells.Fragment <- function(x, ...) {
  cells <- slot(object = x, name = "cells")
  return(names(x = cells))
}

# Re-export Seurat generic
#' @export
Seurat::Cells

#' @param value A vector of cell names to store in the \code{\link{Fragment}}
#' object
#' @rdname Cells
#' @export
#' @concept fragments
#' @method Cells<- Fragment
"Cells<-.Fragment" <- function(x, ..., value) {
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

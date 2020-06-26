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
#' @importFrom tools md5sum file_ext
#' @export
#' @concept fragments
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
  # file must end in gz otherwise data.table::fread fails
  if (file_ext(x = path) != "gz") {
    stop("File must end in .gz")
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
#' @param chunksize Number of fragments to read at each iteration.
#' @param tolerance Fraction of input cells that can be unseen before returning
#' TRUE. For example, \code{tolerance = 0.01} will return TRUE when 99% of cells
#' have observed fragments in the file. This can be useful if there are cells
#' present that have much fewer total counts, and would require extensive
#' searching before a fragment from those cells are found.
#' @param max.iter Maximum number of chunks to read in without finding the
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
  chunksize = 1e7,
  tolerance = 0.2,
  max.iter = 5,
  verbose = TRUE
) {
  cell_barcodes <- GetFragmentData(object = object, slot = "cells")
  cells <- SetIfNull(x = cells, y = cell_barcodes)
  if (is.null(x = cells)) {
    warning("No cells stored in object")
    return(TRUE)
  }
  filepath <- GetFragmentData(object = object, slot = "path")
  is.remote <- grepl(pattern = "^http|^ftp", x = filepath)
  # if remote, return TRUE
  if (is.remote) {
    return(TRUE)
  }
  x <- 0
  min.cells <- length(x = cells) - round(x = tolerance * length(x = cells))
  while (TRUE) {
    if (verbose) {
      message("Reading ", chunksize, " fragments")
    }
    chunk <- readchunk(filepath = filepath, x = x, chunksize = chunksize)
    if (isFALSE(x = chunk) | (nrow(x = chunk) == 0)) {
      return(FALSE)
    }
    cells <- setdiff(x = cells, y = unique(x = chunk$cell))
    if (length(x = cells) <= min.cells) {
      return(TRUE)
    }
    if ((!is.null(x = max.iter)) & (x >= max.iter)) {
      return(FALSE)
    }
    x <- x + 1
  }
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

# Read chunk of a file, return FALSE if reading gives an error,
# otherwise return the chunk of file
#' @importFrom future nbrOfWorkers
readchunk <- function(filepath, x, chunksize) {
  tryCatch(
    expr = fread(
      file = filepath,
      skip = x * chunksize,
      nrows = chunksize,
      col.names = c("chr", "start", "end", "cell", "count"),
      nThread = nbrOfWorkers()
    ),
    error = function(x) return(FALSE)
  )
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
#' @export
#' @concept fragments
#' @method Cells Fragment
#' @importFrom Seurat Cells
Cells.Fragment <- function(x, ...) {
  cells <- slot(object = x, name = "cells")
  return(names(x = cells))
}

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

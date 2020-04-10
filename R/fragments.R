#' @include generics.R
#'
NULL

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
#' @param cells A character vector containing cell barcodes contained in the
#' fragment file. This does not need to be all cells in the fragment file,
#' but there should be no cells in the vector that are not present in the
#' fragment file. A search of the file will be performed until at least one
#' fragment from each cell is found.
#' @param prefix A prefix to attach to cell barcodes in the file. This will be
#' automatically added to cell barcodes that are returned by functions that
#' use the \code{Fragment} class object. If NULL, don't append a prefix.
#' @param suffix The same as the \code{prefix} argument, but for a suffix.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{ValidateCells}
#' @importFrom tools md5sum file_ext
CreateFragmentObject <- function(
  path,
  cells,
  prefix = NULL,
  suffix = NULL,
  verbose = TRUE,
  ...
) {
  # check that file exists and is indexed
  index.file <- paste0(path, ".tbi")
  if (!file.exists(path)) {
    stop("Fragment file does not exist.")
  }
  if (!file.exists(index.file)) {
    stop("Fragment file is not indexed.")
  }
  # file must end in gz otherwise data.table::fread fails
  if (file_ext(x = path) != "gz") {
    stop("File must end in .gz")
  }
  # compute hash of the file and index
  if (verbose) {
    message("Computing hash")
  }
  hashes <- md5sum(files = c(path, index.file))
  # create object
  frags <- new(
    Class = "Fragment",
    path = path,
    hash = unname(obj = hashes),
    cells = cells,
    prefix = SetIfNull(x = prefix, y = ""),
    suffix = SetIfNull(x = suffix, y = "")
  )
  # validate cells
  if (ValidateCells(object = frags, verbose = verbose, ...)) {
    return(frags)
  } else {
    stop("Not all cells requested could be found in the fragment file.")
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
#' @importFrom data.table fread
ValidateCells <- function(
  object,
  cells = NULL,
  chunksize = 500000,
  tolerance = 0.05,
  max.iter = 8,
  verbose = TRUE
) {
  cells <- SetIfNull(x = cells, y = Cells(x = object))
  filepath <- GetFragmentData(object = object, slot = "path")
  x <- 0
  min.cells <- round(x = tolerance * length(x = cells))
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
    if (!is.null(x = max.iter) & max.iter >= x) {
      return(FALSE)
    }
    x <- x + 1
  }
}

#' Validate hashes for Fragment object
#'
#' @param object A Fragment object
#' @param verbose Display messages
#' @export
#' @importFrom tools md5sum
ValidateHash <- function(object, verbose = TRUE) {
  path <- GetFragmentData(object = object, slot = "path")
  index.file <- paste0(path, ".tbi")
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
ValidateFragments <- function(
  object,
  ...
) {
  valid.cells <- ValidateCells(object = object, verbose = verbose, ...)
  valid.hash <- ValidateHash(object = object, verbose = verbose, ...)
  return(valid.cells & valid.hash)
}

# Read chunk of a file, return FALSE if reading gives an error,
# otherwise return the chunk of file
readchunk <- function(filepath, x, chunksize) {
  tryCatch(
    expr = fread(
      file = filepath,
      skip = x * chunksize,
      nrows = chunksize,
      col.names = c("chr", "start", "end", "cell", "count")
    ),
    error = function(x) return(FALSE)
  )
}

#' @rdname Cells
#' @export
#' @method Cells Fragment
#' @importFrom Seurat Cells
Cells.Fragment <- function(x) {
  return(slot(object = x, name = "cells"))
}

globalVariables(names = c("chr", "start"), package = "Signac")
#' FilterFragments
#'
#' Remove cells from a fragments file that are not present in a given list of
#' cells. Note that this reads the whole fragments file into memory, so may
#' require a lot of memory depending on the size of the fragments file.
#'
#' @param fragment.path Path to a tabix-indexed fragments file
#' @param cells A vector of cells to retain
#' @param output.path Name and path for output tabix file. A tabix index file
#' will also be created in the same location, with the .tbi file extension.
#' @param assume.sorted Assume sorted input and don't sort the filtered file.
#' Can save a lot of time, but indexing will fail if assumption is wrong.
#' @param compress Compress filtered fragments using bgzip (default TRUE)
#' @param index Index the filtered tabix file (default TRUE)
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[data.table]{fread}}
#'
#' @importFrom data.table fread fwrite
#' @importFrom Rsamtools indexTabix bgzip
#' @export
#' @return None
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' output.path = file.path(tempdir(), "filtered.tsv")
#'
#' FilterFragments(
#'   fragment.path = fpath,
#'   cells = colnames(atac_small),
#'   output.path = output.path
#' )
#' }
FilterFragments <- function(
  fragment.path,
  cells,
  output.path,
  assume.sorted = FALSE,
  compress = TRUE,
  index = TRUE,
  verbose = TRUE,
  ...
) {
  if (verbose) {
    message("Retaining ", length(x = cells), " cells")
    message("Reading fragments")
  }
  reads <- fread(
    file = fragment.path,
    col.names = c("chr", "start", "end", "cell", "count"),
    showProgress = verbose,
    ...
  )
  reads <- reads[reads$cell %in% cells, ]
  if (!assume.sorted) {
    if (verbose) {
      message("Sorting fragments")
    }
    reads <- reads[with(data = reads, expr = order(chr, start)), ]
  }
  if (verbose) {
    message("Writing output")
  }
  fwrite(
    x = reads,
    file = output.path,
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
  rm(reads)
  invisible(x = gc())
  if (compress) {
    if (verbose) {
      message("Compressing output")
    }
    outf <- bgzip(file = output.path, dest = paste0(output.path, ".gz"))
    if (file.exists(outf)) {
      file.remove(output.path)
    }
    if (index) {
      if (verbose) {
        message("Building index")
      }
      index.file <- indexTabix(
        file = paste0(outf), format = "bed", zeroBased = TRUE
      )
    }
  }
}


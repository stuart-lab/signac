#' Read MGATK output
#'
#' Read output files from MGATK (\url{https://github.com/caleblareau/mgatk}) and create an
#' Assay object.
#'
#' @param dir Path to directory containing MGATK output files
#' @param verbose Display messages
#'
#' @importFrom Matrix sparseMatrix
#'
#' @return Returns a list containing a sparse matrix and a dataframe.
#' The sparse matrix contains read counts for each base at each position
#' and strand, and the dataframe contains the total depth for each cell.
#'
#' @export
ReadMGATK <- function(dir, verbose = TRUE){
  if (!dir.exists(paths = dir)) {
    stop("Directory not found")
  }
  a.path <- list.files(path = dir, pattern = "*.A.txt.gz", full.names = TRUE)
  c.path <- list.files(path = dir, pattern = "*.C.txt.gz", full.names = TRUE)
  t.path <- list.files(path = dir, pattern = "*.T.txt.gz", full.names = TRUE)
  g.path <- list.files(path = dir, pattern = "*.G.txt.gz", full.names = TRUE)
  refallele.path <- list.files(path = dir, pattern = "chrM_refAllele.txt", full.names = TRUE)
  depthfile.path <- list.files(path = dir, pattern = "*.depthTable.txt", full.names = TRUE)
  cellbarcode.path <- list.files(path = dir, pattern = "barcodes.tsv", full.names = TRUE)

  if (verbose) {
    message("Reading allele counts")
  }
  column.names <- c("pos", "cellbarcode", "plus", "minus")
  a.counts <- read.table(file = a.path, sep = ',', header = FALSE, stringsAsFactors = FALSE, col.names = column.names)
  c.counts <- read.table(file = c.path, sep = ',', header = FALSE, stringsAsFactors = FALSE, col.names = column.names)
  t.counts <- read.table(file = t.path, sep = ',', header = FALSE, stringsAsFactors = FALSE, col.names = column.names)
  g.counts <- read.table(file = g.path, sep = ',', header = FALSE, stringsAsFactors = FALSE, col.names = column.names)

  if (verbose) {
    message("Reading metadata")
  }
  refallele <- read.table(file = refallele.path, header = FALSE, stringsAsFactors = FALSE, col.names = c("pos", "ref"))
  depth <- read.table(file = depthfile.path, header = FALSE, stringsAsFactors = FALSE, col.names = c("cellbarcode", "depth"))
  cellbarcodes <- unique(x = readLines(con = cellbarcode.path))
  cb.lookup <- seq_along(along.with = cellbarcodes)
  names(cb.lookup) <- cellbarcodes

  if (verbose) {
    message("Building matrices")
  }
  a.mat <- SparseMatrixFromBaseCounts(basecounts = a.counts, cells = cb.lookup, dna.base = "A")
  c.mat <- SparseMatrixFromBaseCounts(basecounts = c.counts, cells = cb.lookup, dna.base = "C")
  t.mat <- SparseMatrixFromBaseCounts(basecounts = t.counts, cells = cb.lookup, dna.base = "T")
  g.mat <- SparseMatrixFromBaseCounts(basecounts = g.counts, cells = cb.lookup, dna.base = "G")

  counts <- rbind(a.mat[[1]], c.mat[[1]], t.mat[[1]], g.mat[[1]],
                  a.mat[[2]], c.mat[[2]], t.mat[[2]], g.mat[[2]])

  return(list("counts" = counts, "depth" = depth))
}

# Create sparse matrix from base counts
#
# Create a sparse position x cell matrix for
# containing read counts for each base in each
# cell, for + and - strands.
#
# @param basecounts A dataframe containing read counts at each position for each cell
# @param cells A lookup table giving the cell barcode numeric ID
#
# @return Returns a list of two sparse matrices
SparseMatrixFromBaseCounts <- function(basecounts, cells, dna.base) {
  fwd.mat <- sparseMatrix(
    i = basecounts$pos,
    j = cells[basecounts$cellbarcode],
    x = basecounts$plus
  )
  colnames(x = fwd.mat) <- names(x = cells)
  rownames(x = fwd.mat) <- paste(dna.base, 1:nrow(fwd.mat), "fwd", sep = "_")
  rev.mat <- sparseMatrix(
    i = basecounts$pos,
    j = cells[basecounts$cellbarcode],
    x = basecounts$minus
  )
  colnames(x = rev.mat) <- names(x = cells)
  rownames(x = rev.mat) <- paste(dna.base, 1:nrow(rev.mat), "rev", sep = "_")
  return(list(fwd.mat, rev.mat))
}

#' Call mitochondrial variants
#'
#' Identify mitochondrial variants present in single cells.
#'
#' @param object An Assay object
#' @param verbose Display messages
#'
#' @return Returns a dataframe
#' @export
CallVariants <- function(object, verbose = TRUE, ...) {
  return()
}

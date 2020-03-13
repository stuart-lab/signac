#' Read MGATK output
#'
#' Read output files from MGATK (\url{https://github.com/caleblareau/mgatk}).
#'
#' @param dir Path to directory containing MGATK output files
#' @param verbose Display messages
#'
#' @return Returns a list containing a sparse matrix (counts) and two dataframes (depth and refallele).
#' The sparse matrix contains read counts for each base at each position
#' and strand.
#' The depth dataframe contains the total depth for each cell.
#' The refallele dataframe contains the reference genome allele at each position.
#'
#' @export
#' @examples
#' \dontrun{
#' data.dir <- "path/to/data/directory"
#' mgatk <- ReadMGATK(dir = data.dir)
#' }
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
  depth <- read.table(file = depthfile.path, header = FALSE, stringsAsFactors = FALSE, col.names = c("cellbarcode", "mito.depth"), row.names = 1)
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

  return(list("counts" = counts, "depth" = depth, "refallele" = refallele))
}

# Create sparse matrix from base counts
#
# Create a sparse position x cell matrix for
# containing read counts for each base in each
# cell, for + and - strands.
#
# @param basecounts A dataframe containing read counts at each position for each cell
# @param cells A lookup table giving the cell barcode numeric ID
#' @importFrom Matrix sparseMatrix
#
# @return Returns a list of two sparse matrices
SparseMatrixFromBaseCounts <- function(basecounts, cells, dna.base) {
  fwd.mat <- sparseMatrix(
    i = basecounts$pos,
    j = cells[basecounts$cellbarcode],
    x = basecounts$plus
  )
  colnames(x = fwd.mat) <- names(x = cells)
  rownames(x = fwd.mat) <- paste(dna.base, 1:nrow(fwd.mat), "fwd", sep = "-")
  rev.mat <- sparseMatrix(
    i = basecounts$pos,
    j = cells[basecounts$cellbarcode],
    x = basecounts$minus
  )
  colnames(x = rev.mat) <- names(x = cells)
  rownames(x = rev.mat) <- paste(dna.base, 1:nrow(rev.mat), "rev", sep = "-")
  return(list(fwd.mat, rev.mat))
}

#' Call mitochondrial variants
#'
#' Identify mitochondrial variants present in single cells.
#'
#' @param object An Assay object
#' @param refallele A dataframe containing reference alleles for the
#' mitochondrial genome.
#' @param stabilize.variance Stabilize variance
#' @param low.coverage.threshold Low coverage threshold
#' @param verbose Display messages
#'
#' @return Returns a dataframe
#' @export
CallVariants <- function(
  object,
  refallele,
  stabilize.variance = TRUE,
  low.coverage.threshold = 10,
  verbose = TRUE
) {
  # determine key coverage statistics
  coverages <- ComputeTotalCoverage(object = object, verbose = verbose)
  ref_allele <- toupper(x = refallele$ref)

  # TODO
  return()
}


####### Not exported

# Compute total mitochondrial coverage
#
# Sum the counts for each base and each strand to
# find the total coverage at each mitochondrial base
# for each cell
#
# Assumes that the matrix is stacked by row, such that
# the first n rows correspond to the mitochondrial base
# positions in order, for the first base and first strand,
# followed by n:2n rows for the next base, etc.
#
# @param object A matrix containing coverages for each
# base and each strand
# @param verbose Display messages
# @return Returns a matrix
ComputeTotalCoverage <- function(object, verbose = TRUE) {
  if (verbose) {
    message("Computing total coverage per base")
  }
  rowstep <- nrow(x = object) / 8
  mat.list <- list()
  for (i in seq_len(length.out = 8)) {
    mat.list[[i]] <- object[(rowstep * (i-1) + 1):(rowstep*i), ]
  }
  coverage <- Reduce(f = `+`, x = mat.list)
  coverage <- as.matrix(x = coverage)
  rownames(x = coverage) <- seq_along(along.with = rownames(x = coverage))
  return(coverage)
}

# Process mutation for one DNA base
#
# @param letter DNA base to process
# @param verbose Display messages
ProcessLetter <- function(letter, verbose = TRUE) {
  # TODO
  return()
}

# Extract mutation matrix
#
# Subset the mitochondrial count matrix by nucleotide
# and strand.
#
# @param object A matrix containing counts for each
# mitochondrial base and strand
# @param letter Which DNA base to use (A, T, C, G)
# @param strand Which strand to use (fwd, rev)
# @return Returns a sparse matrix
GetMutationMatrix <- function(object, letter, strand) {
  keep.rows <- paste(letter, seq_len(length.out = nrow(x = object) / 8), strand, sep = '-')
  return(object[keep.rows, ])
}


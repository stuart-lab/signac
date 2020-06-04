#' Read MGATK output
#'
#' Read output files from MGATK (\url{https://github.com/caleblareau/mgatk}).
#'
#' @param dir Path to directory containing MGATK output files
#' @param verbose Display messages
#'
#' @return Returns a list containing a sparse matrix (counts) and two dataframes
#' (depth and refallele).
#'
#' The sparse matrix contains read counts for each base at each position
#' and strand.
#'
#' The depth dataframe contains the total depth for each cell.
#' The refallele dataframe contains the reference genome allele at each
#' position.
#'
#' @export
#' @concept mito
#' @examples
#' \dontrun{
#' data.dir <- system.file("extdata", "test_mgatk", package="Signac")
#' mgatk <- ReadMGATK(dir = data.dir)
#' }
ReadMGATK <- function(dir, verbose = TRUE) {
  if (!dir.exists(paths = dir)) {
    stop("Directory not found")
  }
  a.path <- list.files(path = dir, pattern = "*.A.txt.gz", full.names = TRUE)
  c.path <- list.files(path = dir, pattern = "*.C.txt.gz", full.names = TRUE)
  t.path <- list.files(path = dir, pattern = "*.T.txt.gz", full.names = TRUE)
  g.path <- list.files(path = dir, pattern = "*.G.txt.gz", full.names = TRUE)

  refallele.path <- list.files(
    path = dir,
    pattern = "chrM_refAllele.txt*",
    full.names = TRUE
  )

  # The depth file lists all barcodes that were genotyped
  depthfile.path <- list.files(
    path = dir,
    pattern = "*.depthTable.txt",
    full.names = TRUE
  )

  if (verbose) {
    message("Reading allele counts")
  }
  column.names <- c("pos", "cellbarcode", "plus", "minus")
  a.counts <- read.table(
    file = a.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )
  c.counts <- read.table(
    file = c.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )
  t.counts <- read.table(
    file = t.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )
  g.counts <- read.table(
    file = g.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )

  if (verbose) {
    message("Reading metadata")
  }
  refallele <- read.table(
    file = refallele.path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("pos", "ref")
  )
  refallele$ref <- toupper(x = refallele$ref)
  depth <- read.table(
    file = depthfile.path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("cellbarcode", "mito.depth"),
    row.names = 1
  )
  cellbarcodes <- unique(x = rownames(depth))
  cb.lookup <- seq_along(along.with = cellbarcodes)
  names(cb.lookup) <- cellbarcodes

  if (verbose) {
    message("Building matrices")
  }

  maxpos <- dim(refallele)[1]
  a.mat <- SparseMatrixFromBaseCounts(
    basecounts = a.counts, cells = cb.lookup, dna.base = "A", maxpos = maxpos
  )
  c.mat <- SparseMatrixFromBaseCounts(
    basecounts = c.counts, cells = cb.lookup, dna.base = "C", maxpos = maxpos
  )
  t.mat <- SparseMatrixFromBaseCounts(
    basecounts = t.counts, cells = cb.lookup, dna.base = "T", maxpos = maxpos
  )
  g.mat <- SparseMatrixFromBaseCounts(
    basecounts = g.counts, cells = cb.lookup, dna.base = "G", maxpos = maxpos
  )

  counts <- rbind(a.mat[[1]], c.mat[[1]], t.mat[[1]], g.mat[[1]],
                  a.mat[[2]], c.mat[[2]], t.mat[[2]], g.mat[[2]])

  return(list("counts" = counts, "depth" = depth, "refallele" = refallele))
}

#' @param refallele A dataframe containing reference alleles for the
#' mitochondrial genome.
#' @param stabilize.variance Stabilize variance
#' @param low.coverage.threshold Low coverage threshold
#' @param verbose Display messages
#'
#' @return Returns a dataframe
#' @export
#' @concept mito
#' @rdname IdentifyVariants
#' @examples
#' \dontrun{
#' data.dir <- "path/to/data/directory"
#' mgatk <- ReadMGATK(dir = data.dir)
#' variant.df <- IdentifyVariants(
#'   object = mgatk$counts,
#'   refallele = mgatk$refallele
#' )
#' }
IdentifyVariants.default <- function(
  object,
  refallele,
  stabilize_variance = TRUE,
  low_coverage_threshold = 10,
  verbose = TRUE,
  ...
) {
  coverages <- ComputeTotalCoverage(object = object, verbose = verbose)
  a.df <- ProcessLetter(
    object = object,
    letter = "A",
    coverage = coverages,
    ref_alleles = refallele,
    stabilize_variance = stabilize_variance,
    low_coverage_threshold = low_coverage_threshold,
    verbose = verbose
  )
  t.df <- ProcessLetter(
    object = object,
    letter = "T",
    coverage = coverages,
    ref_alleles = refallele,
    stabilize_variance = stabilize_variance,
    low_coverage_threshold = low_coverage_threshold,
    verbose = verbose
  )
  c.df <- ProcessLetter(
    object = object,
    letter = "C",
    coverage = coverages,
    ref_alleles = refallele,
    stabilize_variance = stabilize_variance,
    low_coverage_threshold = low_coverage_threshold,
    verbose = verbose
  )
  g.df <- ProcessLetter(
    object = object,
    letter = "G",
    coverage = coverages,
    ref_alleles = refallele,
    stabilize_variance = stabilize_variance,
    low_coverage_threshold = low_coverage_threshold,
    verbose = verbose
  )
  return(rbind(a.df, t.df, c.df, g.df))
}

#' @importFrom Seurat GetAssayData
#' @rdname IdentifyVariants
#' @method IdentifyVariants Assay
#' @concept mito
#' @export
IdentifyVariants.Assay <- function(
  object,
  refallele,
  ...
) {
  counts <- GetAssayData(object = object, slot = 'counts')
  df <- IdentifyVariants(object = counts, refallele = refallele, ...)
  return(df)
}

#' @importFrom Seurat GetAssay DefaultAssay
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @rdname IdentifyVariants
#' @concept mito
#' @method IdentifyVariants Seurat
#' @export
IdentifyVariants.Seurat <- function(
  object,
  refallele,
  assay = NULL,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  assay.obj <- GetAssay(object = object, assay = assay)
  df <- IdentifyVariants(object = assay.obj, refallele = refallele, ...)
  return(df)
}

####### Not exported

# Create sparse matrix from base counts
#
# Create a sparse position x cell matrix for
# containing read counts for each base in each
# cell, for + and - strands.
#
# @param basecounts A dataframe containing read counts at each position for each
# cell
# @param cells A lookup table giving the cell barcode numeric ID
# @param dna.base Used to specify the alternate allele
# @param maxpos specifies the end of the mtDNA genome (otherwise, the mat will be of wrong dim)
#' @importFrom Matrix sparseMatrix
#
# @return Returns a list of two sparse matrices
SparseMatrixFromBaseCounts <- function(basecounts, cells, dna.base, maxpos) {

  # Vector addition guarantee correct dimension
  fwd.mat <- sparseMatrix(
    i = c(basecounts$pos,maxpos),
    j = c(cells[basecounts$cellbarcode],1),
    x = c(basecounts$plus,0)
  )
  colnames(x = fwd.mat) <- names(x = cells)
  rownames(x = fwd.mat) <- paste(
    dna.base,
    seq_len(length.out = nrow(fwd.mat)),
    "fwd",
    sep = "-"
  )
  rev.mat <- sparseMatrix(
    i = c(basecounts$pos,maxpos),
    j = c(cells[basecounts$cellbarcode],1),
    x = c(basecounts$minus,0)
  )
  colnames(x = rev.mat) <- names(x = cells)
  rownames(x = rev.mat) <- paste(
    dna.base,
    seq_len(length.out = nrow(rev.mat)),
    "rev",
    sep = "-"
  )
  return(list(fwd.mat, rev.mat))
}

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
    mat.list[[i]] <- object[(rowstep * (i - 1) + 1):(rowstep * i), ]
  }
  coverage <- Reduce(f = `+`, x = mat.list)
  coverage <- as.matrix(x = coverage)
  rownames(x = coverage) <- seq_along(along.with = rownames(x = coverage))
  return(coverage)
}

# Process mutation for one DNA base
#
# @param object A matrix containing nucleotide counts for each base and strand
# of the genome. This should be a stacked matrix, such that the number of
# rows = 8 * the length of the genome.
# @param letter DNA base to process
# @param ref_alleles A dataframe containing reference genome alleles and
# position
# @param coverage Total coverage for all bases and strands
# @param verbose Display messages
#' @importFrom Matrix summary rowMeans rowSums
#' @import data.table
ProcessLetter <- function(
  object,
  letter,
  ref_alleles,
  coverage,
  stabilize_variance = TRUE,
  low_coverage_threshold = 10,
  verbose = TRUE
) {
  if (verbose) {
    message("Processing ", letter)
  }
  boo <- ref_alleles$ref != letter & ref_alleles$ref != "N"
  cov <- coverage[boo, ]

  variant_name <- paste0(
    as.character(ref_alleles$pos),
    ref_alleles$ref,
    ">",
    letter
  )[boo]

  nucleotide <- paste0(
    ref_alleles$ref,
    ">",
    letter
  )[boo]

  position_filt <- ref_alleles$pos[boo]

  # get forward and reverse counts
  fwd.counts <- GetMutationMatrix(
    object = object,
    letter = letter,
    strand = "fwd"
  )[boo, ]

  rev.counts <- GetMutationMatrix(
    object = object,
    letter = letter,
    strand = "rev"
  )[boo, ]

  # convert to row, column, value
  fwd.ijx <- summary(fwd.counts)
  rev.ijx <- summary(rev.counts)

  # get bulk coverage stats
  bulk <- (rowSums(fwd.counts + rev.counts) / rowSums(cov))
  # replace NA or NaN with 0
  bulk[is.na(bulk)] <- 0
  bulk[is.nan(bulk)] <- 0

  # find correlation between strands for cells with >0 counts on either strand
  # group by variant (row) and find correlation between strand depth
  both.strand <- data.table(cbind(fwd.ijx, rev.ijx$x))
  both.strand$i <- variant_name[both.strand$i]
  colnames(both.strand) <- c("variant", "cell_idx", "forward", "reverse")

  cor_dt <- suppressWarnings(expr = both.strand[, .(cor = cor(
    x = forward, y = reverse, method = "pearson", use = "pairwise.complete")
  ), by = list(variant)])

  # Put in vector for convenience
  cor_vec_val <- cor_dt$cor
  names(cor_vec_val) <- as.character(cor_dt$variant)

  # Compute the single-cell data
  mat <- (fwd.counts + rev.counts) / cov
  rownames(mat) <- variant_name
  # set NAs and NaNs to zero
  mat@x[!is.finite(mat@x)] <- 0

  # Stablize variances by replacing low coverage cells with mean
  if (stabilize_variance) {
    idx_mat <- which(cov < low_coverage_threshold, arr.ind = TRUE)
    idx_mat_mean <- bulk[idx_mat[, 1]]
    ones <- 1 - sparseMatrix(
      i = c(idx_mat[, 1], dim(x = mat)[1]),
      j = c(idx_mat[, 2], dim(x = mat)[2]),
      x = 1
    )
    means_mat <- sparseMatrix(
      i = c(idx_mat[, 1], dim(x = mat)[1]),
      j = c(idx_mat[, 2], dim(x = mat)[2]),
      x = c(idx_mat_mean, 0)
    )
    mmat2 <- mat * ones + means_mat
    variance <- SparseRowVar(x = mmat2)
  } else {
    variance <- SparseRowVar(x = mat)
  }
  detected <- (fwd.counts >= 2) + (rev.counts >= 2)

  # Compute per-mutation summary statistics
  var_summary_df <- data.frame(
    position = position_filt,
    nucleotide = nucleotide,
    variant = variant_name,
    vmr = variance / (bulk + 0.00000000001),
    mean = round(x = bulk, digits = 7),
    variance = round(x = variance, digits = 7),
    n_cells_conf_detected = rowSums(x = detected == 2),
    n_cells_over_5 = rowSums(x = mat >= 0.05),
    n_cells_over_10 = rowSums(x = mat >= 0.10),
    n_cells_over_20 = rowSums(x = mat >= 0.20),
    strand_correlation = cor_vec_val[variant_name],
    mean_coverage = rowMeans(x = cov),
    stringsAsFactors = FALSE,
    row.names = variant_name
  )
  return(var_summary_df)
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
  keep.rows <- paste(
    letter,
    seq_len(length.out = nrow(x = object) / 8),
    strand,
    sep = "-"
  )
  return(object[keep.rows, ])
}

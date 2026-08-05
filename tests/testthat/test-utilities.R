library(GenomicRanges)
library(SeuratObject)
library(Matrix)

suppressWarnings(RNGversion(vstr = "3.5.3"))

# AverageCounts / CellsPerGroup / SortIdents / BinaryIdentMatrix ---------------

test_that("AverageCounts works", {
  expect_equal(
    object = as.vector(x = AverageCounts(object = atac_small)),
    expected = 30.11,
    tolerance = 1 / 1000
  )
})

test_that("AverageCounts works with group.by", {
  Idents(atac_small) <- "cluster"
  groups <- sort(x = unique(x = as.character(x = atac_small$cluster)))
  res <- AverageCounts(atac_small, group.by = "cluster")
  # one positive mean per cluster, named by cluster
  expect_type(res, "double")
  expect_equal(length(x = res), length(x = groups))
  expect_setequal(names(x = res), groups)
  expect_true(all(res > 0))
})

test_that("AverageCounts without group.by uses Idents", {
  Idents(atac_small) <- "cluster"
  # with idents set to cluster, omitting group.by must reproduce the
  # per-cluster means returned when group.by = "cluster" is given explicitly
  res <- AverageCounts(atac_small, group.by = NULL, verbose = FALSE)
  expect_equal(res, AverageCounts(atac_small, group.by = "cluster"))
})

test_that("AverageCountMatrix returns features x groups matrix", {
  Idents(atac_small) <- "cluster"
  res <- Signac:::AverageCountMatrix(
    object = atac_small, group.by = "cluster"
  )
  expect_true(is.matrix(res) || inherits(res, "Matrix"))
  expect_equal(nrow(res), nrow(atac_small[["peaks"]]))
  # one column per group identity present
  expect_equal(ncol(res), length(unique(atac_small$cluster)))
  expect_equal(rownames(res), rownames(atac_small[["peaks"]]))
  # all values are non-negative (averaged from counts)
  expect_true(all(res >= 0))
})

test_that("AverageCountMatrix errors on missing layer", {
  expect_error(
    Signac:::AverageCountMatrix(
      object = atac_small, layer = "nonexistent_layer"
    ),
    regexp = "layer is not present"
  )
})

test_that("CellsPerGroup works", {
  Idents(atac_small) <- "cluster"
  expect_equal(
    object = as.vector(x = CellsPerGroup(object = atac_small)),
    expected = c(55, 45, 0)
  )
})

test_that("CellsPerGroup with explicit group.by", {
  Idents(atac_small) <- "cluster"
  res <- CellsPerGroup(atac_small, group.by = "cluster")
  expect_true(is.numeric(res))
  # cluster has values 1 and 2; with useNA="always" there's also an NA bucket
  expect_equal(sum(res), ncol(atac_small))
  expect_true("1" %in% names(res))
})

test_that("BinaryIdentMatrix encodes one-hot group membership", {
  Idents(atac_small) <- "cluster"
  m <- Signac:::BinaryIdentMatrix(object = atac_small)
  expect_s4_class(m, "CsparseMatrix")
  expect_equal(ncol(m), ncol(atac_small))
  # rows are unique groups
  expect_equal(nrow(m), length(unique(Idents(atac_small))))
  # each cell belongs to exactly one group (one-hot)
  expect_equal(as.vector(Matrix::colSums(m)), rep(1, ncol(atac_small)))
})

test_that("SortIdents works", {
  set.seed(1)
  atac_small$test <- sample(1:10, ncol(atac_small), replace = TRUE)
  atac_small <- SortIdents(object = atac_small, label = "test")
  expect_equal(
    object = levels(atac_small$test),
    expected = c("6", "1", "8", "2", "10", "5", "9", "7", "3", "4")
  )
  Idents(atac_small) <- sample(1:10, ncol(atac_small), replace = TRUE)
  atac_small <- SortIdents(object = atac_small)
  expect_equal(
    object = levels(Idents(atac_small)),
    expected = c("10", "5", "7", "3", "6", "4", "2", "8", "1", "9")
  )
})

test_that("AccessiblePeaks returns peaks open in the requested identity", {
  Idents(atac_small) <- "cluster"
  pks <- AccessiblePeaks(
    object = atac_small, idents = "1", min.cells = 1
  )
  expect_type(pks, "character")
  expect_true(length(pks) > 0)
  # all returned features exist in the assay
  expect_true(all(pks %in% rownames(atac_small[["peaks"]])))
  # increasing min.cells must return a subset (monotone)
  pks2 <- AccessiblePeaks(
    object = atac_small, idents = "1", min.cells = 10
  )
  expect_true(all(pks2 %in% pks))
  expect_lte(length(pks2), length(pks))
})

test_that("GetGroups works with idents subset", {
  Idents(atac_small) <- "cluster"
  g <- Signac:::GetGroups(atac_small, group.by = "cluster", idents = "1")
  expect_true(all(as.character(g) == "1"))
})

# corSparse / Sparse* utilities ------------------------------------------------

test_that("corSparse computes correlation matrix", {
  X <- as(matrix(rnorm(100), 10, 10), "CsparseMatrix")
  res <- corSparse(X)
  expect_equal(dim(res), c(10, 10))
  # must match the dense Pearson correlation
  expect_equal(
    unname(as.matrix(res)), unname(cor(as.matrix(X))), tolerance = 1e-6
  )
})

test_that("corSparse two-matrix and covariance", {
  X <- as(matrix(rnorm(100), 10, 10), "CsparseMatrix")
  Y <- as(matrix(rnorm(100), 10, 10), "CsparseMatrix")
  res <- corSparse(X, Y)
  expect_equal(
    unname(as.matrix(res)),
    unname(cor(as.matrix(X), as.matrix(Y))), tolerance = 1e-6
  )
  res_cov <- corSparse(X, cov = TRUE)
  expect_equal(
    unname(as.matrix(res_cov)), unname(cov(as.matrix(X))), tolerance = 1e-6
  )
})

test_that("corSparse errors on incompatible dimensions", {
  X <- as(matrix(rnorm(100), 10, 10), "CsparseMatrix")
  Y <- as(matrix(rnorm(50), 5, 10), "CsparseMatrix")
  expect_error(corSparse(X, Y), regexp = "same number of rows")
})

test_that("SparseRowVar matches base::apply variance for sparse input", {
  m <- as(matrix(1:20, 4, 5), "CsparseMatrix")
  expect_equal(Signac:::SparseRowVar(m), apply(as.matrix(m), 1, var))
  expect_equal(Signac:::SparseColVar(m), apply(as.matrix(m), 2, var))
})

test_that("SparseSpearmanCor matches stats::cor(method='spearman') on dense", {
  set.seed(1)
  m <- matrix(rpois(100, 2), 10, 10)
  m_sparse <- as(m, "CsparseMatrix")
  res <- Signac:::SparseSpearmanCor(m_sparse)
  expect_equal(dim(res), c(10, 10))
  # must match base R's Spearman correlation
  expect_equal(
    unname(as.matrix(res)),
    unname(cor(m, method = "spearman")), tolerance = 1e-6
  )
})

# SubsetMatrix / AddMissing ----------------------------------------------------

test_that("SubsetMatrix filters by min.rows/min.cols thresholds", {
  # Construct a matrix with known per-row/col nonzero counts:
  #   row a: 4 nonzero, b: 2, c: 1, d/e: 0
  #   col A: 3, B: 2, C: 1, D: 1, E: 0
  mat <- matrix(0, 5, 5, dimnames = list(letters[1:5], LETTERS[1:5]))
  mat[1, 1:4] <- 1
  mat[2, 1:2] <- 1
  mat[3, 1]   <- 1
  # min.rows is a minimum, so rows with rowcount >= 2 are kept ("a" and "b"),
  # and likewise cols "A" and "B"
  out <- SubsetMatrix(mat, min.rows = 2, min.cols = 2, max.row.val = 10)
  expect_equal(rownames(out), c("a", "b"))
  expect_equal(colnames(out), c("A", "B"))
  # min.rows = 3 keeps only the row and column with 3 or more nonzero entries
  out3 <- SubsetMatrix(mat, min.rows = 3, min.cols = 3, max.row.val = 10)
  expect_equal(out3, 1)
  # max.col.val excludes columns whose max equals or exceeds the cutoff
  out2 <- SubsetMatrix(mat, max.row.val = 10, max.col.val = 0.5)
  # No column has max < 0.5 since each retained column has at least one 1
  expect_equal(length(out2), 0)
})

test_that("AddMissing adds zero columns/rows", {
  X <- sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(1, 2),
    dims = c(3, 3),
    dimnames = list(c("a", "b", "c"), c("x", "y", "z"))
  )
  out <- Signac:::AddMissing(X, cells = c("x", "y", "z", "w"))
  expect_equal(ncol(out), 4)
  out2 <- Signac:::AddMissing(X, features = c("a", "b", "c", "d"))
  expect_equal(nrow(out2), 4)
})

# ChunkGRanges / Extend / FindRegion / Lookup ---------------------------------

test_that("ChunkGRanges works", {
  granges <- GRanges(
    seqnames = c("chr1"),
    ranges = IRanges(start = seq(1, 10), end = seq(1, 10) + 1)
  )
  split_ranges <- ChunkGRanges(granges = granges, nchunk = 3)
  correct_split <- list(
    GRanges(seqnames = "chr1", ranges = IRanges(start = 1:3, end = (1:3) + 1)),
    GRanges(seqnames = "chr1", ranges = IRanges(start = 4:6, end = (4:6) + 1)),
    GRanges(seqnames = "chr1", ranges = IRanges(start = 7:10, end = (7:10) + 1))
  )
  expect_equal(object = split_ranges, expected = correct_split)
})

test_that("Extend works", {
  granges <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(start = c(1000, 1200), end = c(10000, 312100)),
    strand = c("+", "-")
  )
  correct_extended <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(start = c(500, 200), end = c(11000, 312600)),
    strand = c("+", "-")
  )
  extended_granges <- Extend(x = granges, upstream = 500, downstream = 1000)
  expect_equal(object = extended_granges, expected = correct_extended)
})

test_that("Extend respects from.midpoint", {
  gr <- GRanges("chr1", IRanges(1000, 2000), strand = "+")
  e <- Extend(gr, upstream = 100, downstream = 100, from.midpoint = TRUE)
  expect_equal(width(e), 201)
  e2 <- Extend(gr, upstream = 100, downstream = 100)
  expect_equal(width(e2), 1201)
})

test_that("Extend on minus strand swaps upstream/downstream", {
  gr <- GRanges("chr1", IRanges(1000, 2000), strand = "-")
  e <- Extend(gr, upstream = 200, downstream = 100)
  expect_equal(start(e), 900)
  expect_equal(end(e), 2200)
})

test_that("Extend on GRanges preserves seqnames", {
  gr <- GRanges("chr1", IRanges(c(1000, 2000), c(1500, 2500)))
  res <- Extend(gr, upstream = 100, downstream = 100)
  expect_equal(as.character(seqnames(res)), c("chr1", "chr1"))
})

test_that("LookupGeneCoords finds known gene", {
  ann <- Annotation(atac_small)
  gene <- ann$gene_name[1]
  gr <- LookupGeneCoords(object = atac_small, gene = gene)
  expect_s4_class(gr, "GRanges")
})

test_that("LookupGeneCoords returns NULL for unknown gene", {
  gr <- LookupGeneCoords(object = atac_small, gene = "ZZZZNotAGene")
  expect_null(gr)
})

test_that("FindRegion handles strings and GRanges", {
  r <- Signac:::FindRegion(atac_small, region = "chr1:1000-2000")
  expect_s4_class(r, "GRanges")
  r2 <- Signac:::FindRegion(
    atac_small, region = GRanges("chr1", IRanges(1000, 2000))
  )
  expect_s4_class(r2, "GRanges")
})

test_that("FindRegion errors on unknown gene", {
  expect_error(
    Signac:::FindRegion(atac_small, region = "ZZZZNotAGene"),
    regexp = "Gene not found"
  )
})

test_that("GetGeneRanges collapses to one range per gene", {
  ann <- Annotation(atac_small)
  collapsed <- Signac:::GetGeneRanges(ranges = ann)
  expect_s4_class(collapsed, "GRanges")
  expect_equal(length(collapsed), length(unique(ann$gene_id)))
  expect_false(any(duplicated(collapsed$gene_id)))
})

test_that("GetGeneRanges handles unstranded ranges", {
  gr <- GRanges(
    "chr1", IRanges(c(100, 200), c(150, 250)),
    strand = c("*", "*"),
    tx_id = c("t1", "t1"),
    gene_id = c("g1", "g1"),
    gene_name = c("g1", "g1"),
    gene_biotype = c("protein_coding", "protein_coding"),
    type = c("exon", "exon")
  )
  res <- Signac:::GetGeneRanges(gr)
  expect_s4_class(res, "GRanges")
  expect_equal(as.character(strand(res)), "+")
  expect_equal(start(res), 100)
  expect_equal(end(res), 250)
})

test_that("GetGeneRanges spans all transcripts of a gene", {
  # two transcripts of one gene with distant first exons
  gr <- GRanges(
    "chr1", IRanges(c(1000, 5000, 20000, 24000), c(1200, 5200, 20200, 24200)),
    strand = "+",
    tx_id = c("t1", "t1", "t2", "t2"),
    gene_id = "g1", gene_name = "g1",
    gene_biotype = "protein_coding", type = "exon"
  )
  res <- Signac:::GetGeneRanges(gr)
  expect_equal(length(res), 1)
  expect_equal(start(res), 1000)
  expect_equal(end(res), 24200)
})

test_that("GetGeneRanges drops genes spanning several chromosomes", {
  # geneA is confined to chr1; PAR1 is annotated on both chrX and chrY, which
  # would otherwise collapse to a nonsense chrX range absorbing the chrY end
  gr <- GRanges(
    c("chr1", "chr1", "chrX", "chrX", "chrY", "chrY"),
    IRanges(c(10000, 20000, 1000, 2000, 50000, 51000),
            c(10500, 20500, 1200, 2200, 50200, 51200)),
    strand = "+"
  )
  gr$tx_id <- c("txA", "txA", "txX", "txX", "txY", "txY")
  gr$gene_id <- c("ENSGA", "ENSGA", "PAR1", "PAR1", "PAR1", "PAR1")
  gr$gene_name <- c("geneA", "geneA", "PAR1", "PAR1", "PAR1", "PAR1")
  gr$gene_biotype <- "protein_coding"
  gr$type <- "exon"

  expect_warning(
    res <- Signac:::GetGeneRanges(gr),
    regexp = "more than one chromosome"
  )
  expect_equal(res$gene_id, "ENSGA")
  expect_equal(as.character(seqnames(res)), "chr1")
  expect_equal(end(res), 20500)
})

test_that("GetTranscriptRanges drops transcripts spanning several chromosomes", {
  gr <- GRanges(
    c("chr1", "chr1", "chrX", "chrY"),
    IRanges(c(10000, 20000, 1000, 50000), c(10500, 20500, 1200, 50200)),
    strand = "+"
  )
  gr$tx_id <- c("txA", "txA", "txBAD", "txBAD")
  gr$gene_id <- "ENSGA"
  gr$gene_name <- "geneA"
  gr$gene_biotype <- "protein_coding"
  gr$type <- "exon"

  expect_warning(
    res <- Signac:::GetTranscriptRanges(gr),
    regexp = "more than one chromosome"
  )
  expect_equal(res$tx_id, "txA")
})

test_that("no multi-chromosome warning on a well formed annotation", {
  ann <- Annotation(atac_small)
  expect_no_warning(Signac:::GetGeneRanges(ranges = ann))
  expect_no_warning(Signac:::GetTranscriptRanges(ranges = ann))
})

test_that("GetTranscriptRanges returns one range per transcript", {
  gr <- GRanges(
    "chr1", IRanges(c(1000, 5000, 20000, 24000), c(1200, 5200, 20200, 24200)),
    strand = "+",
    tx_id = c("t1", "t1", "t2", "t2"),
    gene_id = "g1", gene_name = "g1",
    gene_biotype = "protein_coding", type = "exon"
  )
  res <- Signac:::GetTranscriptRanges(gr)
  expect_equal(length(res), 2)
  expect_setequal(res$tx_id, c("t1", "t2"))
  expect_equal(start(res)[res$tx_id == "t1"], 1000)
  expect_equal(end(res)[res$tx_id == "t1"], 5200)
  expect_equal(start(res)[res$tx_id == "t2"], 20000)
  expect_true(all(res$gene_name == "g1"))
})

test_that("GetTranscriptRanges gives the correct TSS on both strands", {
  # a minus-strand transcript's TSS is its highest coordinate
  gr <- GRanges(
    "chr1", IRanges(c(1000, 5000, 1000, 5000), c(1200, 5200, 1200, 5200)),
    strand = c("+", "+", "-", "-"),
    tx_id = c("plus", "plus", "minus", "minus"),
    gene_id = c("gp", "gp", "gm", "gm"),
    gene_name = c("gp", "gp", "gm", "gm"),
    gene_biotype = "protein_coding", type = "exon"
  )
  tx <- Signac:::GetTranscriptRanges(gr)
  tss <- GenomicRanges::resize(tx, width = 1, fix = "start")
  expect_equal(start(tss)[tx$tx_id == "plus"], 1000)
  expect_equal(start(tss)[tx$tx_id == "minus"], 5200)
})

test_that("GetTranscriptRanges falls back to gene_id without a transcript column", {
  gr <- GRanges(
    "chr1", IRanges(c(1000, 5000), c(1200, 5200)),
    strand = "+",
    gene_id = c("g1", "g2"), gene_name = c("g1", "g2"),
    gene_biotype = "protein_coding"
  )
  res <- Signac:::GetTranscriptRanges(gr)
  expect_equal(length(res), 2)
  expect_setequal(res$tx_id, c("g1", "g2"))
})

test_that("GetTranscriptRanges uses transcript_id when tx_id is absent", {
  gr <- GRanges(
    "chr1", IRanges(c(1000, 5000), c(1200, 5200)),
    strand = "+",
    transcript_id = c("t1", "t1"),
    gene_id = "g1", gene_name = "g1",
    gene_biotype = "protein_coding"
  )
  res <- Signac:::GetTranscriptRanges(gr)
  expect_equal(length(res), 1)
  expect_equal(res$tx_id, "t1")
  expect_equal(start(res), 1000)
  expect_equal(end(res), 5200)
})

test_that("GetTSSPositions returns TSS GRanges", {
  ann <- Annotation(atac_small)
  tss <- GetTSSPositions(ranges = ann)
  expect_s4_class(tss, "GRanges")
  expect_true(all(width(tss) == 1))
})

# MatchRegionStats -------------------------------------------------------------

test_that("MatchRegionStats matches GC distribution", {
  set.seed(1)
  # GC-rich query against a pool that is mostly GC-poor: the matched set should
  # track the query's GC content, not the pool average
  pool <- data.frame(
    GC.percent = c(runif(60, 0.30, 0.45), runif(40, 0.55, 0.70)),
    row.names = paste0("p", 1:100)
  )
  query <- data.frame(
    GC.percent = runif(8, 0.60, 0.70), row.names = paste0("q", 1:8)
  )
  out <- MatchRegionStats(
    meta.feature = pool, query.feature = query, n = 20, verbose = FALSE
  )
  expect_type(out, "character")
  expect_equal(length(out), 20)
  # selected features come from the pool, without repeats
  expect_true(all(out %in% rownames(pool)))
  expect_false(any(duplicated(out)))
  # the matched set is GC-rich like the query, not like the (GC-poor) pool
  selected.gc <- pool[out, "GC.percent"]
  expect_true(all(selected.gc > 0.5))
  expect_lt(
    abs(mean(selected.gc) - mean(query$GC.percent)),
    abs(mean(pool$GC.percent) - mean(query$GC.percent))
  )
})

# ExtractCell / ExtractField / IsMatrixEmpty / isRemote -----------------------

test_that("ExtractCell works", {
  expect_equal(
    object = ExtractCell(x = "chr1\t1\t300\tTGCA\t1"), expected = "TGCA"
  )
})

test_that("ExtractCell returns NULL for empty input", {
  expect_null(ExtractCell(character(0)))
})

test_that("ExtractCell returns positions when requested", {
  out <- ExtractCell(
    x = c("chr1\t1\t100\tBC1\t1", "chr1\t100\t200\tBC2\t1"),
    positions = TRUE
  )
  expect_named(out, c("cell", "start", "end"))
  expect_equal(out$cell, c("BC1", "BC2"))
})

test_that("ExtractField extracts and joins fields", {
  expect_equal(ExtractField("a_b_c", field = 2), "b")
  expect_equal(ExtractField("a_b_c_d", field = "1,3"), "a_c")
})

test_that("ExtractField handles non-default delimiter", {
  expect_equal(ExtractField("a;b;c", field = 2, delim = ";"), "b")
})

test_that("isRemote detects HTTP/FTP paths", {
  expect_true(Signac:::isRemote("https://example.com/file"))
  expect_true(Signac:::isRemote("http://example.com"))
  expect_true(Signac:::isRemote("ftp://example.com"))
  expect_false(Signac:::isRemote("/local/file.tsv"))
})

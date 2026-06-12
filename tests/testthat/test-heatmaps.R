library(GenomicRanges)
library(SeuratObject)

setup_obj_frag <- function() {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  frags <- CreateFragmentObject(
    path = fpath, cells = colnames(atac_small),
    tolerance = 0.5, verbose = FALSE
  )
  obj <- atac_small
  Fragments(obj) <- frags
  Idents(obj) <- "cluster"
  obj
}

# RegionMatrix -----------------------------------------------------------------

test_that("RegionMatrix produces expected matrix structure", {
  obj <- setup_obj_frag()
  regions <- head(granges(obj), 5)
  rm <- RegionMatrix(
    object = obj, regions = regions,
    upstream = 100, downstream = 100, verbose = FALSE
  )
  expect_named(rm, c("matrix", "parameters"))
  expect_equal(rm$parameters$upstream, 100)
  expect_equal(rm$parameters$downstream, 100)
  # one matrix per group identity
  expect_setequal(
    names(rm$matrix),
    as.character(unique(Idents(obj)))
  )
  # matrix dimensions: nrow = regions, ncol = upstream + downstream + 1
  for (m in rm$matrix) {
    expect_equal(dim(m), c(length(regions), 201))
  }
})

test_that("RegionMatrix errors when no fragments", {
  expect_error(
    RegionMatrix(
      object = atac_small,
      regions = head(granges(atac_small), 5),
      upstream = 100, downstream = 100, verbose = FALSE
    ),
    regexp = "No fragment files present in assay"
  )
})

test_that("RegionMatrix errors on Assay5 (non ChromatinAssay5)", {
  obj <- atac_small
  obj[["normal"]] <- SeuratObject::CreateAssay5Object(
    counts = matrix(0, 3, ncol(atac_small),
      dimnames = list(c("a", "b", "c"), colnames(atac_small)))
  )
  expect_error(
    RegionMatrix(
      object = obj, regions = head(granges(obj), 5),
      assay = "normal", verbose = FALSE
    ),
    regexp = "ChromatinAssay5"
  )
})

test_that("RegionMatrix.default validates list of Fragment objects", {
  expect_error(
    RegionMatrix(
      object = list("not_a_fragment"),
      regions = GRanges("chr1", IRanges(1, 1)),
      verbose = FALSE
    ),
    regexp = "Fragment objects"
  )
})

# RegionPlot -------------------------------------------------------------------

test_that("RegionPlot returns ggplot built from RegionMatrix output", {
  obj <- setup_obj_frag()
  rm <- RegionMatrix(
    object = obj, regions = head(granges(obj), 5),
    upstream = 100, downstream = 100, verbose = FALSE
  )
  p <- RegionPlot(object = rm, window = 10)
  expect_s3_class(p, "ggplot")
  expect_gt(length(p$layers), 0)
  # plot data covers the user-requested window range
  expect_setequal(
    unique(p$data$group),
    names(rm$matrix)
  )
})

# RegionHeatmap ----------------------------------------------------------------

test_that("RegionHeatmap returns plot with normalize = TRUE", {
  obj <- setup_obj_frag()
  rm <- RegionMatrix(
    object = obj, regions = head(granges(obj), 5),
    upstream = 100, downstream = 100, verbose = FALSE
  )
  p <- RegionHeatmap(
    object = rm,
    upstream = 50, downstream = 50, window = 10, min.counts = 0
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))

  # normalize divides each group's matrix by its cell count: verify the
  # plotted values are scaled by cells-per-group relative to normalize = FALSE.
  # max.cutoff = NA disables the per-plot quantile clipping so the comparison
  # is exact.
  pn <- RegionHeatmap(
    object = rm, upstream = 50, downstream = 50, window = 10,
    min.counts = 0, normalize = TRUE, max.cutoff = NA
  )
  pu <- RegionHeatmap(
    object = rm, upstream = 50, downstream = 50, window = 10,
    min.counts = 0, normalize = FALSE, max.cutoff = NA
  )
  dn <- pn[[1]]$data
  du <- pu[[1]]$data
  # there should be at least one nonzero, normalized value to compare
  expect_gt(sum(dn$value > 0), 0)
  merged <- merge(
    du, dn,
    by = c("bin", "name", "group"), suffixes = c(".u", ".n")
  )
  merged <- merged[merged$value.n != 0, ]
  expect_gt(nrow(merged), 0)
  # unnormalized / normalized should equal the number of cells in each group
  expected.factor <- as.numeric(rm$parameters$cells[merged$group])
  expect_equal(merged$value.u / merged$value.n, expected.factor)
  # normalized values are strictly smaller (every group has > 1 cell here)
  expect_lt(max(dn$value), max(du$value))
})

test_that("RegionHeatmap with normalize = FALSE", {
  obj <- setup_obj_frag()
  rm <- RegionMatrix(
    object = obj, regions = head(granges(obj), 5),
    upstream = 100, downstream = 100, verbose = FALSE
  )
  p <- RegionHeatmap(
    object = rm,
    upstream = 50, downstream = 50, window = 10,
    min.counts = 0, normalize = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

# H6/H7 regression tests: coincident-insertion counting and strand mirroring
# Build a small bgzipped + tabix-indexed fragment file from a data frame of
# fragment coordinates (chr, start, end), all assigned to one cell.
make_indexed_fragments <- function(starts, ends) {
  frag_df <- data.frame(
    chr = "chr1", start = starts, end = ends, bc = "cell1", n = 1
  )
  frag_df <- frag_df[order(frag_df$start), ]
  tsv <- tempfile(fileext = ".tsv")
  write.table(
    x = frag_df, file = tsv, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  gz <- tempfile(fileext = ".tsv.gz")
  Rsamtools::bgzip(file = tsv, dest = gz, overwrite = TRUE)
  invisible(Rsamtools::indexTabix(file = gz, format = "bed"))
  unlink(tsv)
  CreateFragmentObject(
    path = gz, cells = c(cell1 = "cell1"), validate = FALSE, verbose = FALSE
  )
}

region_matrix <- function(frags, strand) {
  region <- GenomicRanges::GRanges(
    seqnames = "chr1", ranges = IRanges::IRanges(start = 5000, end = 5000),
    strand = strand
  )
  RegionMatrix(
    object = frags, regions = region, group.by = c(cell1 = "a"),
    upstream = 1000, downstream = 1000, verbose = FALSE
  )$matrix$a
}

test_that("RegionMatrix sums coincident Tn5 insertions", {
  # two fragments share start 4500; all positions lie within the region
  frags <- make_indexed_fragments(
    starts = c(4500, 4500, 4900, 5200, 5500),
    ends = c(4600, 4700, 5100, 5300, 5800)
  )
  plus <- region_matrix(frags, "+")
  minus <- region_matrix(frags, "-")

  # 5 fragments x 2 ends = 10 insertions; the shared start must be counted
  # twice (a vectorised increment collapses duplicates to give 9 / max 1)
  expect_equal(sum(plus), 10)
  expect_equal(max(plus), 2)
  expect_equal(sum(minus), 10)
  expect_equal(max(minus), 2)
})

test_that("RegionMatrix keeps + and - strands centre-aligned", {
  # insertion positions symmetric about the region centre (5000): a correct,
  # centre-aligned strand reversal maps this set onto itself, so the plus and
  # minus matrices must be identical
  frags <- make_indexed_fragments(
    starts = c(4500, 4800),
    ends = c(5500, 5200)
  )
  plus <- region_matrix(frags, "+")
  minus <- region_matrix(frags, "-")
  expect_equal(minus, plus)
})

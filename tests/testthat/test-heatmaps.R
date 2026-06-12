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

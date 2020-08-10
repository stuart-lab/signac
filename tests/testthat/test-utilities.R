library(GenomicRanges)

test_that("AverageCounts works", {
  expect_equal(
    object = as.vector(x = AverageCounts(object = atac_small)),
    expected = c(162.80, 118.74),
    tolerance = 1 / 1000
  )
})

test_that("CellsPerGroup works", {
  expect_equal(
    object = as.vector(x = CellsPerGroup(object = atac_small)),
    expected = c(50, 50, 0)
  )
})

test_that("GRanges conversion works", {
  correct_string <- c("chr1-1-10", "chr2-12-3121")
  correct_ranges <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(start = c(1, 12), end = c(10, 3121))
  )
  granges <- StringToGRanges(regions = correct_string)
  string_ranges <- GRangesToString(grange = correct_ranges)
  expect_equal(object = granges, expected = correct_ranges)
  expect_equal(object = string_ranges, expected = correct_string)
})

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
    ranges = IRanges(
      start = c(1000, 1200),
      end = c(10000, 312100)
    ),
    strand = c("+", "-")
  )
  correct_extended <- GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges(
      start = c(500, 200),
      end = c(11000, 312600)
    ),
    strand = c("+", "-")
  )
  extended_granges <- Extend(x = granges, upstream = 500, downstream = 1000)
  expect_equal(object = extended_granges, expected = correct_extended)
})

test_that("ExtractCell works", {
  expect_equal(
    object = ExtractCell(x = "chr1\t1\t300\tTGCA\t1"), expected = "TGCA"
  )
})

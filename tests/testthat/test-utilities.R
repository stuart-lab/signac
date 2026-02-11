library(GenomicRanges)

test_that("AverageCounts works", {
  expect_equal(
    object = as.vector(x = AverageCounts(object = atac_small)),
    expected = 30.11,
    tolerance = 1 / 1000
  )
})

test_that("CellsPerGroup works", {
  Idents(atac_small) <- "cluster"
  expect_equal(
    object = as.vector(x = CellsPerGroup(object = atac_small)),
    expected = c(55, 45, 0)
  )
})

test_that("SortIdents works",{
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

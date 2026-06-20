library(GenomicRanges)
library(SeuratObject)

test_that("granges works on GRangesAssay", {
  gr <- granges(x = atac_small[["peaks"]])
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), nrow(atac_small[["peaks"]]))
})

test_that("granges works on Seurat", {
  gr <- granges(x = atac_small)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), nrow(atac_small[["peaks"]]))
})

test_that("granges rejects use.mcols on GRangesAssay", {
  expect_error(
    granges(x = atac_small[["peaks"]], use.mcols = TRUE),
    regexp = "use.mcols"
  )
})

test_that("granges rejects use.mcols on Seurat", {
  expect_error(
    granges(x = atac_small, use.mcols = TRUE),
    regexp = "use.mcols"
  )
})

test_that("seqinfo and seqnames work on GRangesAssay", {
  si <- seqinfo(x = atac_small[["peaks"]])
  expect_s4_class(si, "Seqinfo")
  sn <- seqnames(x = atac_small[["peaks"]])
  expect_s4_class(sn, "Rle")
  expect_true(all(as.character(sn) == "chr1"))
})

test_that("seqlevels, seqlengths, genome, isCircular work on GRangesAssay", {
  # atac_small uses an unspecified genome with one chromosome and no seqlengths
  expect_equal(seqlevels(x = atac_small[["peaks"]]), "chr1")
  expect_true(is.na(seqlengths(x = atac_small[["peaks"]])[["chr1"]]))
  expect_true(is.na(genome(x = atac_small[["peaks"]])[["chr1"]]))
  expect_true(is.na(isCircular(x = atac_small[["peaks"]])[["chr1"]]))
})

test_that("seqinfo and friends on Seurat match assay-level values", {
  expect_s4_class(seqinfo(x = atac_small), "Seqinfo")
  expect_equal(seqlevels(x = atac_small), "chr1")
  expect_s4_class(seqnames(x = atac_small), "Rle")
  # delegation: Seurat method must return identical to the assay method
  expect_identical(
    seqinfo(x = atac_small),
    seqinfo(x = atac_small[["peaks"]])
  )
  expect_identical(
    seqlengths(x = atac_small),
    seqlengths(x = atac_small[["peaks"]])
  )
  expect_identical(
    genome(x = atac_small),
    genome(x = atac_small[["peaks"]])
  )
  expect_identical(
    isCircular(x = atac_small),
    isCircular(x = atac_small[["peaks"]])
  )
})

test_that("findOverlaps query=assay against its own granges hits each peak once", {
  peaks <- atac_small[["peaks"]]
  gr <- granges(x = peaks)
  hits <- findOverlaps(query = peaks, subject = gr)
  expect_s4_class(hits, "Hits")
  # peaks are non-overlapping, so each query hits exactly one subject
  expect_equal(length(hits), length(gr))
  co <- countOverlaps(query = peaks, subject = gr)
  expect_equal(co, rep(1L, length(gr)))
})

test_that("findOverlaps on Seurat delegates to assay", {
  hits1 <- findOverlaps(query = atac_small, subject = granges(atac_small))
  hits2 <- findOverlaps(
    query = atac_small[["peaks"]],
    subject = granges(atac_small[["peaks"]])
  )
  expect_identical(hits1, hits2)
})

test_that("nearest/precede/follow/distance methods delegate to GRanges", {
  peaks <- atac_small[["peaks"]]
  gr <- granges(peaks)
  # GRangesAssay methods should produce results identical to running the
  # underlying GenomicRanges function on the extracted granges.
  expect_identical(nearest(x = peaks, subject = gr), nearest(x = gr, subject = gr))
  expect_identical(precede(x = peaks, subject = gr), precede(x = gr, subject = gr))
  expect_identical(follow(x = peaks, subject = gr),  follow(x = gr, subject = gr))
  expect_identical(distance(x = peaks, y = gr), distance(x = gr, y = gr))
  expect_identical(
    distanceToNearest(x = peaks, subject = gr),
    distanceToNearest(x = gr, subject = gr)
  )
})

test_that("nearest with single argument delegates to granges", {
  peaks <- atac_small[["peaks"]]
  expect_identical(nearest(x = peaks), nearest(x = granges(peaks)))
})

test_that("range/reduce/disjoin/gaps/coverage delegate to underlying granges", {
  peaks <- atac_small[["peaks"]]
  gr <- granges(peaks)
  expect_identical(range(x = peaks), range(x = gr))
  expect_identical(reduce(x = peaks), reduce(x = gr))
  expect_identical(disjoin(x = peaks), disjoin(x = gr))
  expect_identical(gaps(x = peaks), gaps(x = gr))
  expect_identical(isDisjoint(x = peaks), isDisjoint(x = gr))
  expect_identical(disjointBins(x = peaks), disjointBins(x = gr))
  expect_identical(coverage(x = peaks), coverage(x = gr))
})

test_that("Seurat-level dispatch matches GRangesAssay dispatch", {
  gr <- granges(atac_small)
  # delegation: Seurat method == calling on the default assay
  expect_identical(
    nearest(x = atac_small, subject = gr),
    nearest(x = atac_small[["peaks"]], subject = gr)
  )
  expect_identical(range(x = atac_small), range(x = atac_small[["peaks"]]))
  expect_identical(coverage(x = atac_small), coverage(x = atac_small[["peaks"]]))
})

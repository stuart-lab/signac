test_that("FeatureMatrix works on grange on diff seqnames", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  features <- suppressWarnings(c(
    granges(atac_small),
    GenomicRanges::GRanges(
      seqnames = "fake_chr",
      ranges = IRanges::IRanges(start = 1, end = 1000)
    )
  ))
  expect_warning(mat <- FeatureMatrix(
    fragments = fragments,
    features = features,
    verbose = FALSE,
    fragtk = FALSE
  ))
  expect_equal(dim(mat), c(100, 76))
  mat <- FeatureMatrix(
    fragments = fragments,
    features = features,
    keep_all_features = TRUE,
    fragtk = FALSE,
    verbose = FALSE
  )
  expect_equal(dim(mat), c(101, 76))
})

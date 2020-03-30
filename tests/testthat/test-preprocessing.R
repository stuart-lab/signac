suppressWarnings(RNGversion(vstr = "3.5.3"))

test_that("BinarizeCounts works", {
  set.seed(1)
  # matrix
  mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
  bin_mat <- BinarizeCounts(object = mat)

  # sparse matrix
  mat_sparse <- as(object = mat, Class = 'dgCMatrix')
  bin_mat_sparse <- BinarizeCounts(object = mat_sparse)

  expect_equal(object = as.vector(bin_mat[1,]), expected = c(0,1,0,1,1))
  expect_equal(object = as.vector(bin_mat_sparse[1,]), expected = c(0,1,0,1,1))
})

test_that("DownsampleFeatures works", {
  set.seed(1)
  atac_ds <- DownsampleFeatures(object = atac_small, n = 5, verbose = FALSE)
  expect_equal(
    object = VariableFeatures(object = atac_ds),
    expected = c("chr1:1804025-1804468","chr1:2221250-2223390","chr1:6074646-6075340",
                 "chr1:8931013-8932013","chr1:1562519-1567986")
  )
})

test_that("FindTopFeatures works", {
  VariableFeatures(atac_small) <- NULL
  atac_small <- FindTopFeatures(object = atac_small)
  expect_equal(
    object = head(VariableFeatures(object = atac_small)),
    expected = c("chr1:1549446-1552535",
                 "chr1:1051006-1053102",
                 "chr1:1240091-1245762",
                 "chr1:1333514-1336003",
                 "chr1:1309645-1311492",
                 "chr1:928630-937949")
  )
})

test_that("FRiP works", {
  atac_small <- FRiP(
    object = atac_small,
    bin.assay = 'bins',
    peak.assay = 'peaks',
    chromosome = NULL
  )
  expect_equal(
    object = as.vector(x = head(atac_small$FRiP)),
    expected = c(1.4090909,1.4444444,1.1363636,0.8798077,1.0000000,0.2166667),
    tolerance = 1/100000
  )
})

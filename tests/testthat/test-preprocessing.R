suppressWarnings(RNGversion(vstr = "3.5.3"))

test_that("BinarizeCounts works", {
  set.seed(1)
  # matrix
  mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
  bin_mat <- BinarizeCounts(object = mat)

  # sparse matrix
  mat_sparse <- as.sparse(x = mat)
  bin_mat_sparse <- BinarizeCounts(object = mat_sparse)

  expect_equal(object = as.vector(bin_mat[1,]), expected = c(0,1,0,1,1))
  expect_equal(object = as.vector(bin_mat_sparse[1,]), expected = c(0,1,0,1,1))
})

test_that("DownsampleFeatures works", {
  set.seed(1)
  atac_ds <- DownsampleFeatures(object = atac_small, n = 5, verbose = FALSE)
  expect_equal(
    object = VariableFeatures(object = atac_ds),
    expected = c("chr21:35013350-35016192","chr16:637973-640868",
                 "chr5:81045135-81048185","chr7:128863917-128866113","chr1:37938558-37953752")
  )
})

test_that("FindTopFeatures works", {
  VariableFeatures(atac_small) <- NULL
  atac_small <- FindTopFeatures(object = atac_small)
  expect_equal(
    object = head(VariableFeatures(object = atac_small)),
    expected = c("chr19:39887431-39927397","chr17:7736488-7763641","chr19:13935996-13962675",
                 "chr19:1247604-1277130","chr19:45969649-45989808","chr19:12892327-12913486")
  )
})

test_that("FRiP works", {
  atac_small$FRiP <- NULL
  atac_small <- FRiP(object = atac_small, bin.assay = 'bins', peak.assay = 'peaks', chromosome = NULL)
  expect_equal(
    object = as.vector(x = head(atac_small$FRiP)),
    expected = c(0.8992537,1.0020060,0.8639241,0.8753388,0.9545455,0.9678530),
    tolerance = 1/100000
  )
})

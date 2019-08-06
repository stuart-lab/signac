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
    expected = c("chr1:1812400-1813494","chr1:2236104-2237259",
                 "chr1:6074646-6075340","chr1:8933459-8934461",
                 "chr1:1553343-1553743")
  )
})

test_that("FindTopFeatures works", {
  VariableFeatures(atac_small) <- NULL
  atac_small <- FindTopFeatures(object = atac_small)
  expect_equal(
    object = head(VariableFeatures(object = atac_small)),
    expected = c("chr1:2157847-2188813","chr1:6843960-6846894","chr1:2471903-2481288",
                 "chr1:3815928-3820356","chr1:2515241-2519350","chr1:6051145-6055407")
  )
})

test_that("FRiP works", {
  atac_small <- FRiP(object = atac_small, bin.assay = 'bins', peak.assay = 'peaks', chromosome = NULL)
  expect_equal(
    object = as.vector(x = head(atac_small$FRiP)),
    expected = c(1.4090909,1.4814815,1.1666667,0.9759615,1.0384615,0.2666667),
    tolerance = 1/100000
  )
})

suppressWarnings(RNGversion(vstr = "3.5.3"))

test_that("BinarizeCounts works", {
  set.seed(1)
  # matrix
  mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
  bin_mat <- BinarizeCounts(object = mat)

  # sparse matrix
  mat_sparse <- as(object = mat, Class = "dgCMatrix")
  bin_mat_sparse <- BinarizeCounts(object = mat_sparse)

  expect_equal(
    object = as.vector(bin_mat[1, ]), expected = c(0, 1, 0, 1, 1)
  )
  expect_equal(
    object = as.vector(bin_mat_sparse[1, ]), expected = c(0, 1, 0, 1, 1)
  )
})

test_that("DownsampleFeatures works", {
  set.seed(1)
  atac_ds <- DownsampleFeatures(object = atac_small, n = 5, verbose = FALSE)
  expect_equal(
    object = Seurat::VariableFeatures(object = atac_ds),
    expected = c(
      "chr1-1804025-1804468",
      "chr1-2221250-2223390",
      "chr1-6074646-6075340",
      "chr1-8931013-8932013",
      "chr1-1562519-1567986"
      )
  )
})

test_that("FindTopFeatures works", {
  VariableFeatures(atac_small) <- NULL
  atac_small <- FindTopFeatures(object = atac_small)
  expect_equal(
    object = head(Seurat::VariableFeatures(object = atac_small)),
    expected = c("chr1-2157847-2188813",
                 "chr1-2471903-2481288",
                 "chr1-6843960-6846894",
                 "chr1-3815928-3820356",
                 "chr1-8935313-8940649",
                 "chr1-2515241-2519350")
  )
})

test_that("FRiP works", {
  atac_small <- FRiP(
    object = atac_small,
    assay = "peaks",
    total.fragments = "fragments"
  )
  expect_equal(
    object = as.vector(x = head(atac_small$FRiP)),
    expected = c(
      0.0069896591, 0.0086400864, 0.0102987567,
      0.0093475841, 0.0081182600, 0.0001859341
    ),
    tolerance = 1 / 100000
  )
})

test_that("FeatureMatrix works", {
  computed_fmat <- readRDS("../testdata/featurematrix.rds")
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fragments <- CreateFragmentObject(
    path = fpath,
    cells = colnames(x = atac_small),
    tolerance = 0.5,
    verbose = FALSE
  )
  fm <- FeatureMatrix(
    fragments = fragments,
    features = granges(atac_small),
    verbose = FALSE
  )
  expect_identical(object = fm, expected = computed_fmat)
})

test_that("NucleosomeSignal works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fragments <- CreateFragmentObject(
    path = fpath,
    cells = colnames(x = atac_small),
    tolerance = 0.5,
    verbose = FALSE
  )
  Fragments(object = atac_small) <- fragments
  ns <- NucleosomeSignal(object = atac_small, verbose = FALSE)
  expect_equal(
    object = as.numeric(x = head(x = ns$nucleosome_signal)),
    expected = c(NaN, NaN, 0.0, 2.5, NaN, NaN)
  )
})

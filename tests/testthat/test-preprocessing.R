suppressWarnings(RNGversion(vstr = "3.5.3"))

test_that("BinarizeCounts works", {
  set.seed(1)
  # matrix
  mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
  bin_mat <- BinarizeCounts(object = mat)

  # sparse matrix
  mat_sparse <- as(object = mat, Class = "CsparseMatrix")
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
    object = SeuratObject::VariableFeatures(object = atac_ds),
    expected = c(
      "chr1:898350-899223",
      "chr1:955190-956101",
      "chr1:1059203-1060060",
      "chr1:1250624-1251529",
      "chr1:860139-860923"
    )
  )
})

test_that("FindTopFeatures works", {
  VariableFeatures(atac_small) <- NULL
  atac_small <- FindTopFeatures(object = atac_small)
  expect_equal(
    object = head(SeuratObject::VariableFeatures(object = atac_small)),
    expected = c(
      "chr1:191183-192084",
      "chr1:270850-271755",
      "chr1:273946-274792",
      "chr1:854732-855551",
      "chr1:877256-878073",
      "chr1:897006-897867"
    )
  )
})

test_that("FRiP works", {
  atac_small <- FRiP(
    object = atac_small,
    assay = "peaks",
    total.fragments = "nCount_peaks"
  )
  expect_equal(
    object = as.vector(x = head(atac_small$FRiP)),
    expected = c(
      1, 1, 1, 1, 1, 1
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
    fragtk = FALSE,
    verbose = FALSE
  )
  expect_identical(object = fm, expected = computed_fmat)
})

test_that("CreateMotifMatrix works", {
  pwm <- readRDS("../testdata/pwm_2motifs.rds")
  genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
  genome <- Rsamtools::FaFile(genome.fasta)
  skip_if_not_installed("motifmatchr")
  motif.matrix <- suppressWarnings(CreateMotifMatrix(
    features = granges(atac_small),
    pwm = pwm,
    genome = genome
  ))
  expect_equal(dim(motif.matrix), c(100, 2))
})

test_that("PearsonResidualVar calculates variance correctly", {
  set.seed(42)
  # Create small test matrix
  mat <- Matrix::rsparsematrix(nrow = 50, ncol = 30, density = 0.2)
  
  result <- PearsonResidualVar(mat, theta = 100, verbose = FALSE)
  
  # Basic checks
  expect_true(is.data.frame(result))
  expect_true("ResidualVariance" %in% colnames(result))
  expect_true("mean" %in% colnames(result))
  expect_true("count" %in% colnames(result))
  expect_true("rank" %in% colnames(result))
})

test_that("PearsonResidualVar batch processing accumulates correctly", {
  set.seed(123)
  mat <- Matrix::rsparsematrix(nrow = 20, ncol = 50, density = 0.15)
  
  # Run with different batch sizes - results should match
  result_batch5 <- PearsonResidualVar(mat, ncell.batch = 5, verbose = FALSE)
  result_batch25 <- PearsonResidualVar(mat, ncell.batch = 25, verbose = FALSE)
  result_batch50 <- PearsonResidualVar(mat, ncell.batch = 50, verbose = FALSE)
  
  # Results should be identical regardless of batch size
  expect_equal(result_batch5$ResidualVariance, 
               result_batch25$ResidualVariance, 
               tolerance = 1e-10)
  expect_equal(result_batch25$ResidualVariance, 
               result_batch50$ResidualVariance, 
               tolerance = 1e-10)
})

test_that("PearsonResidualVar handles edge cases", {
  # Matrix with some all-zero rows
  mat <- Matrix::Matrix(0, nrow = 5, ncol = 10, sparse = TRUE)
  mat[1, c(1, 3, 5)] <- c(2, 4, 6)
  mat[3, c(2, 4)] <- c(1, 3)
  # Rows 2, 4, 5 are all zeros
  
  result <- PearsonResidualVar(mat, verbose = FALSE)
  
  # Should only return rows with non-zero mean (2 rows)
  expect_equal(nrow(result), 5)
  expect_equal(nrow(result[result$mean > 0, ]), 2)
  
  # Check the means are correct for the non-zero rows
  # Row 1 has values [2,0,4,0,6,0,0,0,0,0] so mean = 12/10 = 1.2
  # Row 3 has values [0,1,0,3,0,0,0,0,0,0] so mean = 4/10 = 0.4
  expect_equal(result$mean[1], 1.2, tolerance = 1e-10)
  expect_equal(result$mean[3], 0.4, tolerance = 1e-10)
})

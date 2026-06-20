library(SeuratObject)
library(Matrix)

suppressWarnings(RNGversion(vstr = "3.5.3"))

# LSI / RunTFIDF reference values ---------------------------------------------

test_that("LSI works", {
  set.seed(seed = 1)
  mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
  method1 <- RunTFIDF(object = mat, method = 1)
  method2 <- RunTFIDF(object = mat, method = 2)
  method3 <- RunTFIDF(object = mat, method = 3)
  method4 <- RunTFIDF(object = mat, method = 4)

  expect_equal(
    object = method1[1, ],
    expected = c(0.000000, 7.957927, 0.000000, 7.131699, 8.805025),
    tolerance = 1 / 1000
  )
  expect_equal(
    object = method2[1, ],
    expected = c(0.0000000, 0.1980421, 0.0000000, 0.0866434, 0.4620981),
    tolerance = 1 / 1000
  )
  expect_equal(
    object = method3[1, ],
    expected = c(0.000000, 5.516015, 0.000000, 4.943317, 6.103178),
    tolerance = 1 / 1000
  )
  expect_equal(object = method4[1, ], expected = c(0, 2, 0, 1, 2))

  lsi <- suppressWarnings(RunSVD(object = mat))
  embeddings <- SeuratObject::Embeddings(object = lsi)
  loadings <- SeuratObject::Loadings(object = lsi)

  expect_equal(
    object = as.vector(embeddings[1, ]),
    expected = c(0.51255352, 0.08674426, -1.33604004, 1.18108240),
    tolerance = 1 / 1000
  )
  expect_equal(
    object = as.vector(loadings[1, ]),
    expected = c(-0.4024075, 0.4292469, 0.6463644, 0.1740785),
    tolerance = 1 / 1000
  )
})

# Jaccard ----------------------------------------------------------------------

test_that("Jaccard works", {
  set.seed(1)
  mat <- matrix(
    data = sample(x = c(0, 1), size = 25, replace = TRUE),
    nrow = 5
  )
  jm <- Jaccard(x = mat, y = mat)
  expect_equal(object = jm[1, ], expected = c(1, 1 / 3, 2 / 5, 1 / 3, 0))
})

test_that("Jaccard computes a matrix between two matrices", {
  set.seed(1)
  X <- matrix(sample(c(0, 1), 100, replace = TRUE), 10, 10)
  Y <- matrix(sample(c(0, 1), 100, replace = TRUE), 10, 10)
  res <- Jaccard(X, Y)
  expect_equal(dim(res), c(10, 10))
})

test_that("Jaccard warns on values > 1", {
  X <- matrix(sample(0:3, 100, replace = TRUE), 10, 10)
  Y <- matrix(sample(0:3, 100, replace = TRUE), 10, 10)
  expect_warning(Jaccard(X, Y), regexp = "binarize")
})

# RunSVD -----------------------------------------------------------------------

test_that("RunSVD on matrix returns DimReduc", {
  set.seed(1)
  m <- matrix(rnorm(500), nrow = 50, ncol = 10)
  res <- RunSVD(m, n = 5, verbose = FALSE)
  expect_s4_class(res, "DimReduc")

  embeddings <- Embeddings(object = res)
  # one row per column (cell) of the input, one column per requested component
  expect_equal(object = dim(embeddings), expected = c(10, 5))
  expect_equal(object = Key(object = res), expected = "LSI_")
  expect_equal(
    object = colnames(embeddings),
    expected = paste0("LSI_", 1:5)
  )
  # default scale.embeddings = TRUE => each component has mean 0, SD 1
  expect_equal(
    object = unname(apply(X = embeddings, MARGIN = 2, FUN = mean)),
    expected = rep(0, 5),
    tolerance = 1 / 1000
  )
  expect_equal(
    object = unname(apply(X = embeddings, MARGIN = 2, FUN = sd)),
    expected = rep(1, 5),
    tolerance = 1 / 1000
  )
  expect_equal(
    object = as.vector(embeddings[1, ]),
    expected = c(-1.2330585, 0.2963341, -0.0124394, -0.9656036, 0.5545593),
    tolerance = 1 / 1000
  )
})

test_that("RunSVD with pca = TRUE works", {
  set.seed(1)
  m <- matrix(rnorm(500), nrow = 50, ncol = 10)
  res <- RunSVD(m, n = 5, pca = TRUE, verbose = FALSE)
  expect_s4_class(res, "DimReduc")

  # pca = TRUE uses the PCA_ key (not LSI_)
  expect_equal(object = Key(object = res), expected = "PCA_")

  emb.pca <- Embeddings(object = res)
  expect_equal(object = dim(emb.pca), expected = c(10, 5))
  expect_equal(
    object = colnames(emb.pca),
    expected = paste0("PCA_", 1:5)
  )

  # The pca branch weights embeddings by eigenvalues and does NOT scale them
  # to unit variance, so the result must differ from pca = FALSE. If pca were
  # ignored these embeddings would equal the scaled (unit-variance) LSI result.
  set.seed(1)
  m2 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  res.lsi <- RunSVD(m2, n = 5, pca = FALSE, verbose = FALSE)
  emb.lsi <- Embeddings(object = res.lsi)

  expect_false(
    isTRUE(all.equal(
      target = unname(emb.pca),
      current = unname(emb.lsi),
      tolerance = 1 / 1000
    ))
  )
  # unscaled, eigenvalue-weighted: component SDs are not all 1
  expect_false(
    isTRUE(all.equal(
      target = unname(apply(X = emb.pca, MARGIN = 2, FUN = sd)),
      current = rep(1, 5),
      tolerance = 1 / 1000
    ))
  )
  expect_equal(
    object = as.vector(emb.pca[1, ]),
    expected = c(-1.2143982, 0.4486555, -0.2778468, -0.3807493, -0.6286513),
    tolerance = 1 / 1000
  )
})

test_that("RunSVD with scale.max clipping works", {
  set.seed(1)
  m <- matrix(rnorm(500), nrow = 50, ncol = 10)

  # without clipping the embeddings exceed 0.5 in magnitude, so a scale.max
  # of 0.5 must actually clip them
  res.unclipped <- RunSVD(m, n = 5, verbose = FALSE)
  expect_gt(max(abs(Embeddings(object = res.unclipped))), 0.5)

  res <- RunSVD(m, n = 5, scale.max = 0.5, verbose = FALSE)
  expect_s4_class(res, "DimReduc")

  embeddings <- Embeddings(object = res)
  expect_equal(object = dim(embeddings), expected = c(10, 5))
  # clipping bound is respected
  expect_lte(max(abs(embeddings)), 0.5)
  # and clipping actually bit (some values were pushed to the bound)
  expect_equal(object = max(abs(embeddings)), expected = 0.5)
})

test_that("RunSVD on Seurat returns updated Seurat", {
  obj <- atac_small
  VariableFeatures(obj) <- rownames(obj[["peaks"]])
  res <- RunSVD(obj, n = 5, verbose = FALSE)
  expect_s4_class(res, "Seurat")
  expect_true("lsi" %in% names(res@reductions))

  embeddings <- Embeddings(object = res, reduction = "lsi")
  # one row per cell, one column per requested component
  expect_equal(object = dim(embeddings), expected = c(ncol(obj), 5))
  expect_equal(object = Key(object = res[["lsi"]]), expected = "LSI_")
  expect_equal(
    object = as.vector(embeddings[1, ]),
    expected = c(-0.9426501, 1.0001435, 0.7316286, 0.4466094, 0.1253871),
    tolerance = 1 / 1000
  )
})

test_that("RunSVD on assay returns DimReduc", {
  a <- atac_small[["peaks"]]
  VariableFeatures(a) <- rownames(a)
  res <- RunSVD(a, n = 5, verbose = FALSE)
  expect_s4_class(res, "DimReduc")

  embeddings <- Embeddings(object = res)
  expect_equal(object = dim(embeddings), expected = c(ncol(a), 5))
  expect_equal(object = Key(object = res), expected = "LSI_")
  expect_equal(
    object = as.vector(embeddings[1, ]),
    expected = c(-0.9426501, 1.0001435, 0.7316286, 0.4466094, 0.1253871),
    tolerance = 1 / 1000
  )
})

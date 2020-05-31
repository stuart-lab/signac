suppressWarnings(RNGversion(vstr = "3.5.3"))

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
  expect_equal(
    object = method4[1, ],
    expected = c(0, 2, 0, 1, 2)
  )

  lsi <- suppressWarnings(RunSVD(object = mat))
  embeddings <- Seurat::Embeddings(object = lsi)
  loadings <- Seurat::Loadings(object = lsi)

  expect_equal(
    object = as.vector(embeddings[1, ]),
    expected = c(0.51255352, -0.08674426, 1.33604004, 1.18108240),
    tolerance = 1 / 1000
  )
  expect_equal(
    object = as.vector(loadings[1, ]),
    expected = c(-0.4024075, -0.4292469, -0.6463644, 0.1740785),
    tolerance = 1 / 1000
  )
})

test_that("Jaccard works", {
  set.seed(1)
  mat <- matrix(data = sample(x = c(0, 1), size = 25, replace = TRUE), nrow = 5)
  jm <- Jaccard(x = mat, y = mat)
  expect_equal(object = jm[1, ], expected = c(1, 1 / 3, 2 / 5, 1 / 3, 0))
})

test_that("FoldChange dispatches on ChromatinAssay5", {
  skip_if_not_installed("Seurat")
  cells.1 <- colnames(atac_small)[atac_small$cluster == 1]
  cells.2 <- colnames(atac_small)[atac_small$cluster == 2]
  res <- Seurat::FoldChange(
    object = atac_small[["peaks"]],
    cells.1 = cells.1, cells.2 = cells.2
  )
  expect_s3_class(res, "data.frame")
  expect_true("avg_log2FC" %in% colnames(res))
  expect_equal(nrow(res), nrow(atac_small[["peaks"]]))
  expect_equal(rownames(res), rownames(atac_small[["peaks"]]))
  expect_true(is.numeric(res$avg_log2FC))
  expect_false(any(is.nan(res$avg_log2FC)))
})

test_that("FoldChange ChromatinAssay5 respects fc.name", {
  skip_if_not_installed("Seurat")
  cells.1 <- colnames(atac_small)[atac_small$cluster == 1]
  cells.2 <- colnames(atac_small)[atac_small$cluster == 2]
  res <- Seurat::FoldChange(
    object = atac_small[["peaks"]],
    cells.1 = cells.1,
    cells.2 = cells.2,
    fc.name = "my_fc"
  )
  expect_true("my_fc" %in% colnames(res))
})

test_that("FoldChange ChromatinAssay5 supports custom mean function and base", {
  skip_if_not_installed("Seurat")
  cells.1 <- colnames(atac_small)[atac_small$cluster == 1]
  cells.2 <- colnames(atac_small)[atac_small$cluster == 2]
  res <- Seurat::FoldChange(
    object = atac_small[["peaks"]],
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = function(x) Matrix::rowMeans(x),
    base = exp(1)
  )
  expect_s3_class(res, "data.frame")
  expect_true("avg_logFC" %in% colnames(res))
})

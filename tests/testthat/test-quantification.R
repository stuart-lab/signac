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
    object = fragments,
    features = features,
    verbose = FALSE,
    fragtk = FALSE
  ))
  expect_equal(dim(mat), c(100, 76))
  mat <- FeatureMatrix(
    object = fragments,
    features = features,
    keep_all_features = TRUE,
    fragtk = FALSE,
    verbose = FALSE
  )
  expect_equal(dim(mat), c(101, 76))
})

test_that("FeatureMatrix returns a BPCells IterableMatrix when bpcells = TRUE", {
  skip_if_not_installed("BPCells")
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  features <- granges(atac_small)
  baseline <- FeatureMatrix(
    object = fragments,
    features = features,
    fragtk = FALSE,
    verbose = FALSE
  )
  persist.dir <- tempfile(pattern = "signac_bpcells_test_")
  bp <- FeatureMatrix(
    object = fragments,
    features = features,
    fragtk = FALSE,
    bpcells = TRUE,
    bpcells.dir = persist.dir,
    verbose = FALSE
  )
  expect_s4_class(bp, "IterableMatrix")
  expect_true(dir.exists(persist.dir))
  expect_equal(dim(bp), dim(baseline))
  expect_equal(rownames(bp), rownames(baseline))
  expect_equal(colnames(bp), colnames(baseline))
  # materialize and compare values
  bp.sparse <- as(bp, "dgCMatrix")
  expect_equal(
    as.matrix(bp.sparse[, colnames(baseline)]),
    as.matrix(baseline)
  )
})

test_that("FeatureMatrix dispatches through ChromatinAssay with bpcells = TRUE", {
  skip_if_not_installed("BPCells")
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  frags <- CreateFragmentObject(
    path = fpath, cells = colnames(x = atac_small),
    tolerance = 0.5, verbose = FALSE
  )
  obj <- atac_small
  Fragments(obj) <- NULL
  Fragments(obj) <- frags
  persist.dir <- tempfile(pattern = "signac_bpcells_assay_")
  bp <- FeatureMatrix(
    object = obj,
    features = granges(atac_small),
    fragtk = FALSE,
    bpcells = TRUE,
    bpcells.dir = persist.dir,
    verbose = FALSE
  )
  expect_s4_class(bp, "IterableMatrix")
  expect_true(dir.exists(persist.dir))
  expect_equal(ncol(bp), ncol(atac_small))
})

test_that("GenomeBinMatrix returns IterableMatrix when bpcells = TRUE", {
  skip_if_not_installed("BPCells")
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(
    path = fpath, cells = colnames(x = atac_small), verbose = FALSE
  )
  genome <- c(chr1 = 780007)
  persist.dir <- tempfile(pattern = "signac_bpcells_bin_")
  bp <- GenomeBinMatrix(
    fragments = fragments,
    genome = genome,
    binsize = 100000,
    bpcells = TRUE,
    bpcells.dir = persist.dir,
    verbose = FALSE
  )
  expect_s4_class(bp, "IterableMatrix")
  expect_true(dir.exists(persist.dir))
})

test_that("GeneActivity returns IterableMatrix when bpcells = TRUE", {
  skip_if_not_installed("BPCells")
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  frags <- CreateFragmentObject(
    path = fpath, cells = colnames(x = atac_small),
    tolerance = 0.5, verbose = FALSE
  )
  obj <- atac_small
  Fragments(obj) <- NULL
  Fragments(obj) <- frags
  persist.dir <- tempfile(pattern = "signac_bpcells_ga_")
  bp <- suppressWarnings(GeneActivity(
    object = obj,
    fragtk = FALSE,
    bpcells = TRUE,
    bpcells.dir = persist.dir,
    verbose = FALSE
  ))
  expect_s4_class(bp, "IterableMatrix")
  expect_true(dir.exists(persist.dir))
  # the intermediate tempdir should have been cleaned up
  expect_equal(
    length(list.files(tempdir(), pattern = "^signac_geneactivity_")),
    0
  )
})

test_that("FeatureMatrix validates bpcells arguments", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  features <- granges(atac_small)
  # bpcells = TRUE without bpcells.dir should error
  expect_error(
    FeatureMatrix(
      object = fragments, features = features, fragtk = FALSE,
      bpcells = TRUE, verbose = FALSE
    ),
    regexp = "bpcells.dir.*must be supplied"
  )
  # non-empty existing dir for bpcells.dir should be rejected
  existing <- tempfile(pattern = "signac_bpcells_existing_")
  dir.create(existing)
  file.create(file.path(existing, "placeholder"))
  expect_error(
    FeatureMatrix(
      object = fragments, features = features, fragtk = FALSE,
      bpcells = TRUE, bpcells.dir = existing, verbose = FALSE
    ),
    regexp = "already exists and is not empty"
  )
  # warning when bpcells.dir supplied without bpcells = TRUE
  expect_warning(
    FeatureMatrix(
      object = fragments, features = features, fragtk = FALSE,
      bpcells = FALSE, bpcells.dir = tempfile(), verbose = FALSE
    ),
    regexp = "bpcells.dir` is ignored"
  )
})

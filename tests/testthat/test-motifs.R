pwm <- readRDS("../testdata/pwm_2motifs.rds")

test_that("ReadPWM works", {
  skip_if_not_installed("TFBSTools")
  result <- ReadPWM("../testdata/pwm_dir")
  expect_s4_class(result, "PWMatrixList")
  expect_equal(length(result), 2)
  # short_names = TRUE by default
  expect_true(all(c("AHR", "ALX1") %in% names(result)))
  # ID should remain the full motif ID
  expect_equal(TFBSTools::ID(result[["AHR"]]), "AHR.H14CORE.0.P.B")
  # check matrix dimensions (4 nucleotides x n positions)
  expect_equal(nrow(TFBSTools::Matrix(result[["AHR"]])), 4)
  expect_equal(ncol(TFBSTools::Matrix(result[["AHR"]])), 10)
  expect_equal(ncol(TFBSTools::Matrix(result[["ALX1"]])), 20)
})

test_that("ReadPWM short_names FALSE works", {
  skip_if_not_installed("TFBSTools")
  result <- ReadPWM("../testdata/pwm_dir", short_names = FALSE)
  expect_true(all(
    c("AHR.H14CORE.0.P.B", "ALX1.H14CORE.0.SM.B") %in% names(result)
  ))
})

test_that("ReadJASPAR works", {
  skip_if_not_installed("TFBSTools")
  result <- ReadJASPAR("../testdata/test_jaspar.txt")
  expect_s4_class(result, "PWMatrixList")
  expect_equal(length(result), 2)
  expect_true(all(c("MA0004", "MA0069") %in% names(result)))
  # check matrix dimensions
  expect_equal(nrow(TFBSTools::Matrix(result[["MA0004"]])), 4)
  expect_equal(ncol(TFBSTools::Matrix(result[["MA0004"]])), 6)
  expect_equal(ncol(TFBSTools::Matrix(result[["MA0069"]])), 14)
})

test_that("ReadJASPAR errors on malformed input", {
  skip_if_not_installed("TFBSTools")
  bad_file <- tempfile()
  writeLines(c(">BAD_MOTIF", "0.1 0.2 0.3 0.4", "0.4 0.3 0.2 0.1"), bad_file)
  expect_error(ReadJASPAR(bad_file), "expected 4")
  unlink(bad_file)
})
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)


test_that("AddMotifs works", {
  skip_on_cran()
  skip_if_not_installed("motifmatchr")
  expect_warning(motif <- AddMotifs(
    atac_small[["peaks"]],
    genome,
    pwm,
    verbose = FALSE
  ))
  expect_equal(dim(Motifs(motif)), c(100, 2))
})

test_that("AddMotifs works with fakechr", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  features <- suppressWarnings(c(
    granges(atac_small),
    GenomicRanges::GRanges(
      seqnames = "fake_chr",
      ranges = IRanges::IRanges(start = 1, end = 1000)
    )
  ))
  mat <- FeatureMatrix(
    fragments = fragments,
    features = features,
    keep_all_features = TRUE,
    fragtk = FALSE,
    verbose = FALSE
  )
  expect_warning(object <- CreateGRangesAssay(
    counts = mat,
    fragments = fragments,
    verbose = FALSE
  ))
  skip_on_cran()
  skip_if_not_installed("motifmatchr")
  motif <- suppressWarnings(AddMotifs(
    granges(object)[1:20],
    genome,
    pwm,
    verbose = FALSE
  ))
  expect_equal(dim(motif), c(20, 2))
  test <- SetAssayData(
    object = object,
    layer = "motifs",
    new.data = motif,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(test)), c(20, 2))
})

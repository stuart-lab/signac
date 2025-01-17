pwm <- readRDS("../testdata/pwm_2motifs.rds")
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)


test_that("AddMotifs works", {
  skip_if_not_installed("motifmatchr")
  motif <- AddMotifs(
    atac_small[["peaks"]],
    genome,
    pwm,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(motif)), c(323, 2))
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
    verbose = FALSE
  )
  expect_warning(object <- CreateChromatinAssay(
    counts = mat,
    sep = c(":", "-"),
    fragments = fragments,
    verbose = FALSE
  ))
  skip_if_not_installed("motifmatchr")
  motif <- suppressWarnings(AddMotifs(
    granges(object)[1:20],
    genome,
    pwm,
    verbose = FALSE
  ))
  expect_equal(dim(motif), c(20, 2))
  expect_warning(test <- SetAssayData(
    object = object,
    layer = "motifs",
    new.data = motif,
    verbose = FALSE
  ))
  expect_equal(dim(Motifs(test)), c(324, 2))
  # If some features are not in the ChromatinAssay
  # They should be removed
  expect_warning(test <- SetAssayData(
    object = atac_small,
    layer = "motifs",
    new.data = motif,
    verbose = FALSE
  ))
  expect_equal(dim(Motifs(test)), c(323, 2))
})

pwm <- readRDS("../testdata/pwm_2motifs.rds")
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)


test_that("AddMotifs works", {
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
    sep = c("-", "-"),
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
  test <- SetAssayData(
    object = object,
    layer = "motifs",
    new.data = motif,
    verbose = FALSE
  )
  expect_equal(dim(Motifs(test)), c(20, 2))
})

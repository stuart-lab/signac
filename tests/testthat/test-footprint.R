suppressWarnings(RNGversion(vstr = "3.5.3"))

fpath_headered <- system.file(
  "extdata", "fragments_header.tsv.gz", package = "Signac"
)
cells <- colnames(x = atac_small)
frags <- CreateFragmentObject(
  path = fpath_headered,
  cells = cells,
  verbose = FALSE,
  tolerance = 0.5
)
Fragments(atac_small) <- NULL
Fragments(atac_small) <- frags
pwm <- readRDS("../testdata/pwm_2motifs.rds")
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)

test_that("Footprint works", {
  skip_if_not_installed("motifmatchr")
  expect_warning(test_atac <- AddMotifs(
    object = atac_small,
    genome = genome,
    pfm = pwm,
    verbose = FALSE
  ))
  test_atac <- Footprint(
    object = test_atac,
    motif.name = names(x = pwm),
    compute.expected = FALSE,
    genome = genome
  )
  test_df <- GetFootprintData(object = test_atac, features = c("MA0031.1"))
  expect_equal(sum(test_df$count), 1013)

  test_width <- length(pwm[[1]]) +
    RegionAggr(test_atac)[[1]]@upstream +
    RegionAggr(test_atac)[[1]]@downstream
  expect_equal(dim(x = RegionAggr(object = test_atac)[[1]]@matrix),
               c(100, test_width))
  expect_equal(RegionAggNames(object = test_atac), c("MA0030.1", "MA0031.1"))
})

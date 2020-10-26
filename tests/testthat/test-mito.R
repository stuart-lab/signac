test_that("Importing data from mgatk works", {
  set.seed(1)
  data.dir <- system.file("extdata", "test_mgatk", package="Signac")
  mgatk <- ReadMGATK(dir = data.dir)

  # We imported 11 cells
  expect_equal(object = dim(mgatk$depth), expected = c(11,1))

  # We get + and - counts for every position
  expect_equal(object = dim(mgatk$counts), expected = c(8*dim(mgatk$refallele)[1],11))
})

test_that("Variant calling is operational", {
  data.dir <- system.file("extdata", "test_mgatk", package="Signac")
  mgatk <- ReadMGATK(dir = data.dir)
  var_df <- IdentifyVariants(mgatk$counts, mgatk$refallele)
  possible_vars_df <- var_df[var_df$strand_correlation > 0.65 & var_df$n_cells_conf_detected > 3,]
  homoplasmic_vars <- possible_vars_df[possible_vars_df$mean > 0.9,"variant"]
  heteroplasmic_vars <- possible_vars_df[possible_vars_df$mean < 0.5,"variant"]

  # We get the expected # of homoplasmic variants
  expect_equal(object = length(homoplasmic_vars), expected = 42)

  # We get the expected # of heteroplasmic variants
  expect_equal(object = length(heteroplasmic_vars), expected = 11)

})

test_that("Allele frequency calculation works", {
  data.dir <- system.file("extdata", "test_mgatk", package="Signac")
  mgatk <- ReadMGATK(dir = data.dir)
  vars.compute <- c("627G>A", "1888G>A", "1888G>C")
  expected <- new("dgCMatrix",
                  i = c(1L, 0L, 2L, 1L, 0L, 2L, 1L, 2L, 0L, 0L,1L, 1L, 0L, 2L),
                  p = c(0L, 1L, 3L, 3L, 4L, 6L, 8L, 9L, 10L, 11L, 12L, 14L),
                  Dim = c(3L, 11L), Dimnames = list(c("627G>A", "1888G>A", "1888G>C"),
                                                    c("ACAGGCCGTGGTCGAA-1", "AAACTGCAGAGTCCGA-1", "AAAGATGCAATGGCTT-1",
                                                      "ACTAGGTAGTGTCGGA-1", "ACCATCCTCTTGTCGC-1", "ACAGAAATCTATCTCA-1",
                                                      "AATGGCTCAACGTAGG-1", "ACCAAACGTTGTGAGG-1", "AACCGATTCAAGAGAT-1",
                                                      "ACTAGGTCAATGCCAT-1", "AAAGGGCAGCTACGCC-1")),
                  x = c(0.0947368421052632,
                        0.0955631399317406, 0.00273972602739726, 0.108108108108108, 0.0833333333333333,
                        0.00409836065573771, 0.0256410256410256, 0.00233100233100233,
                        0.087431693989071, 0.13768115942029, 0.0308788598574822, 0.072,
                        0.0254237288135593, 0.00689655172413793),
                  factors = list())
  alleles <- AlleleFreq(
    object = mgatk$counts, variants = vars.compute
  )
  expect_equal(object = alleles, expected = expected)
})

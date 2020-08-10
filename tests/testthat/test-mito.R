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

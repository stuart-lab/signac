suppressWarnings(RNGversion(vstr = "3.5.3"))

count_vec <- c(4, 3, 2, 3, 2, 2, 2, 8, 1, 1, 1, 4, 3, 1, 5, 2, 2, 2, 1, 1, 2,
               1, 1, 3, 1, 1, 2, 1, 3, 2, 1, 2, 1, 2, 2, 1, 3, 1, 2, 1, 1, 2,
               1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1)

mononucleosome <- c(1, 1, 0, 2, 0, 0, 1, 5, 0, 0, 1, 4, 1, 1, 4, 1, 1, 1, 0, 1, 1,
                    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1,
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)

nucleosome_free <- c(3, 1, 2, 0, 2, 2, 0, 2, 0, 0, 0, 0, 2, 0, 1, 1, 1, 1, 1, 0, 1,
                     0, 1, 2, 1, 1, 2, 1, 3, 2, 1, 2, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1,
                     0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1)

read_vec <- c(5, 4, 2, 5, 3, 2, 3, 13, 1, 2, 6, 8, 5, 1, 8, 4, 6, 3, 1, 1, 2,
              4, 1, 4, 1, 8, 12, 1, 4, 4, 2, 19, 4, 3, 3, 1, 4, 3, 4, 1, 1,
              7, 1, 1, 1, 1, 1, 3, 2, 1, 3, 5, 2, 1)

test_that("CountFragments works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  counts <- CountFragments(fragments = fpath)

  expect_equal(
    object = counts$frequency_count,
    expected = count_vec
  )
  expect_equal(
    object = counts$mononucleosomal,
    expected = mononucleosome
  )

  expect_equal(
    object = counts$nucleosome_free,
    expected = nucleosome_free
  )

  expect_equal(
    object = counts$reads_count,
    expected = read_vec
  )
})

test_that("ExtractFragments works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- paste0("test_", cells)
  frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
  counts <- ExtractFragments(fragments = frags, verbose = FALSE)

  expect_equal(
    object = counts$CB,
    expected = names(x = cells)
  )

  expect_equal(
    object = head(x = counts$frequency_count),
    expected = c(0, 0, 3, 8, 0, 0)
  )
  expect_equal(
    object = head(x = counts$mononucleosomal),
    expected = c(0, 0, 0, 5, 0, 0)
  )

  expect_equal(
    object = head(x = counts$nucleosome_free),
    expected = c(0, 0, 3, 2, 0, 0)
  )

  expect_equal(
    object = head(x = counts$reads_count),
    expected = c(0, 0, 4, 13, 0, 0)
  )
})

test_that("ValidateCells works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- paste0("test_", cells)
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  valid <- Signac:::ValidateCells(
    object = frags,
    verbose = FALSE,
    tolerance = 0.5
  )
  invalid <- Signac:::ValidateCells(
    object = frags,
    verbose = FALSE,
    tolerance = 0
  )
  expect_true(object = valid)
  expect_false(object = invalid)
})

test_that("ValidateHash works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- paste0("test_", cells)
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  valid <- Signac:::ValidateHash(
    object = frags,
    verbose = FALSE
  )
  expect_true(object = valid)
})

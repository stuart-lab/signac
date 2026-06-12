test_that("LoadGWAS reads valid file", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tbase_pair_location\tp_value\tvariant_id",
    "chr1\t1000\t1e-8\trs1",
    "chr1\t2000\t1e-4\trs2"
  ), tf)
  gwas <- LoadGWAS(gwas.file = tf)
  expect_s3_class(gwas, "data.frame")
  expect_equal(nrow(gwas), 2)
  expect_setequal(
    colnames(gwas),
    c("chromosome", "base_pair_location", "p_value", "variant_id")
  )
  expect_equal(gwas$chromosome, c("chr1", "chr1"))
  expect_equal(gwas$base_pair_location, c(1000, 2000))
  expect_equal(gwas$p_value, c(1e-8, 1e-4))
  expect_equal(gwas$variant_id, c("rs1", "rs2"))
})

test_that("LoadGWAS errors on missing required columns", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition",
    "chr1\t1000"
  ), tf)
  expect_error(LoadGWAS(gwas.file = tf), regexp = "Missing required")
})

test_that("LoadGWAS fills variant_id from effect/other alleles", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tbase_pair_location\tp_value\teffect_allele\tother_allele",
    "chr1\t1000\t1e-8\tA\tG"
  ), tf)
  gwas <- LoadGWAS(gwas.file = tf)
  expect_true("variant_id" %in% colnames(gwas))
  expect_equal(gwas[["variant_id"]][1], "chr1_1000_A_G")
})

test_that("LoadLDData reads valid file and renames position to base_pair_location", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tr2",
    "chr1\t1000\t0.9",
    "chr1\t2000\t0.5"
  ), tf)
  ld <- LoadLDData(ld.file = tf)
  expect_s3_class(ld, "data.frame")
  expect_setequal(colnames(ld), c("chromosome", "base_pair_location", "r2"))
  expect_equal(ld$r2, c(0.9, 0.5))
  expect_equal(ld$base_pair_location, c(1000L, 2000L))
  expect_type(ld$base_pair_location, "integer")
})

test_that("LoadLDData handles mixed-case headers", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "Chromosome\tPosition\tR2",
    "chr1\t1000\t0.9",
    "chr1\t2000\t0.5"
  ), tf)
  ld <- LoadLDData(ld.file = tf)
  expect_s3_class(ld, "data.frame")
  expect_equal(nrow(ld), 2)
  expect_setequal(colnames(ld), c("chromosome", "base_pair_location", "r2"))
  expect_equal(ld$chromosome, c("chr1", "chr1"))
  expect_equal(ld$base_pair_location, c(1000L, 2000L))
  expect_equal(ld$r2, c(0.9, 0.5))
})

test_that("LoadLDData errors on missing required columns", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c("chromosome\tposition", "chr1\t1000"), tf)
  expect_error(LoadLDData(ld.file = tf), regexp = "Missing required")
})

test_that("LoadCredibleSets reads valid file with PIP filter", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tpip\tcs",
    "chr1\t1000\t0.5\t1",
    "chr1\t2000\t0.001\t1",  # below threshold
    "chr1\t3000\t0.5\t-1"     # not in credible set
  ), tf)
  cs <- LoadCredibleSets(credset.file = tf, credset.threshold = 0.01)
  expect_s3_class(cs, "data.frame")
  # only chr1:1000 passes both PIP threshold and cs != -1
  expect_equal(nrow(cs), 1)
  expect_equal(cs$base_pair_location, 1000)
  expect_equal(cs$pip, 0.5)
  expect_equal(cs$credset_id, 1)
})

test_that("LoadCredibleSets handles mixed-case headers", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "Chromosome\tPosition\tPIP\tCS",
    "chr1\t1000\t0.5\t1",
    "chr1\t2000\t0.001\t1",  # below threshold
    "chr1\t3000\t0.5\t-1"     # not in credible set
  ), tf)
  cs <- LoadCredibleSets(credset.file = tf, credset.threshold = 0.01)
  expect_s3_class(cs, "data.frame")
  # only chr1:1000 passes both PIP threshold and cs != -1
  expect_equal(nrow(cs), 1)
  expect_setequal(
    colnames(cs),
    c("chromosome", "base_pair_location", "pip", "credset_id")
  )
  expect_equal(cs$chromosome, "chr1")
  expect_equal(cs$base_pair_location, 1000L)
  expect_equal(cs$pip, 0.5)
  expect_equal(cs$credset_id, 1)
})

test_that("LoadCredibleSets warns when nothing passes", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tpip\tcs",
    "chr1\t1000\t0.001\t1"
  ), tf)
  expect_warning(
    LoadCredibleSets(credset.file = tf, credset.threshold = 0.01),
    regexp = "No variants"
  )
})

test_that("LoadCredibleSets errors on missing required columns", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c("chromosome\tposition", "chr1\t1000"), tf)
  expect_error(LoadCredibleSets(credset.file = tf), regexp = "Missing required")
})

test_that("LoadCredibleSets errors on no header", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c("chr1\t1000\t0.5\t1"), tf)
  # data.table will read this with V1, V2 headers
  expect_error(LoadCredibleSets(credset.file = tf), regexp = "header|Missing")
})

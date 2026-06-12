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
    object = fragments,
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

# Build a ChromatinAssay5 with a synthetic motif matrix (no motifmatchr /
# TFBSTools dependency) so FindMotifs can be tested in CI.
make_motif_object <- function(hits) {
  # hits: named list of motif -> character vector of features with the motif
  obj <- atac_small
  peaks <- obj[["peaks"]]
  feats <- rownames(x = peaks)
  motif.ids <- names(x = hits)
  mat <- Matrix::Matrix(
    data = 0, nrow = length(x = feats), ncol = length(x = motif.ids),
    sparse = TRUE, dimnames = list(feats, motif.ids)
  )
  for (m in motif.ids) {
    mat[hits[[m]], m] <- 1
  }
  motif.obj <- CreateMotifObject(
    data = mat,
    pwm = setNames(as.list(motif.ids), motif.ids),
    motif.names = setNames(as.list(motif.ids), motif.ids)
  )
  peaks <- SetAssayData(
    object = peaks, layer = "motifs", new.data = motif.obj, verbose = FALSE
  )
  # meta features required by the automatic (matched) background selection
  peaks[["GC.percent"]] <- setNames(
    seq(40, 60, length.out = length(x = feats)), feats
  )
  peaks[["count"]] <- setNames(
    Matrix::rowSums(x = GetAssayData(object = peaks, layer = "counts")), feats
  )
  obj[["peaks"]] <- peaks
  obj
}

test_that("FindMotifs computes valid hypergeometric enrichment p-values", {
  feats <- rownames(x = atac_small[["peaks"]])
  query <- feats[1:20]
  bg <- setdiff(x = feats, y = query) # disjoint background of 80
  # ENRICHED: 15/20 query and 3/80 background; FLAT: spread evenly
  obj <- suppressWarnings(make_motif_object(hits = list(
    ENRICHED = c(query[1:15], bg[1:3]),
    FLAT = feats[seq(1, length(x = feats), by = 4)]
  )))

  res <- suppressWarnings(FindMotifs(
    object = obj, features = query, background = bg, verbose = FALSE
  ))

  # p-values are valid probabilities
  expect_true(all(res$pvalue >= 0 & res$pvalue <= 1))
  expect_false(any(is.nan(x = res$pvalue)))
  expect_gt(res["ENRICHED", "pvalue"], 0)

  # exact match to a one-sided Fisher's exact test on the 2x2 table
  # query: 15 hit / 5 miss ; background: 3 hit / 77 miss
  tab <- matrix(data = c(15, 3, 5, 77), nrow = 2)
  expect_equal(
    object = res["ENRICHED", "pvalue"],
    expected = fisher.test(x = tab, alternative = "greater")$p.value
  )

  # enriched motif is significant; flat motif is not
  expect_lt(res["ENRICHED", "pvalue"], 1e-5)
  expect_gt(res["FLAT", "pvalue"], 0.05)
  expect_gt(res["ENRICHED", "fold.enrichment"], res["FLAT", "fold.enrichment"])
})

test_that("FindMotifs p-value is invariant to background containing the query", {
  feats <- rownames(x = atac_small[["peaks"]])
  query <- feats[1:20]
  bg <- setdiff(x = feats, y = query)
  obj <- suppressWarnings(make_motif_object(hits = list(
    ENRICHED = c(query[1:15], bg[1:3]),
    FLAT = feats[seq(1, length(x = feats), by = 4)]
  )))

  # background disjoint from query vs. a superset that contains the query
  res.disjoint <- suppressWarnings(FindMotifs(
    object = obj, features = query, background = bg, verbose = FALSE
  ))
  res.superset <- suppressWarnings(FindMotifs(
    object = obj, features = query, background = feats, verbose = FALSE
  ))
  # the union population is identical, so the p-value must match
  expect_equal(
    object = res.disjoint["ENRICHED", "pvalue"],
    expected = res.superset["ENRICHED", "pvalue"]
  )
})

test_that("FindMotifs returns finite p-values when query exceeds background", {
  feats <- rownames(x = atac_small[["peaks"]])
  obj <- suppressWarnings(make_motif_object(hits = list(
    M1 = feats[1:40],
    M2 = feats[20:44]
  )))
  # query (70) larger than the disjoint background (30)
  query <- feats[1:70]
  bg <- feats[71:100]
  res <- suppressWarnings(FindMotifs(
    object = obj, features = query, background = bg, verbose = FALSE
  ))
  expect_false(any(is.nan(x = res$pvalue)))
  expect_true(all(res$pvalue >= 0 & res$pvalue <= 1))
})

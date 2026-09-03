library(GenomicRanges)
library(SeuratObject)
library(Matrix)

suppressWarnings(RNGversion(vstr = "3.5.3"))

# BinarizeCounts ---------------------------------------------------------------

test_that("BinarizeCounts works", {
  set.seed(1)
  mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
  bin_mat <- BinarizeCounts(object = mat)

  mat_sparse <- as(object = mat, Class = "CsparseMatrix")
  bin_mat_sparse <- BinarizeCounts(object = mat_sparse)

  expect_equal(object = as.vector(bin_mat[1, ]), expected = c(0, 1, 0, 1, 1))
  expect_equal(
    object = as.vector(bin_mat_sparse[1, ]), expected = c(0, 1, 0, 1, 1)
  )
})

test_that("BinarizeCounts.Assay binarizes counts", {
  m <- matrix(sample(0:5, 100, replace = TRUE), 10, 10)
  rownames(m) <- paste0("p", 1:10)
  colnames(m) <- paste0("c", 1:10)
  a <- SeuratObject::CreateAssayObject(counts = m)
  res <- BinarizeCounts(a, verbose = FALSE)
  expect_s4_class(res, "Assay")
  vals <- as.vector(GetAssayData(res, layer = "counts"))
  expect_setequal(unique(vals), c(0, 1))
  expect_equal(sum(vals == 1), sum(m > 0))
})

test_that("BinarizeCounts.Assay5 binarizes counts", {
  m <- matrix(sample(0:5, 100, replace = TRUE), 10, 10)
  rownames(m) <- paste0("p", 1:10)
  colnames(m) <- paste0("c", 1:10)
  m_sparse <- as(m, "CsparseMatrix")
  a <- SeuratObject::CreateAssay5Object(counts = m_sparse)
  res <- BinarizeCounts(a, verbose = FALSE)
  expect_s4_class(res, "Assay5")
  vals <- as.vector(GetAssayData(res, layer = "counts"))
  expect_setequal(unique(vals), c(0, 1))
  expect_equal(sum(vals == 1), sum(m > 0))
})

test_that("BinarizeCounts works on Seurat", {
  b <- BinarizeCounts(atac_small, verbose = FALSE)
  expect_s4_class(b, "Seurat")
  mat <- GetAssayData(b[["peaks"]], layer = "counts")
  expect_true(all(mat@x %in% c(0, 1)))
})

test_that("BinarizeCounts on dense matrix produces 0/1 values", {
  m <- matrix(sample(0:3, 100, replace = TRUE), 10, 10)
  res <- BinarizeCounts(m, verbose = FALSE)
  expect_true(all(res %in% c(0, 1)))
})

# DownsampleFeatures -----------------------------------------------------------

test_that("DownsampleFeatures works", {
  set.seed(1)
  atac_ds <- DownsampleFeatures(object = atac_small, n = 5, verbose = FALSE)
  expect_equal(
    object = SeuratObject::VariableFeatures(object = atac_ds),
    expected = c(
      "chr1:898350-899223",
      "chr1:955190-956101",
      "chr1:1059203-1060060",
      "chr1:1250624-1251529",
      "chr1:860139-860923"
    )
  )
})


# FindTopFeatures --------------------------------------------------------------

test_that("FindTopFeatures works", {
  VariableFeatures(atac_small) <- NULL
  atac_small <- FindTopFeatures(object = atac_small)
  expect_equal(
    object = head(SeuratObject::VariableFeatures(object = atac_small)),
    expected = c(
      "chr1:1115790-1116694",
      "chr1:1307720-1308738",
      "chr1:778263-779184",
      "chr1:1231645-1232553",
      "chr1:1012999-1013896",
      "chr1:1068591-1069593"
    )
  )
})

test_that("FindTopFeatures returns the most accessible features first", {
  obj <- atac_small
  VariableFeatures(obj) <- NULL
  obj <- FindTopFeatures(object = obj, min.cutoff = "q50", verbose = FALSE)
  counts <- Matrix::rowSums(
    x = LayerData(object = obj[["peaks"]], layer = "counts")
  )
  vf <- VariableFeatures(object = obj)
  # ordered by decreasing total count
  expect_false(is.unsorted(x = rev(x = counts[vf])))
  # and drawn from the upper half of the count distribution, not the lower
  expect_gt(median(x = counts[vf]), median(x = counts))
})

test_that("FindTopFeatures sums counts over layers", {
  obj <- atac_small[["peaks"]]
  counts <- LayerData(object = obj, layer = "counts")
  split.obj <- CreateAssay5Object(
    counts = list(counts.1 = counts[, 1:50], counts.2 = counts[, 51:100])
  )
  split.obj <- FindTopFeatures(object = split.obj, verbose = FALSE)
  meta <- split.obj[[]]
  # the summary is stored under every layer, and is the same for each
  expect_true("vf_topfeatures_counts.1_count" %in% colnames(x = meta))
  expect_true("vf_topfeatures_counts.2_count" %in% colnames(x = meta))
  expect_equal(
    object = meta$vf_topfeatures_counts.1_count,
    expected = meta$vf_topfeatures_counts.2_count
  )
  # counts are summed over the layers, not taken from one of them
  expect_equal(
    object = meta$vf_topfeatures_counts.1_count,
    expected = unname(obj = Matrix::rowSums(x = counts))
  )
})

test_that("FindTopFeatures.default works", {
  m <- as(matrix(rpois(100, lambda = 2), 10, 10), "CsparseMatrix")
  rownames(m) <- paste0("p", 1:10)
  res <- FindTopFeatures(m, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_true("count" %in% colnames(res))
})

test_that("FindTopFeatures.Assay5 populates VariableFeatures", {
  obj <- atac_small[["peaks"]]
  VariableFeatures(obj) <- NULL
  res <- FindTopFeatures(obj, verbose = FALSE)
  expect_s4_class(res, "Assay5")
  expect_gt(length(VariableFeatures(res)), 0)
  expect_true(all(VariableFeatures(res) %in% rownames(obj)))
})

test_that("FindTopFeatures Seurat with min.cutoff = q0 retains all features", {
  obj <- atac_small
  VariableFeatures(obj) <- NULL
  res <- FindTopFeatures(obj, min.cutoff = "q0", verbose = FALSE)
  # q0 means features above 0th percentile — all features
  expect_equal(length(VariableFeatures(res)), nrow(obj[["peaks"]]))
})

test_that("FindTopFeatures Seurat with min.cutoff = q50 retains ~half", {
  obj <- atac_small
  VariableFeatures(obj) <- NULL
  res <- FindTopFeatures(obj, min.cutoff = "q50", verbose = FALSE)
  # q50 means features above 50th percentile — roughly half
  vf <- length(VariableFeatures(res))
  expect_gt(vf, 0)
  expect_lt(vf, nrow(obj[["peaks"]]))
})

test_that("FindTopFeatures with numeric min.cutoff filters by count", {
  obj <- atac_small
  VariableFeatures(obj) <- NULL
  res <- FindTopFeatures(obj, min.cutoff = 5, verbose = FALSE)
  feature_meta <- res[["peaks"]][[]]
  count_col <- grep("count$", colnames(feature_meta), value = TRUE)[1]
  var_col <- grep("variable$", colnames(feature_meta), value = TRUE)[1]
  is_var <- feature_meta[[var_col]]
  # all variable features have count > 5; all non-variable have count <= 5
  expect_true(all(feature_meta[[count_col]][is_var] > 5))
  expect_true(all(feature_meta[[count_col]][!is_var] <= 5))
})

test_that("FindTopFeatures with NULL min.cutoff sets all features as variable", {
  obj <- atac_small
  VariableFeatures(obj) <- NULL
  res <- FindTopFeatures(obj, min.cutoff = NULL, verbose = FALSE)
  expect_equal(length(VariableFeatures(res)), nrow(obj[["peaks"]]))
})

test_that("FindTopFeatures Assay5 with numeric min.cutoff = 0 filters out zero-count", {
  obj <- atac_small[["peaks"]]
  res <- FindTopFeatures(obj, min.cutoff = 0, verbose = FALSE)
  feature_meta <- res[[]]
  count_col <- grep("count$", colnames(feature_meta), value = TRUE)[1]
  var_col <- grep("variable$", colnames(feature_meta), value = TRUE)[1]
  is_var <- feature_meta[[var_col]]
  expect_true(all(feature_meta[[count_col]][is_var] > 0))
})

# FRiP -------------------------------------------------------------------------

test_that("FRiP works", {
  # use a total that is genuinely larger than the in-peak counts; with
  # total.fragments = "nCount_peaks" the fraction is trivially 1 for every cell
  obj <- atac_small
  obj$total <- obj$nCount_peaks * 2
  obj <- FRiP(object = obj, assay = "peaks", total.fragments = "total")
  # FRiP is the in-peak count over the supplied total
  expect_equal(
    object = as.vector(x = obj$FRiP),
    expected = as.vector(x = obj$nCount_peaks / obj$total)
  )
  # with total = 2 * in-peak counts every (non-empty) cell is exactly 0.5
  frip <- obj$FRiP[!is.nan(x = obj$FRiP)]
  expect_true(all(abs(frip - 0.5) < 1e-8))
})

test_that("FRiP with custom col.name", {
  res <- FRiP(
    object = atac_small, assay = "peaks",
    total.fragments = "nCount_peaks", col.name = "myFRiP"
  )
  expect_true("myFRiP" %in% colnames(res[[]]))
})

test_that("FRiP with missing column produces NA values", {
  # Missing column results in NA FRiP values (not an error)
  res <- FRiP(
    object = atac_small, assay = "peaks", total.fragments = "nonexistent"
  )
  expect_s4_class(res, "Seurat")
  expect_true(all(is.na(res$FRiP)))
})

# PearsonResidualVar -----------------------------------------------------------

test_that("PearsonResidualVar calculates variance correctly", {
  set.seed(42)
  mat <- Matrix::rsparsematrix(nrow = 50, ncol = 30, density = 0.2)
  result <- PearsonResidualVar(mat, theta = 100, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_true("ResidualVariance" %in% colnames(result))
  expect_true("mean" %in% colnames(result))
  expect_true("count" %in% colnames(result))
  expect_true("rank" %in% colnames(result))
})

test_that("PearsonResidualVar batch processing accumulates correctly", {
  set.seed(123)
  mat <- Matrix::rsparsematrix(nrow = 20, ncol = 50, density = 0.15)
  result_batch5 <- PearsonResidualVar(mat, ncell.batch = 5, verbose = FALSE)
  result_batch25 <- PearsonResidualVar(mat, ncell.batch = 25, verbose = FALSE)
  result_batch50 <- PearsonResidualVar(mat, ncell.batch = 50, verbose = FALSE)
  expect_equal(result_batch5$ResidualVariance,
    result_batch25$ResidualVariance, tolerance = 1e-10
  )
  expect_equal(result_batch25$ResidualVariance,
    result_batch50$ResidualVariance, tolerance = 1e-10
  )
})

test_that("PearsonResidualVar handles edge cases", {
  mat <- Matrix::Matrix(0, nrow = 5, ncol = 10, sparse = TRUE)
  mat[1, c(1, 3, 5)] <- c(2, 4, 6)
  mat[3, c(2, 4)] <- c(1, 3)
  result <- PearsonResidualVar(mat, verbose = FALSE)
  expect_equal(nrow(result), 5)
  expect_equal(nrow(result[result$mean > 0, ]), 2)
  expect_equal(result$mean[1], 1.2, tolerance = 1e-10)
  expect_equal(result$mean[3], 0.4, tolerance = 1e-10)
})

test_that("PearsonResidualVar.default returns variance per feature", {
  set.seed(1)
  m <- as(matrix(rpois(1000, lambda = 5), 100, 10), "CsparseMatrix")
  rownames(m) <- paste0("p", 1:100)
  colnames(m) <- paste0("c", 1:10)
  res <- suppressWarnings(PearsonResidualVar(
    object = m, nfeatures = 50, min.counts = 1, verbose = FALSE
  ))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(m))
  expect_true(all(c("mean", "ResidualVariance", "count", "rank") %in%
    colnames(res)))
  expect_true(all(res$ResidualVariance >= 0))
})

test_that("PearsonResidualVar.Seurat stores ranks in feature metadata", {
  set.seed(1)
  res <- PearsonResidualVar(atac_small, min.counts = 1, verbose = FALSE)
  expect_s4_class(res, "Seurat")
  feature_meta <- res[["peaks"]][[]]
  prv_cols <- grep("pearson", colnames(feature_meta), value = TRUE)
  expect_true(length(prv_cols) >= 4)
  expect_true(any(grepl("ResidualVariance", prv_cols)))
  expect_true(any(grepl("rank$", prv_cols)))
})

test_that("PearsonResidualVar.StdAssay stores ranks", {
  set.seed(1)
  res <- PearsonResidualVar(
    atac_small[["peaks"]], min.counts = 1, verbose = FALSE
  )
  expect_s4_class(res, "Assay5")
  feature_meta <- res[[]]
  expect_true(any(grepl("ResidualVariance", colnames(feature_meta))))
})

# RunTFIDF ---------------------------------------------------------------------

test_that("RunTFIDF works on matrix with each method", {
  set.seed(42)
  m <- matrix(sample(0:5, 100, replace = TRUE), 10, 10)
  rownames(m) <- paste0("p", 1:10)
  colnames(m) <- paste0("c", 1:10)
  m_sparse <- as(m, "CsparseMatrix")
  results <- lapply(X = 1:4, FUN = function(meth) {
    as.matrix(RunTFIDF(object = m_sparse, method = meth, verbose = FALSE))
  })
  # the four methods produce four distinct normalizations of the same shape
  for (res in results) expect_equal(dim(res), dim(m_sparse))
  expect_equal(length(unique(lapply(X = results, FUN = round, 6))), 4L)
  # method 1 (default) must equal the documented log(TF * IDF) transform
  tf <- t(x = t(x = m) / colSums(x = m))
  idf <- ncol(x = m) / rowSums(x = m)
  expect_equal(results[[1]], log1p(x = tf * idf * 1e4), tolerance = 1e-6)
})

test_that("RunTFIDF accepts a precomputed IDF vector", {
  m <- matrix(sample(1:5, 100, replace = TRUE), 10, 10)
  rownames(m) <- paste0("p", 1:10)
  colnames(m) <- paste0("c", 1:10)
  m_sparse <- as(m, "CsparseMatrix")
  idf <- rep(2, 10)
  res <- RunTFIDF(object = m_sparse, idf = idf, verbose = FALSE)
  expect_equal(dim(res), dim(m_sparse))
})

test_that("RunTFIDF errors on bad IDF vector", {
  m <- matrix(sample(1:5, 100, replace = TRUE), 10, 10)
  m_sparse <- as(m, "CsparseMatrix")
  expect_error(
    RunTFIDF(m_sparse, idf = "wrong", verbose = FALSE), regexp = "numeric"
  )
  expect_error(
    RunTFIDF(m_sparse, idf = c(0, rep(1, 9)), verbose = FALSE),
    regexp = "cannot be zero"
  )
  expect_error(
    RunTFIDF(m_sparse, idf = c(1, 2), verbose = FALSE), regexp = "Length"
  )
})

test_that("RunTFIDF Seurat populates data layer with normalized values", {
  res <- RunTFIDF(atac_small, verbose = FALSE)
  expect_s4_class(res, "Seurat")
  d <- GetAssayData(res[["peaks"]], layer = "data")
  expect_equal(dim(d), dim(GetAssayData(atac_small[["peaks"]], layer = "counts")))
  # method 1 produces strictly positive log-scale values where counts > 0
  expect_true(all(d@x > 0))
})

test_that("RunTFIDF gives the same result for dense and sparse input", {
  m_sparse <- LayerData(object = atac_small[["peaks"]], layer = "counts")
  m_dense <- as.matrix(x = m_sparse)
  for (method in 1:4) {
    sparse.res <- suppressWarnings(
      RunTFIDF(object = m_sparse, method = method, verbose = FALSE)
    )
    dense.res <- suppressWarnings(
      RunTFIDF(object = m_dense, method = method, verbose = FALSE)
    )
    expect_equal(
      object = as.matrix(x = dense.res),
      expected = as.matrix(x = sparse.res)
    )
    # zero-count features and cells must not produce Inf or NaN
    expect_false(anyNA(x = as.matrix(x = dense.res)))
  }
})

test_that("RunTFIDF handles features and cells with no counts", {
  m <- LayerData(object = atac_small[["peaks"]], layer = "counts")
  m[1, ] <- 0
  m[, 1] <- 0
  res <- suppressWarnings(RunTFIDF(object = m, verbose = FALSE))
  expect_false(anyNA(x = as.matrix(x = res)))
  expect_true(all(as.matrix(x = res)[1, ] == 0))
  expect_true(all(as.matrix(x = res)[, 1] == 0))
})

test_that("RunTFIDF method 3 produces different values than method 1", {
  res1 <- RunTFIDF(atac_small, method = 1, verbose = FALSE)
  res3 <- RunTFIDF(atac_small, method = 3, verbose = FALSE)
  d1 <- GetAssayData(res1[["peaks"]], layer = "data")
  d3 <- GetAssayData(res3[["peaks"]], layer = "data")
  expect_equal(dim(d1), dim(d3))
  expect_false(isTRUE(all.equal(as.vector(d1), as.vector(d3))))
})

test_that("RunTFIDF.StdAssay produces normalized assay", {
  res <- RunTFIDF(atac_small[["peaks"]], verbose = FALSE)
  expect_s4_class(res, "Assay5")
  d <- GetAssayData(res, layer = "data")
  expect_equal(dim(d), dim(GetAssayData(atac_small[["peaks"]], layer = "counts")))
  expect_true(all(d@x > 0))
})

test_that("RunTFIDF default with empty cells warning", {
  m <- matrix(0, 5, 5)
  rownames(m) <- paste0("p", 1:5)
  colnames(m) <- paste0("c", 1:5)
  m_sparse <- as(m, "CsparseMatrix")
  expect_warning(
    RunTFIDF(m_sparse, verbose = FALSE), regexp = "0 total"
  )
})

# FitMeanVar -------------------------------------------------------------------

test_that("FitMeanVar fits successfully on a larger sparse matrix", {
  set.seed(1)
  mat <- as(matrix(rpois(10000, lambda = 3), 1000, 10), "CsparseMatrix")
  rownames(mat) <- paste0("p", 1:1000)
  res <- suppressWarnings(FitMeanVar(object = mat, verbose = FALSE))
  expect_s3_class(res, "data.frame")
})

test_that("FitMeanVar.data.frame validates input", {
  expect_error(
    FitMeanVar(object = data.frame(x = 1:10), verbose = FALSE),
    regexp = "Mean and variance"
  )
})

# RegionStats ------------------------------------------------------------------

test_that("RegionStats keeps stats aligned when a region is on an absent seqlevel", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("BSgenome")

  # tiny on-disk genome: chr1/chr2/chr3 (chrUn is deliberately NOT present)
  genome_seqs <- Biostrings::DNAStringSet(c(
    chr1 = paste(rep("A", 100), collapse = ""),       # GC% = 0
    chr2 = paste(rep("C", 100), collapse = ""),       # GC% = 100
    chr3 = paste(rep(c("G", "A"), 50), collapse = "") # GC% = 50
  ))
  fa <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(genome_seqs, filepath = fa)
  Rsamtools::indexFa(fa)
  on.exit(unlink(c(fa, paste0(fa, ".fai"))), add = TRUE)
  genome <- Rsamtools::FaFile(fa)

  # the absent-contig region (chrUn) sits *before* kept regions, so a positional
  # reorder would silently misalign the per-region stats
  gr <- GRanges(c("chr1", "chrUn", "chr2", "chr3"), IRanges(1, 100))
  GenomeInfoDb::seqlevels(gr) <- c("chr1", "chr2", "chr3", "chrUn")

  res <- suppressWarnings(RegionStats(object = gr, genome = genome))

  # rows must come back in the original region order
  expect_equal(rownames(res), as.character(1:4))
  # GC.percent must stay matched to the correct region (NA for the absent one)
  expect_equal(unname(res[, "GC.percent"]), c(0, NA, 100, 50))
  # the absent-seqlevel region is the only NA row
  expect_equal(which(is.na(res[, "GC.percent"])), 2L)
})

# ATACqc -----------------------------------------------------------------------

test_that("ATACqc errors when fragtk not installed", {
  if (nchar(Sys.which("fragtk")) == 0) {
    expect_error(
      ATACqc(
        object = "fake.tsv.gz",
        annotations = GRanges("chr1", IRanges(1, 10)),
        verbose = FALSE
      ),
      regexp = "fragtk"
    )
  } else {
    skip("fragtk present, skipping not-installed test")
  }
})

test_that("ATACqc does not overwrite existing metadata columns", {
  if (nchar(Sys.which("fragtk")) == 0) {
    skip("fragtk not installed")
  }
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  frags <- CreateFragmentObject(
    path = fpath,
    cells = colnames(x = atac_small),
    tolerance = 0.5,
    verbose = FALSE
  )
  obj <- atac_small
  Fragments(obj[["peaks"]]) <- frags
  # pre-existing column that ATACqc also computes, as would be left behind by
  # an earlier ATACqc run
  obj$TSS_enrichment <- 0.5

  res <- expect_warning(
    object = ATACqc(object = obj, assay = "peaks", verbose = FALSE),
    regexp = "TSS_enrichment"
  )
  # existing column is preserved, not overwritten
  expect_true(all(res$TSS_enrichment == 0.5))
  # the newly computed value is stored under the suffixed name
  expect_true("TSS_enrichment.atacqc" %in% colnames(res[[]]))

  # suffix = NULL restores overwriting behavior
  res2 <- suppressWarnings(
    ATACqc(object = obj, assay = "peaks", suffix = NULL, verbose = FALSE)
  )
  expect_false("TSS_enrichment.atacqc" %in% colnames(res2[[]]))
  expect_true("TSS_enrichment" %in% colnames(res2[[]]))
  expect_false(all(res2$TSS_enrichment == 0.5))
})

test_that("FeatureMatrix works", {
  computed_fmat <- readRDS("../testdata/featurematrix.rds")
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(
    path = fpath,
    cells = colnames(x = atac_small),
    tolerance = 0.5,
    verbose = FALSE
  )
  fm <- FeatureMatrix(
    object = fragments,
    features = granges(atac_small),
    fragtk = FALSE,
    verbose = FALSE
  )
  expect_identical(object = fm, expected = computed_fmat)
})

test_that("FitMeanVar does not select features below min.cutoff", {
  set.seed(42)
  cm <- as(matrix(rpois(400 * 80, 1.2), nrow = 400), "CsparseMatrix")
  rownames(cm) <- paste0("p", 1:400)
  colnames(cm) <- paste0("c", 1:80)

  # a threshold no feature can reach: nothing is eligible
  expect_error(
    suppressWarnings(FitMeanVar(
      cm, min.cutoff = 1e6, nfeatures = 400, bins = 20,
      sample_per_bin = 50, verbose = FALSE
    )),
    regexp = "No features remain after filtering by min.cutoff"
  )

  # a reachable threshold: only eligible features are marked variable
  thresh <- quantile(Matrix::rowSums(cm), 0.9)
  res <- suppressWarnings(FitMeanVar(
    cm, min.cutoff = thresh, nfeatures = 400, bins = 20,
    sample_per_bin = 50, verbose = FALSE
  ))
  expect_gt(sum(res$variable), 0)
  expect_equal(sum(res$variable & res$total.counts < thresh), 0)
  expect_equal(sum(res$variable), sum(res$total.counts >= thresh))
})

test_that("PearsonResidualVar honors min.counts", {
  set.seed(7)
  m <- as(matrix(rpois(100 * 40, 5), nrow = 100), "CsparseMatrix")
  rownames(m) <- paste0("p", 1:100)
  colnames(m) <- paste0("c", 1:40)
  counts <- Matrix::rowSums(m)
  thresh <- quantile(counts, 0.5)

  res <- PearsonResidualVar(
    object = m, min.counts = thresh, verbose = FALSE
  )
  # features below the threshold are not eligible for selection
  expect_true(all(is.na(res$rank[res$count < thresh])))
  expect_false(any(is.na(res$rank[res$count >= thresh])))

  # a permissive threshold leaves every feature eligible
  res.all <- PearsonResidualVar(object = m, min.counts = 1, verbose = FALSE)
  expect_false(anyNA(res.all$rank))
})

test_that("FindTopFeatures handles layers with different features", {
  counts <- LayerData(object = atac_small[["peaks"]], layer = "counts")
  # the two layers share only a subset of their features
  obj <- CreateAssay5Object(counts = list(
    counts.1 = counts[1:80, 1:50],
    counts.2 = counts[21:100, 51:100]
  ))
  res <- FindTopFeatures(object = obj, verbose = FALSE)
  meta <- res[[]]
  count.col <- grep("counts.1_count$", colnames(x = meta), value = TRUE)
  expect_equal(nrow(x = meta), 100)
  expect_false(anyNA(x = meta[[count.col]]))
  # counts are summed over both layers, treating an absent feature as zero
  expected <- Matrix::rowSums(x = counts[, 1:50])
  expected[81:100] <- 0
  expected <- expected + c(rep(0, 20), Matrix::rowSums(x = counts[21:100, 51:100]))
  expect_equal(
    object = unname(obj = meta[[count.col]]),
    expected = unname(obj = expected[rownames(x = meta)])
  )
})

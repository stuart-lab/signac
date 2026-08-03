library(GenomicRanges)
library(SeuratObject)
library(Matrix)

pwm <- readRDS("../testdata/pwm_2motifs.rds")
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)

# ReadPWM ----------------------------------------------------------------------

test_that("ReadPWM works", {
  skip_if_not_installed("TFBSTools")
  result <- ReadPWM("../testdata/pwm_dir")
  expect_s4_class(result, "PWMatrixList")
  expect_equal(length(result), 2)
  expect_true(all(c("AHR", "ALX1") %in% names(result)))
  expect_equal(TFBSTools::ID(result[["AHR"]]), "AHR.H14CORE.0.P.B")
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

# ReadJASPAR -------------------------------------------------------------------

test_that("ReadJASPAR works", {
  skip_if_not_installed("TFBSTools")
  result <- ReadJASPAR("../testdata/test_jaspar.txt")
  expect_s4_class(result, "PWMatrixList")
  expect_equal(length(result), 2)
  expect_true(all(c("MA0004", "MA0069") %in% names(result)))
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

# AddMotifs --------------------------------------------------------------------

test_that("AddMotifs works", {
  skip_on_cran()
  skip_if_not_installed("motifmatchr")
  expect_warning(motif <- AddMotifs(
    atac_small[["peaks"]], genome, pwm, verbose = FALSE
  ))
  expect_equal(dim(Motifs(motif)), c(100, 2))
})

test_that("AddMotifs works with fakechr", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  fragments <- CreateFragmentObject(fpath)
  features <- suppressWarnings(c(
    granges(atac_small),
    GRanges(seqnames = "fake_chr", ranges = IRanges(start = 1, end = 1000))
  ))
  mat <- FeatureMatrix(
    object = fragments, features = features,
    keep_all_features = TRUE, fragtk = FALSE, verbose = FALSE
  )
  expect_warning(object <- CreateGRangesAssay(
    counts = mat, fragments = fragments, verbose = FALSE
  ))
  skip_on_cran()
  skip_if_not_installed("motifmatchr")
  motif <- suppressWarnings(AddMotifs(
    granges(object)[1:20], genome, pwm, verbose = FALSE
  ))
  expect_equal(dim(motif), c(20, 2))
  test <- SetAssayData(
    object = object, layer = "motifs", new.data = motif, verbose = FALSE
  )
  expect_equal(dim(Motifs(test)), c(20, 2))
})

test_that("AddMotifs Assay/StdAssay errors", {
  expect_error(
    AddMotifs(
      object = SeuratObject::CreateAssayObject(
        counts = matrix(0, 3, 3,
          dimnames = list(c("a", "b", "c"), c("x", "y", "z")))
      ),
      genome = "fake",
      pfm = list()
    ),
    regexp = "standard Assay"
  )
  expect_error(
    AddMotifs(
      object = SeuratObject::CreateAssay5Object(
        counts = matrix(0, 3, 3,
          dimnames = list(c("a", "b", "c"), c("x", "y", "z")))
      ),
      genome = "fake",
      pfm = list()
    ),
    regexp = "Assay5"
  )
})

# CreateMotifMatrix ------------------------------------------------------------

test_that("CreateMotifMatrix works", {
  skip_on_cran()
  skip_if_not_installed("motifmatchr")
  motif.matrix <- suppressWarnings(CreateMotifMatrix(
    features = granges(atac_small), pwm = pwm, genome = genome
  ))
  expect_equal(dim(motif.matrix), c(100, 2))
})

# ConvertMotifID ---------------------------------------------------------------

test_that("ConvertMotifID converts via name <-> ID", {
  motif.names <- c(id1 = "TF1", id2 = "TF2", id3 = "TF3")
  expect_equal(ConvertMotifID(motif.names, name = "TF1"), "id1")
  expect_equal(
    unname(ConvertMotifID(motif.names, id = c("id1", "id2"))),
    c("TF1", "TF2")
  )
})

test_that("ConvertMotifID validates input", {
  motif.names <- c(id1 = "TF1")
  expect_error(
    ConvertMotifID(motif.names, name = "TF1", id = "id1"),
    regexp = "either name or ID"
  )
  expect_error(ConvertMotifID(motif.names), regexp = "Supply vector")
})

test_that("ConvertMotifID Motif object dispatches", {
  mat <- matrix(
    0, nrow = 4, ncol = 2,
    dimnames = list(c("p1", "p2", "p3", "p4"), c("M1", "M2"))
  )
  m <- CreateMotifObject(data = mat)
  m@motif.names <- list(M1 = "TF1", M2 = "TF2")
  expect_equal(ConvertMotifID(m, name = "TF1"), "M1")
})

test_that("ConvertMotifID Assay errors", {
  a <- SeuratObject::CreateAssayObject(
    counts = matrix(0, 3, 3,
      dimnames = list(c("a", "b", "c"), c("x", "y", "z")))
  )
  expect_error(ConvertMotifID(a), regexp = "standard Assay")
})

test_that("ConvertMotifID GRangesAssay with motifs works", {
  skip_if_not_installed("TFBSTools")
  obj <- atac_small[["peaks"]]
  nrow_data <- nrow(obj)
  pwm <- list(
    M1 = TFBSTools::PWMatrix(
      ID = "M1", name = "TF1",
      profileMatrix = matrix(rnorm(16), 4, 4,
        dimnames = list(c("A", "C", "G", "T"), NULL))
    )
  )
  motif_data <- as(
    matrix(0, nrow_data, 1, dimnames = list(rownames(obj), "M1")),
    "CsparseMatrix"
  )
  m <- CreateMotifObject(
    data = motif_data, pwm = pwm, motif.names = list(M1 = "TF1")
  )
  Motifs(obj) <- m
  expect_equal(unname(ConvertMotifID(obj, name = "TF1")), "M1")
})

# MotifPlot --------------------------------------------------------------------

test_that("MotifPlot errors on non-GRangesAssay", {
  # MotifPlot checks for ggseqlogo before checking the assay class, so the
  # GRangesAssay error is only reachable when ggseqlogo is installed
  skip_if_not_installed("ggseqlogo")
  obj <- atac_small
  obj[["normal"]] <- SeuratObject::CreateAssay5Object(
    counts = matrix(0, 3, ncol(atac_small),
      dimnames = list(c("a", "b", "c"), colnames(atac_small)))
  )
  expect_error(
    MotifPlot(obj, motifs = "M1", assay = "normal"),
    regexp = "GRangesAssay"
  )
})

test_that("MotifPlot gives informative error for missing motifs", {
  skip_if_not_installed("ggseqlogo")
  skip_if_not_installed("TFBSTools")
  ids <- c("MA0030.1", "MA0031.1")
  names(pwm) <- ids
  npeak <- nrow(atac_small[["peaks"]])
  set.seed(1)
  mat <- as(
    matrix(
      sample(c(0, 1), npeak * 2, replace = TRUE), ncol = 2,
      dimnames = list(rownames(atac_small[["peaks"]]), ids)
    ),
    "CsparseMatrix"
  )
  motif <- CreateMotifObject(data = mat, pwm = pwm)
  obj <- atac_small
  obj[["peaks"]] <- SetAssayData(
    object = obj[["peaks"]], layer = "motifs", new.data = motif
  )
  # a present motif still plots (by ID and by name)
  expect_s3_class(
    MotifPlot(obj, motifs = "MA0030.1", assay = "peaks"), "ggplot"
  )
  expect_s3_class(
    MotifPlot(obj, motifs = "FOXF2", assay = "peaks"), "ggplot"
  )
  # a partially-missing request warns and plots the motifs that were found
  expect_warning(
    p <- MotifPlot(obj, motifs = c("MA0030.1", "NOPE"), assay = "peaks"),
    regexp = "NOPE"
  )
  expect_s3_class(p, "ggplot")
  # an all-missing request errors, not the opaque ggseqlogo error
  expect_error(
    MotifPlot(obj, motifs = c("NOPE", "NOPE2"), assay = "peaks"),
    regexp = "None of the requested motifs"
  )
})

# GetMotifData / SetMotifData --------------------------------------------------

test_that("GetMotifData returns the underlying matrix", {
  set.seed(42)
  mat <- matrix(
    sample(c(0, 1), 30, replace = TRUE), nrow = 10,
    dimnames = list(paste0("p", 1:10), c("M1", "M2", "M3"))
  )
  m <- CreateMotifObject(data = mat)
  d <- GetMotifData(m, slot = "data")
  expect_s4_class(d, "CsparseMatrix")
  expect_equal(dim(d), dim(mat))
  expect_equal(rownames(d), rownames(mat))
  expect_equal(colnames(d), colnames(mat))
  expect_equal(as.vector(as.matrix(d)), as.vector(mat))
})

test_that("SetMotifData Seurat dispatches", {
  skip_if_not_installed("TFBSTools")
  obj <- atac_small
  nrow_data <- nrow(obj[["peaks"]])
  motif_data <- as(
    matrix(1, nrow_data, 1,
      dimnames = list(rownames(obj[["peaks"]]), "M1")),
    "CsparseMatrix"
  )
  m <- CreateMotifObject(data = motif_data)
  Motifs(obj) <- m
  new_data <- motif_data
  new_data[1, 1] <- 5
  obj <- SetMotifData(obj, slot = "data", new.data = new_data)
  expect_equal(GetMotifData(obj, slot = "data")[1, 1], 5)
})

# FindMotifs -------------------------------------------------------------------

test_that("FindMotifs errors on non-ChromatinAssay5", {
  obj <- atac_small
  obj[["normal"]] <- SeuratObject::CreateAssay5Object(
    counts = matrix(0, 3, ncol(atac_small),
      dimnames = list(c("a", "b", "c"), colnames(atac_small)))
  )
  DefaultAssay(obj) <- "normal"
  expect_error(
    FindMotifs(obj, features = c("a", "b"), verbose = FALSE),
    regexp = "Cannot run FindMotifs"
  )
})

test_that("FindMotifs returns enrichment data.frame for each motif", {
  skip_if_not_installed("TFBSTools")
  obj <- atac_small
  nrow_data <- nrow(obj[["peaks"]])
  set.seed(1)
  motif_data <- sparseMatrix(
    i = sample(1:nrow_data, 50, replace = TRUE),
    j = sample(1:5, 50, replace = TRUE),
    x = 1,
    dims = c(nrow_data, 5),
    dimnames = list(rownames(obj[["peaks"]]), paste0("M", 1:5))
  )
  motif_data <- (motif_data > 0) * 1
  motif_data <- as(motif_data, "CsparseMatrix")
  pwm <- lapply(1:5, function(i) {
    TFBSTools::PWMatrix(
      ID = paste0("M", i), name = paste0("TF", i),
      profileMatrix = matrix(rnorm(16), 4, 4,
        dimnames = list(c("A", "C", "G", "T"), NULL))
    )
  })
  names(pwm) <- paste0("M", 1:5)
  motif_obj <- CreateMotifObject(
    data = motif_data, pwm = pwm,
    motif.names = setNames(as.list(paste0("TF", 1:5)), paste0("M", 1:5))
  )
  Motifs(obj) <- motif_obj
  features <- head(rownames(obj[["peaks"]]), 20)
  background <- setdiff(rownames(obj[["peaks"]]), features)[1:30]
  res <- FindMotifs(
    object = obj, features = features,
    background = background, verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  # one row per motif in the matrix
  expect_equal(nrow(res), 5)
  # standard enrichment columns
  expect_true(all(c("motif", "fold.enrichment", "pvalue", "p.adjust") %in%
    colnames(res)))
  # p-values are well-formed probabilities
  expect_true(all(res$pvalue >= 0 & res$pvalue <= 1))
})

test_that("CreateMotifObject populates motif.names from a named PWM list", {
  motif.matrix <- matrix(
    data = sample(x = c(0, 1), size = 20, replace = TRUE), ncol = 2
  )
  rownames(x = motif.matrix) <- paste0("feature", 1:10)
  colnames(x = motif.matrix) <- c("MA1", "MA2")
  pwm <- list(MA1 = matrix(0, 4, 4), MA2 = matrix(0, 4, 4))

  mo <- CreateMotifObject(data = motif.matrix, pwm = pwm, motif.names = NULL)
  names.stored <- GetMotifData(object = mo, slot = "motif.names")
  expect_length(names.stored, 2)
  expect_equal(names(x = names.stored), c("MA1", "MA2"))
  expect_equal(unlist(x = names.stored, use.names = FALSE), c("MA1", "MA2"))

  # ConvertMotifID can now round-trip
  expect_equal(ConvertMotifID(object = mo, id = "MA1"), "MA1")

  # no pwm at all still yields an empty (but valid) slot
  mo.nopwm <- CreateMotifObject(data = motif.matrix)
  expect_length(GetMotifData(object = mo.nopwm, slot = "motif.names"), 0)
})

test_that("FindMotifs works when motif.names come from the PWM list", {
  # this combination previously failed with
  # "arguments imply differing number of rows: 3, 0"
  obj <- atac_small
  nrow_data <- nrow(obj[["peaks"]])
  set.seed(1)
  motif_data <- sparseMatrix(
    i = sample(1:nrow_data, 60, replace = TRUE),
    j = sample(1:3, 60, replace = TRUE),
    x = 1,
    dims = c(nrow_data, 3),
    dimnames = list(rownames(obj[["peaks"]]), paste0("M", 1:3))
  )
  motif_data <- as((motif_data > 0) * 1, "CsparseMatrix")
  pwm <- setNames(
    object = lapply(X = 1:3, FUN = function(i) matrix(0, 4, 4)),
    nm = paste0("M", 1:3)
  )
  Motifs(obj) <- CreateMotifObject(
    data = motif_data, pwm = pwm, motif.names = NULL
  )
  features <- head(rownames(obj[["peaks"]]), 20)
  background <- setdiff(rownames(obj[["peaks"]]), features)[1:40]
  res <- FindMotifs(
    object = obj, features = features, background = background, verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_equal(sort(res$motif.name), paste0("M", 1:3))
  expect_false(anyNA(res$motif.name))
})

test_that("Motif meta.data is validated and subset per motif", {
  motif.matrix <- matrix(
    data = sample(x = c(0, 1), size = 30, replace = TRUE), ncol = 3
  )
  rownames(x = motif.matrix) <- paste0("feature", 1:10)
  colnames(x = motif.matrix) <- paste0("M", 1:3)
  md <- data.frame(
    tf = c("TFA", "TFB", "TFC"), row.names = paste0("M", 1:3)
  )
  mo <- CreateMotifObject(data = motif.matrix, meta.data = md)
  expect_equal(nrow(GetMotifData(object = mo, slot = "meta.data")), 3)

  # subsetting motifs carries the matching metadata rows
  sub <- subset(x = mo, motifs = c("M1", "M3"))
  expect_equal(
    rownames(GetMotifData(object = sub, slot = "meta.data")), c("M1", "M3")
  )
  expect_equal(GetMotifData(object = sub, slot = "meta.data")$tf,
               c("TFA", "TFC"))

  # metadata that does not match the motifs is rejected
  bad <- data.frame(tf = "x", row.names = "not_a_motif")
  expect_error(
    CreateMotifObject(data = motif.matrix, meta.data = bad),
    regexp = "inconsistent"
  )
})

test_that("FindMotifs warns on missing query features", {
  skip_if_not_installed("TFBSTools")
  obj <- atac_small
  nrow_data <- nrow(obj[["peaks"]])
  motif_data <- as(
    matrix(1, nrow_data, 2,
      dimnames = list(rownames(obj[["peaks"]]), c("M1", "M2"))),
    "CsparseMatrix"
  )
  motif_obj <- CreateMotifObject(
    data = motif_data,
    pwm = list(
      M1 = TFBSTools::PWMatrix(ID = "M1", name = "TF1",
        profileMatrix = matrix(rnorm(16), 4, 4,
          dimnames = list(c("A", "C", "G", "T"), NULL))),
      M2 = TFBSTools::PWMatrix(ID = "M2", name = "TF2",
        profileMatrix = matrix(rnorm(16), 4, 4,
          dimnames = list(c("A", "C", "G", "T"), NULL)))
    ),
    motif.names = list(M1 = "TF1", M2 = "TF2")
  )
  Motifs(obj) <- motif_obj
  bg <- tail(rownames(obj[["peaks"]]), 30)
  features <- c(head(rownames(obj[["peaks"]]), 15), "ZZZNotFeature")
  expect_warning(
    FindMotifs(
      object = obj, features = features,
      background = bg, verbose = FALSE
    ),
    regexp = "not in the motif matrix"
  )
})

test_that("FindMotifs errors on zero query features", {
  skip_if_not_installed("TFBSTools")
  obj <- atac_small
  nrow_data <- nrow(obj[["peaks"]])
  motif_data <- as(
    matrix(1, nrow_data, 2,
      dimnames = list(rownames(obj[["peaks"]]), c("M1", "M2"))),
    "CsparseMatrix"
  )
  motif_obj <- CreateMotifObject(
    data = motif_data,
    pwm = list(
      M1 = TFBSTools::PWMatrix(ID = "M1", name = "TF1",
        profileMatrix = matrix(rnorm(16), 4, 4,
          dimnames = list(c("A", "C", "G", "T"), NULL))),
      M2 = TFBSTools::PWMatrix(ID = "M2", name = "TF2",
        profileMatrix = matrix(rnorm(16), 4, 4,
          dimnames = list(c("A", "C", "G", "T"), NULL)))
    ),
    motif.names = list(M1 = "TF1", M2 = "TF2")
  )
  Motifs(obj) <- motif_obj
  expect_error(
    suppressWarnings(FindMotifs(
      object = obj, features = "NotAFeature",
      background = head(rownames(obj[["peaks"]]), 30), verbose = FALSE
    )),
    regexp = "No query features"
  )
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

test_that("CreateMotifMatrix pads missing features with all-zero rows", {
  m <- Matrix::Matrix(
    data = c(1, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE, sparse = TRUE
  )
  rownames(m) <- c("chr1:1-100", "chr1:200-300")
  colnames(m) <- c("MA1", "MA2", "MA3")
  feature_order <- c("chr1:1-100", "chrZ:1-100", "chr1:200-300")

  out <- PadMissingFeatures(
    motif.matrix = m,
    feature_order = feature_order,
    missing.features = "chrZ:1-100"
  )

  expect_equal(dim(out), c(3, 3))
  # original feature order is restored
  expect_equal(rownames(out), feature_order)
  # the missing feature carries NO motif hits (the bug planted a spurious 1)
  expect_equal(sum(out["chrZ:1-100", ]), 0)
  # scored features are unchanged
  expect_equal(as.numeric(out["chr1:1-100", ]), c(1, 0, 1))
  expect_equal(as.numeric(out["chr1:200-300", ]), c(0, 1, 0))

  # works for (and preserves) a logical match matrix
  out.lgl <- PadMissingFeatures(
    motif.matrix = as(m > 0, "lgCMatrix"),
    feature_order = feature_order,
    missing.features = "chrZ:1-100"
  )
  expect_equal(sum(out.lgl["chrZ:1-100", ]), 0)
  expect_true(is(out.lgl, "lsparseMatrix"))
})

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

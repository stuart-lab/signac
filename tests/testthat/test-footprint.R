library(GenomicRanges)
library(SeuratObject)
library(Matrix)

suppressWarnings(RNGversion(vstr = "3.5.3"))

# Helper: build a valid RegionAggregation backed by a random matrix
build_fake_agg <- function(name = "TF1", upstream = 100L, downstream = 100L) {
  set.seed(1)
  ncells <- ncol(atac_small)
  total_pos <- upstream + downstream + 1
  mat <- matrix(
    runif(ncells * total_pos),
    nrow = ncells, ncol = total_pos
  )
  rownames(mat) <- colnames(atac_small)
  regions <- GRanges("chr1", IRanges(c(100, 1000), c(100, 1000)))
  CreateRegionAggregationObject(
    mat = mat, regions = regions,
    upstream = upstream, downstream = downstream,
    name = name, cells = rownames(mat)
  )
}

# Footprint end-to-end ---------------------------------------------------------

fpath_headered <- system.file(
  "extdata", "fragments_header.tsv.gz", package = "Signac"
)
cells <- colnames(x = atac_small)
frags <- CreateFragmentObject(
  path = fpath_headered, cells = cells, verbose = FALSE, tolerance = 0.5
)
Fragments(atac_small) <- NULL
Fragments(atac_small) <- frags
pwm <- readRDS("../testdata/pwm_2motifs.rds")
genome.fasta <- system.file("extdata", "chr1_start.fa", package = "Signac")
genome <- Rsamtools::FaFile(genome.fasta)

test_that("Footprint works end-to-end", {
  skip_if_not_installed("motifmatchr")
  expect_warning(test_atac <- AddMotifs(
    object = atac_small, genome = genome, pfm = pwm, verbose = FALSE
  ))
  test_atac <- Footprint(
    object = test_atac, motif.name = names(x = pwm),
    compute.expected = FALSE, genome = genome
  )
  test_df <- GetFootprintData(object = test_atac, features = c("MA0031.1"))
  expect_equal(sum(test_df$count), 1008, tolerance = 1e-3)

  test_width <- length(pwm[[1]]) +
    RegionAggr(test_atac)[[1]]@upstream +
    RegionAggr(test_atac)[[1]]@downstream
  expect_equal(
    dim(x = RegionAggr(object = test_atac)[[1]]@matrix),
    c(100, test_width)
  )
  expect_equal(RegionAggNames(object = test_atac), c("MA0030.1", "MA0031.1"))
})

# BackgroundMeanNorm -----------------------------------------------------------

test_that("BackgroundMeanNorm divides matrix by flank mean", {
  m <- matrix(1:1000, nrow = 5, ncol = 200)
  res <- Signac:::BackgroundMeanNorm(m, background = 50)
  expect_equal(dim(res), dim(m))
  flank_mean <- mean(m[, c(1:50, (ncol(m) - 49):ncol(m))])
  expect_equal(as.vector(res), as.vector(m) / flank_mean)
})

test_that("BackgroundMeanNorm errors on a window narrower than the flanks", {
  m <- matrix(1:400, nrow = 5, ncol = 80)
  expect_error(
    Signac:::BackgroundMeanNorm(m, background = 50),
    regexp = "Cannot compute flanking positions"
  )
})

# Cut matrix cell handling -----------------------------------------------------

test_that("MultiRegionCutMatrix restricts to the requested cells", {
  obj <- atac_small
  regions <- GenomicRanges::resize(
    granges(obj[["peaks"]])[1:5], width = 200, fix = "center"
  )
  want <- colnames(obj)[1:5]
  cm <- Signac:::MultiRegionCutMatrix(
    object = obj[["peaks"]], regions = regions, cells = want
  )
  expect_equal(nrow(cm), length(want))
  expect_setequal(rownames(cm), want)

  # NULL still means all cells in the fragment file
  cm.all <- Signac:::MultiRegionCutMatrix(
    object = obj[["peaks"]], regions = regions, cells = NULL
  )
  expect_equal(nrow(cm.all), ncol(obj))
})

test_that("Footprint labels aggregation rows with the matching cells", {
  skip_if_not_installed("motifmatchr")
  expect_warning(test_atac <- AddMotifs(
    object = atac_small, genome = genome, pfm = pwm, verbose = FALSE
  ))
  test_atac <- Footprint(
    object = test_atac, motif.name = names(x = pwm)[1],
    compute.expected = FALSE, genome = genome
  )
  agg <- RegionAggr(object = test_atac)[[1]]
  expect_equal(nrow(agg@matrix), length(agg@cells))
  expect_equal(agg@cells, Cells(test_atac[["peaks"]]))
  expect_false(anyNA(agg@cells))
})

test_that("Footprint requires a key for each set of supplied regions", {
  regions <- list(
    GRanges("chr1", IRanges(c(100, 200), width = 10)),
    GRanges("chr1", IRanges(c(300, 400), width = 10))
  )
  expect_error(
    Footprint(
      object = atac_small[["peaks"]], genome = genome, regions = regions,
      key = "onlyone", compute.expected = FALSE, verbose = FALSE
    ),
    regexp = "key needs to be supplied for each"
  )
})

# GetMotifSize -----------------------------------------------------------------

test_that("GetMotifSize returns motif width from RegionAggregation", {
  obj <- atac_small
  agg <- build_fake_agg(name = "TF1")
  RegionAggr(obj[["peaks"]], overwrite = TRUE) <- list(TF1 = agg)
  sizes <- Signac:::GetMotifSize(obj, features = "TF1")
  expect_equal(sizes[["TF1"]], 1L)
})

test_that("GetMotifSize errors when feature missing", {
  obj <- atac_small
  expect_error(
    Signac:::GetMotifSize(obj, features = "NoSuchTF"),
    regexp = "No footprinting data"
  )
})

test_that("GetMotifSize handles multiple features", {
  obj <- atac_small
  agg1 <- build_fake_agg(name = "TF1")
  agg2 <- build_fake_agg(name = "TF2")
  RegionAggr(obj[["peaks"]], overwrite = TRUE) <- list(agg1, agg2)
  sizes <- Signac:::GetMotifSize(obj, features = c("TF1", "TF2"))
  expect_equal(length(sizes), 2)
  expect_equal(names(sizes), c("TF1", "TF2"))
  # build_fake_agg uses 1-bp regions, so motif size should be 1
  expect_equal(unname(sizes), c(1L, 1L))
})

# GetFootprintData -------------------------------------------------------------

test_that("GetFootprintData errors on non-ChromatinAssay5", {
  obj <- atac_small
  obj[["normal"]] <- SeuratObject::CreateAssay5Object(
    counts = matrix(0, 3, ncol(atac_small),
      dimnames = list(c("a", "b", "c"), colnames(atac_small)))
  )
  expect_error(
    GetFootprintData(obj, features = "TF1", assay = "normal"),
    regexp = "ChromatinAssay5"
  )
})

test_that("GetFootprintData warns on missing feature", {
  obj <- atac_small
  Idents(obj) <- "cluster"
  expect_warning(
    res <- GetFootprintData(obj, features = "NoSuchTF"),
    regexp = "not found"
  )
})

test_that("GetFootprintData returns expected columns and class types", {
  obj <- atac_small
  Idents(obj) <- "cluster"
  agg <- build_fake_agg(name = "TF1")
  RegionAggr(obj[["peaks"]], overwrite = TRUE) <- list(TF1 = agg)
  res <- GetFootprintData(obj, features = "TF1")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("feature", "class", "position") %in% colnames(res)))
  expect_equal(unique(res$feature), "TF1")
  expect_setequal(unique(res$class), c("Observed", "Expected"))
  expect_gt(nrow(res), 0)
})

# GetFootprintRegions ----------------------------------------------------------

test_that("GetFootprintRegions finds motif positions by ID and common name", {
  nrow_data <- nrow(atac_small[["peaks"]])
  pos_list <- GRangesList(
    M1 = GRanges("chr1", IRanges(c(100, 200), c(110, 210))),
    M2 = GRanges("chr1", IRanges(300, 310))
  )
  motif_data <- as(
    matrix(0, nrow_data, 2,
      dimnames = list(rownames(atac_small[["peaks"]]), c("M1", "M2"))),
    "CsparseMatrix"
  )
  m <- new(
    Class = "Motif",
    data = motif_data, pwm = list(),
    motif.names = list(M1 = "TF1", M2 = "TF2"),
    meta.data = data.frame(row.names = c("M1", "M2")),
    positions = pos_list
  )
  # By motif ID
  by_id <- Signac:::GetFootprintRegions(m, motif.name = "M1")
  expect_s4_class(by_id, "GRanges")
  expect_equal(length(by_id), 2)
  expect_equal(start(by_id), c(100, 200))
  expect_equal(end(by_id), c(110, 210))
  # By common name should match the same positions
  by_name <- Signac:::GetFootprintRegions(m, motif.name = "TF1")
  expect_identical(by_name, by_id)
})

test_that("GetFootprintRegions errors on missing", {
  motif_data <- as(
    matrix(0, 2, 1, dimnames = list(c("a", "b"), "M1")),
    "CsparseMatrix"
  )
  m <- new(
    Class = "Motif",
    data = motif_data, pwm = list(),
    motif.names = list(M1 = "TF1"),
    meta.data = data.frame(row.names = "M1"),
    positions = GRangesList(M1 = GRanges("chr1", IRanges(1, 10)))
  )
  expect_error(
    Signac:::GetFootprintRegions(m, motif.name = "ZZZ"),
    regexp = "not found"
  )
})

test_that("GetFootprintRegions errors when positions missing", {
  motif_data <- as(
    matrix(0, 2, 1, dimnames = list(c("a", "b"), "M1")),
    "CsparseMatrix"
  )
  m <- CreateMotifObject(data = motif_data)
  expect_error(
    Signac:::GetFootprintRegions(m, motif.name = "M1"),
    regexp = "positions"
  )
})

# Footprint dispatch validation ------------------------------------------------

test_that("Footprint.Seurat dispatches and errors without motif/regions", {
  obj <- atac_small
  expect_error(
    Footprint(object = obj, genome = "fake"),
    regexp = "motif|region"
  )
})

test_that("Footprint.ChromatinAssay5 errors on motif and regions both", {
  obj <- atac_small[["peaks"]]
  expect_error(
    Footprint(
      object = obj, genome = "fake",
      motif.name = "TF1", regions = GRanges("chr1", IRanges(1, 1))
    ),
    regexp = "Choose one"
  )
})

test_that("Footprint.ChromatinAssay5 errors on key length mismatch", {
  obj <- atac_small[["peaks"]]
  nrow_data <- nrow(obj)
  pos_list <- GRangesList(
    M1 = GRanges("chr1", IRanges(c(100, 200), c(110, 210))),
    M2 = GRanges("chr1", IRanges(c(300), c(310)))
  )
  motif_data <- as(
    matrix(0, nrow_data, 2,
      dimnames = list(rownames(obj), c("M1", "M2"))),
    "CsparseMatrix"
  )
  m <- new(
    Class = "Motif",
    data = motif_data, pwm = list(),
    motif.names = list(M1 = "TF1", M2 = "TF2"),
    meta.data = data.frame(row.names = c("M1", "M2")),
    positions = pos_list
  )
  Motifs(obj) <- m
  expect_error(
    Footprint(
      object = obj, genome = "fake",
      motif.name = c("M1", "M2"), key = "single_key"
    ),
    regexp = "Key needs to be supplied for each"
  )
})

test_that("Footprint.ChromatinAssay5 errors on inconsistent widths", {
  obj <- atac_small[["peaks"]]
  expect_error(
    Footprint(
      object = obj, genome = "fake",
      regions = list(GRanges("chr1", IRanges(c(1, 100), c(50, 150)))),
      key = "test"
    ),
    regexp = "same width"
  )
})

# RegionAggr<- -----------------------------------------------------------------

test_that("RegionAggr<- assignment works", {
  obj <- atac_small
  agg <- build_fake_agg(name = "TF1")
  RegionAggr(obj[["peaks"]], overwrite = TRUE) <- list(TF1 = agg)
  expect_equal(length(RegionAggr(obj[["peaks"]])), 1)
  expect_equal(RegionAggNames(obj[["peaks"]]), "TF1")
  RegionAggr(obj, overwrite = TRUE) <- list(TF2 = build_fake_agg(name = "TF2"))
  expect_true("TF2" %in% RegionAggNames(obj))
})


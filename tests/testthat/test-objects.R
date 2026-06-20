library(GenomicRanges)
library(SeuratObject)
library(Matrix)

# Helper: build a Fragment2 object pointing at the test fragments file
make_fragments <- function() {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  CreateFragmentObject(
    path = fpath, cells = colnames(atac_small),
    tolerance = 0.5, verbose = FALSE
  )
}

# Create constructors -----------------------------------------------------------

test_that("CreateFragmentObject errors when fragment file missing", {
  expect_error(
    CreateFragmentObject(path = "/nonexistent/file.tsv.gz", verbose = FALSE),
    regexp = "does not exist"
  )
})

test_that("CreateFragmentObject populates the expected slots", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  frags <- CreateFragmentObject(path = fpath, verbose = FALSE)
  expect_s4_class(frags, "Fragment2")
  expect_equal(GetFragmentData(frags, "file.path"), normalizePath(fpath))
  expect_equal(
    GetFragmentData(frags, "file.index"),
    normalizePath(paste0(fpath, ".tbi"))
  )
})

test_that("CreateMotifObject stores the supplied data matrix", {
  mat <- matrix(
    sample(c(0, 1), 30, replace = TRUE), ncol = 3,
    dimnames = list(paste0("p", 1:10), c("M1", "M2", "M3"))
  )
  m <- CreateMotifObject(data = mat)
  expect_s4_class(m, "Motif")
  stored <- GetMotifData(m, slot = "data")
  expect_equal(dim(stored), dim(mat))
  expect_equal(as.matrix(stored), mat)
})

test_that("CreateMotifObject validates data type", {
  expect_error(
    CreateMotifObject(data = list(1, 2, 3)),
    regexp = "matrix or sparse"
  )
})

test_that("CreateMotifObject errors when row names missing", {
  mat <- matrix(0, nrow = 4, ncol = 2)
  colnames(mat) <- c("M1", "M2")
  expect_error(CreateMotifObject(data = mat), regexp = "Row names")
})

test_that("CreateMotifObject errors when PWM list is unnamed", {
  mat <- matrix(
    0, nrow = 2, ncol = 2,
    dimnames = list(c("p1", "p2"), c("M1", "M2"))
  )
  expect_error(
    CreateMotifObject(data = mat, pwm = list(matrix(0, 4, 4), matrix(0, 4, 4))),
    regexp = "named list"
  )
})

test_that("CreateMotifObject errors when names mismatch data", {
  mat <- matrix(
    0, nrow = 2, ncol = 2,
    dimnames = list(c("p1", "p2"), c("M1", "M2"))
  )
  pwm <- list(M3 = matrix(0, 4, 4), M2 = matrix(0, 4, 4))
  expect_error(CreateMotifObject(data = mat, pwm = pwm), regexp = "inconsistent")
})

test_that("CreateGRangesAssay warns on underscores in feature names", {
  starts <- seq(1, 100000, by = 10000)
  ends <- starts + 1000
  fnames <- paste0("chr1_", starts, "_", ends)
  counts <- matrix(
    sample(c(0, 1), 100, replace = TRUE),
    nrow = 10,
    dimnames = list(fnames, paste0("cell", 1:10))
  )
  # the warning is the behaviour under test; the subsequent GRanges parse
  # failure (dotted names are not chr:start-end) is incidental, so swallow it
  expect_warning(
    try(CreateGRangesAssay(counts = counts, verbose = FALSE), silent = TRUE),
    regexp = "underscores"
  )
})

test_that("CreateGRangesAssay attaches the supplied ranges to the assay", {
  starts <- seq(1, 10000, by = 1000)
  ends <- starts + 500
  fnames <- paste0("chr1:", starts, "-", ends)
  m <- matrix(
    sample(0:3, 100, replace = TRUE), 10, 10,
    dimnames = list(fnames, paste0("c", 1:10))
  )
  gr <- GRanges("chr1", IRanges(starts, ends))
  res <- CreateGRangesAssay(counts = m, ranges = gr, verbose = FALSE)
  expect_s4_class(res, "GRangesAssay")
  expect_equal(length(granges(res)), 10)
  expect_equal(start(granges(res)), starts)
  expect_equal(end(granges(res)), ends)
})

test_that("CreateChromatinAssay5 with annotation validates", {
  m <- matrix(
    sample(0:3, 100, replace = TRUE), 10, 10,
    dimnames = list(
      paste0("chr1-", 1:10 * 1000, "-", 1:10 * 1000 + 500),
      paste0("c", 1:10)
    )
  )
  expect_error(
    CreateChromatinAssay5(counts = m, annotation = "not granges", verbose = FALSE),
    regexp = "GRanges"
  )
  gr <- GRanges("chr1", IRanges(1, 100))
  expect_error(
    CreateChromatinAssay5(counts = m, annotation = gr, verbose = FALSE),
    regexp = "tx_id"
  )
})

# CreateRegionAggregationObject ------------------------------------------------

test_that("CreateRegionAggregationObject works", {
  set.seed(1)
  mat <- matrix(0, nrow = 5, ncol = 11)
  rownames(mat) <- paste0("cell", 1:5)
  regions <- GRanges("chr1", IRanges(c(100, 200, 300), c(100, 200, 300)))
  agg <- CreateRegionAggregationObject(
    mat = mat, regions = regions,
    upstream = 5, downstream = 5, name = "TF1"
  )
  expect_s4_class(agg, "RegionAggregation")
})

test_that("CreateRegionAggregationObject validates matrix dimensions", {
  expect_error(
    CreateRegionAggregationObject(
      mat = matrix(0, 0, 0),
      regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = 5, name = "TF"
    ),
    regexp = "zero dimensions"
  )
})

test_that("CreateRegionAggregationObject validates type", {
  expect_error(
    CreateRegionAggregationObject(
      mat = list(),
      regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = 5, name = "TF"
    ),
    regexp = "matrix or sparse"
  )
})

test_that("CreateRegionAggregationObject validates regions type", {
  expect_error(
    CreateRegionAggregationObject(
      mat = matrix(0, 5, 11, dimnames = list(paste0("c", 1:5), NULL)),
      regions = "not a GRanges",
      upstream = 5, downstream = 5, name = "TF"
    ),
    regexp = "GRanges object"
  )
})

test_that("CreateRegionAggregationObject validates upstream/downstream", {
  mat <- matrix(0, 5, 11, dimnames = list(paste0("c", 1:5), NULL))
  expect_error(
    CreateRegionAggregationObject(
      mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
      upstream = "bad", downstream = 5, name = "TF"
    ),
    regexp = "upstream"
  )
  expect_error(
    CreateRegionAggregationObject(
      mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = "bad", name = "TF"
    ),
    regexp = "downstream"
  )
})

test_that("CreateRegionAggregationObject validates expected vector", {
  mat <- matrix(0, 5, 11, dimnames = list(paste0("c", 1:5), NULL))
  expect_error(
    CreateRegionAggregationObject(
      mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = 5, name = "TF",
      expected = c(1, 2, 3)
    ),
    regexp = "match the number of positions"
  )
  expect_error(
    CreateRegionAggregationObject(
      mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = 5, name = "TF",
      expected = "bad"
    ),
    regexp = "numeric vector"
  )
})

test_that("CreateRegionAggregationObject validates cells", {
  mat <- matrix(0, 5, 11)
  expect_error(
    CreateRegionAggregationObject(
      mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = 5, name = "TF",
      cells = c("a", "b")
    ),
    regexp = "Number of cells"
  )
  expect_error(
    CreateRegionAggregationObject(
      mat = matrix(0, 5, 11),
      regions = GRanges("chr1", IRanges(1, 1)),
      upstream = 5, downstream = 5, name = "TF"
    ),
    regexp = "cells information not provided"
  )
})

test_that("IsCompatibleRegionAggregation works", {
  mat <- matrix(0, 5, 11, dimnames = list(paste0("c", 1:5), NULL))
  agg1 <- CreateRegionAggregationObject(
    mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
    upstream = 5L, downstream = 5L, name = "TF"
  )
  agg2 <- CreateRegionAggregationObject(
    mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
    upstream = 5L, downstream = 5L, name = "TF"
  )
  expect_true(Signac:::IsCompatibleRegionAggregation(agg1, agg2))
  mat2 <- matrix(0, 5, 16, dimnames = list(paste0("c", 1:5), NULL))
  agg3 <- CreateRegionAggregationObject(
    mat = mat2, regions = GRanges("chr1", IRanges(1, 1)),
    upstream = 10L, downstream = 5L, name = "TF"
  )
  expect_false(Signac:::IsCompatibleRegionAggregation(agg1, agg3))
})

test_that("SetAssayData merges same-name RegionAggregations without duplicating overlapping cells", {
  regions <- GRanges("chr1", IRanges(100, 100))
  make_agg <- function(cells, marker) {
    mat <- matrix(marker, nrow = length(cells), ncol = 11,
                  dimnames = list(cells, NULL))
    CreateRegionAggregationObject(
      mat = mat, regions = regions,
      upstream = 5L, downstream = 5L, name = "TF1"
    )
  }
  # agg1 has c1,c2,c3 ; agg2 shares c2,c3 and adds c4,c5
  agg1 <- make_agg(c("c1", "c2", "c3"), marker = c(1, 2, 3))
  agg2 <- make_agg(c("c2", "c3", "c4", "c5"), marker = c(20, 30, 40, 50))

  obj <- atac_small[["peaks"]]
  obj <- SetAssayData(obj, layer = "region.aggregation", new.data = agg1)
  obj <- suppressWarnings(
    SetAssayData(obj, layer = "region.aggregation", new.data = agg2)
  )

  res <- GetAssayData(object = obj, layer = "region.aggregation")[[1]]

  # overlapping cells (c2, c3) must not be duplicated
  expect_false(any(duplicated(res@cells)))
  expect_equal(sort(res@cells), c("c1", "c2", "c3", "c4", "c5"))
  # matrix rows must stay in sync with the cells vector
  expect_equal(nrow(res@matrix), length(res@cells))
  expect_no_error(validObject(res))
  # overlapping cells keep agg1's values (not recomputed); new cells come from agg2
  vals <- setNames(res@matrix[, 1], res@cells)
  expect_equal(vals[c("c1", "c2", "c3", "c4", "c5")],
               c(c1 = 1, c2 = 2, c3 = 3, c4 = 40, c5 = 50))
})

# show methods -----------------------------------------------------------------

test_that("show methods print the expected class banners", {
  expect_output(show(atac_small[["peaks"]]), regexp = "GRangesAssay")
  expect_output(show(make_fragments()), regexp = "Fragment v2")
  expect_output(
    show(as(atac_small[["peaks"]], "ChromatinAssay5")),
    regexp = "ChromatinAssay"
  )
  m <- CreateMotifObject(
    data = matrix(0, 4, 2,
      dimnames = list(c("p1", "p2", "p3", "p4"), c("M1", "M2")))
  )
  expect_output(show(m), regexp = "Motif")
})

# Accessor functions -----------------------------------------------------------

test_that("RegionAggr is empty by default on atac_small", {
  expect_equal(length(RegionAggr(atac_small)), 0)
  expect_equal(RegionAggNames(atac_small), character(0))
})

test_that("GetFragmentData returns the correct fragment file path", {
  frags <- make_fragments()
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  expect_equal(GetFragmentData(frags, "file.path"), normalizePath(fpath))
  expect_true(file.exists(GetFragmentData(frags, "file.index")))
})

test_that("Cells<-.Fragment2 rejects unnamed vectors", {
  frags <- make_fragments()
  expect_error(
    Cells(frags) <- c("BC1", "BC2"),
    regexp = "named vector"
  )
})

# GetAssayData / SetAssayData --------------------------------------------------

test_that("GetAssayData GRangesAssay returns layer matrices and granges", {
  obj <- atac_small[["peaks"]]
  m <- GetAssayData(obj, layer = "counts")
  expect_equal(dim(m), c(nrow(obj), ncol(atac_small)))
  expect_equal(rownames(m), rownames(obj))
  # ranges slot is the GRanges stored on the assay
  expect_identical(GetAssayData(obj, layer = "ranges"), granges(obj))
})

test_that("GetAssayData GRangesAssay errors on bad layer", {
  expect_error(
    GetAssayData(atac_small[["peaks"]], layer = "nonexistent_slot"),
    regexp = "layer must be one of"
  )
})

test_that("SetAssayData fragments accepts a single Fragment2 object", {
  frags <- make_fragments()
  obj <- atac_small[["peaks"]]
  obj <- SetAssayData(obj, layer = "fragments", new.data = frags)
  res <- Fragments(obj)
  expect_equal(length(res), 1)
  expect_identical(
    GetFragmentData(res[[1]], "file.path"),
    GetFragmentData(frags, "file.path")
  )
})

test_that("SetAssayData ChromatinAssay5 fragments validates list", {
  obj <- atac_small[["peaks"]]
  expect_error(
    SetAssayData(obj, layer = "fragments", new.data = list("not_a_fragment")),
    regexp = "Fragment"
  )
})

test_that("SetAssayData ChromatinAssay5 annotation validates", {
  obj <- atac_small[["peaks"]]
  expect_error(
    SetAssayData(obj, layer = "annotation", new.data = "not granges"),
    regexp = "GRanges"
  )
  gr <- GRanges("chr1", IRanges(1, 100))
  expect_error(
    SetAssayData(obj, layer = "annotation", new.data = gr),
    regexp = "tx_id"
  )
  gr$tx_id <- "tx1"
  expect_error(
    SetAssayData(obj, layer = "annotation", new.data = gr),
    regexp = "gene_name"
  )
})

test_that("SetAssayData ChromatinAssay5 annotation with transcript_id works", {
  obj <- atac_small[["peaks"]]
  gr <- GRanges("chr1", IRanges(1, 100))
  gr$transcript_id <- "tx1"
  gr$gene_name <- "g1"
  gr$gene_id <- "g1"
  gr$gene_biotype <- "protein_coding"
  gr$type <- "exon"
  obj <- SetAssayData(obj, layer = "annotation", new.data = gr)
  expect_s4_class(Annotation(obj), "GRanges")
  expect_true("tx_id" %in% colnames(mcols(Annotation(obj))))
})

test_that("SetAssayData ChromatinAssay5 bias validates", {
  obj <- atac_small[["peaks"]]
  # a non-vector triggers the setter's own check (a list is a vector, so it
  # would slip past it and only fail later in slot validation)
  expect_error(
    SetAssayData(
      obj, layer = "bias", new.data = GenomicRanges::GRanges("chr1:1-2")
    ),
    regexp = "Bias must be provided as a vector"
  )
})

test_that("SetAssayData ChromatinAssay5 motifs validates and clears", {
  obj <- atac_small[["peaks"]]
  expect_error(
    SetAssayData(obj, layer = "motifs", new.data = "not motif"),
    regexp = "Motif"
  )
  obj <- SetAssayData(obj, layer = "motifs", new.data = NULL)
  expect_null(Motifs(obj))
})

test_that("SetAssayData ChromatinAssay5 links validates", {
  obj <- atac_small[["peaks"]]
  expect_error(
    SetAssayData(obj, layer = "links", new.data = "not list"),
    regexp = "GInteractions"
  )
  expect_error(
    SetAssayData(obj, layer = "links", new.data = list("not_GI")),
    regexp = "GInteractions"
  )
  obj <- SetAssayData(obj, layer = "links", new.data = NULL)
  expect_equal(length(Links(obj)), 0)
})

test_that("SetAssayData ChromatinAssay5 links with single GI and key", {
  obj <- atac_small[["peaks"]]
  gi <- InteractionSet::GInteractions(
    GRanges("chr1", IRanges(c(100, 200), c(150, 250))),
    GRanges("chr1", IRanges(c(500, 600), c(550, 650)))
  )
  obj <- SetAssayData(obj, layer = "links", new.data = gi, key = "mylink")
  expect_true("mylink" %in% names(Links(obj)))
  expect_error(
    SetAssayData(obj, layer = "links", new.data = gi),
    regexp = "Key"
  )
})

test_that("SetAssayData GRangesAssay errors on bad ranges length", {
  obj <- atac_small[["peaks"]]
  expect_error(
    SetAssayData(
      obj, layer = "ranges",
      new.data = GRanges("chr1", IRanges(1, 100))
    ),
    regexp = "Number of ranges"
  )
})

test_that("SetAssayData GRangesAssay rejects non-GRanges for ranges", {
  obj <- atac_small[["peaks"]]
  expect_error(
    SetAssayData(obj, layer = "ranges", new.data = "not granges"),
    regexp = "GRanges"
  )
})

# Setter functions -------------------------------------------------------------

test_that("Annotation<- sets and updates annotation", {
  # Annotation<- round-trip on both Seurat and assay
  obj <- atac_small
  orig <- Annotation(obj)
  Annotation(obj) <- orig
  expect_equal(length(Annotation(obj)), length(orig))
  expect_equal(Annotation(obj)$gene_name, orig$gene_name)

  pk <- atac_small[["peaks"]]
  Annotation(pk) <- orig
  expect_equal(length(Annotation(pk)), length(orig))
})

test_that("Fragments<- sets the same fragments on Seurat and assay, and NULL clears", {
  frags <- make_fragments()
  expected_path <- GetFragmentData(frags, "file.path")

  obj <- atac_small
  Fragments(obj) <- frags
  expect_equal(length(Fragments(obj)), 1)
  expect_identical(
    GetFragmentData(Fragments(obj)[[1]], "file.path"), expected_path
  )

  pk <- atac_small[["peaks"]]
  Fragments(pk) <- frags
  expect_equal(length(Fragments(pk)), 1)
  expect_identical(
    GetFragmentData(Fragments(pk)[[1]], "file.path"), expected_path
  )

  Fragments(obj) <- NULL
  expect_equal(length(Fragments(obj)), 0)
})

test_that("Motifs<- via Seurat method stores Motif object", {
  obj <- atac_small
  mat <- matrix(
    0, nrow(obj[["peaks"]]), 2,
    dimnames = list(rownames(obj[["peaks"]]), c("M1", "M2"))
  )
  mo <- CreateMotifObject(data = mat)
  Motifs(obj) <- mo
  res <- Motifs(obj)
  expect_s4_class(res, "Motif")
  expect_equal(dim(GetMotifData(res, slot = "data")),
               c(nrow(obj[["peaks"]]), 2))
  expect_equal(colnames(GetMotifData(res, slot = "data")), c("M1", "M2"))
})

test_that("Motifs<- with NULL clears motifs", {
  obj <- atac_small[["peaks"]]
  Motifs(obj) <- NULL
  expect_null(Motifs(obj))
})

test_that("Links<- stores GInteractions under key and NULL clears them", {
  gi <- InteractionSet::GInteractions(
    GRanges("chr1", IRanges(c(100, 200), c(150, 250))),
    GRanges("chr1", IRanges(c(500, 600), c(550, 650)))
  )
  obj <- atac_small
  Links(obj) <- list(test = gi)
  expect_equal(names(Links(obj)), "test")
  expect_equal(length(Links(obj)[["test"]]), 2)

  pk <- atac_small[["peaks"]]
  Links(pk) <- list(other = gi)
  expect_equal(names(Links(pk)), "other")

  Links(pk) <- NULL
  expect_equal(length(Links(pk)), 0)
})

test_that("RegionAggr<- on Seurat stores the aggregation object", {
  obj <- atac_small
  agg <- CreateRegionAggregationObject(
    mat = matrix(0, ncol(obj), 11,
      dimnames = list(colnames(obj), NULL)),
    regions = GRanges("chr1", IRanges(1, 1)),
    upstream = 5L, downstream = 5L, name = "TF"
  )
  RegionAggr(obj) <- list(TF = agg)
  res <- RegionAggr(obj)
  expect_equal(length(res), 1)
  expect_s4_class(res[[1]], "RegionAggregation")
  expect_equal(res[[1]]@name, "TF")
  expect_equal(RegionAggNames(obj), "TF")
})

# SetMotifData / GetMotifData --------------------------------------------------

test_that("SetMotifData works on Motif", {
  mat <- matrix(
    0, nrow = 4, ncol = 2,
    dimnames = list(c("p1", "p2", "p3", "p4"), c("M1", "M2"))
  )
  m <- CreateMotifObject(data = mat)
  new_data <- mat
  new_data[1, 1] <- 1
  m2 <- SetMotifData(m, slot = "data", new.data = new_data)
  expect_equal(GetMotifData(m2, slot = "data")[1, 1], 1)
})

# subset and merge -------------------------------------------------------------

test_that("subset.GRangesAssay subsets features and cells by character", {
  obj <- atac_small[["peaks"]]
  feats <- head(rownames(obj), 5)
  cells <- head(colnames(obj), 10)
  s <- subset(obj, features = feats, cells = cells)
  expect_equal(rownames(s), feats)
  expect_equal(colnames(s), cells)
  # granges slot is subset in lockstep with the features
  expect_equal(length(granges(s)), 5)
  expect_equal(as.character(granges(s)), as.character(granges(obj)[1:5]))
})

test_that("subset.GRangesAssay accepts logical and Rle selections", {
  obj <- atac_small[["peaks"]]
  log_feats <- c(rep(TRUE, 5), rep(FALSE, nrow(obj) - 5))
  expect_equal(rownames(subset(obj, features = log_feats)),
               head(rownames(obj), 5))
  rle_feats <- Rle(head(rownames(obj), 5))
  expect_equal(rownames(subset(obj, features = rle_feats)),
               head(rownames(obj), 5))
})

test_that("subset.ChromatinAssay5 errors on wrong-length logical selectors", {
  obj <- atac_small[["peaks"]]
  expect_error(
    subset(obj, cells = c(TRUE, FALSE)), regexp = "Incorrect number"
  )
  expect_error(
    suppressWarnings(subset(obj, features = c(TRUE, FALSE))),
    regexp = "Incorrect number"
  )
})

test_that("subset.Fragment2 retains only requested cells", {
  frags <- make_fragments()
  keep <- head(names(GetFragmentData(frags, "cells")), 10)
  s <- subset(frags, cells = keep)
  expect_s4_class(s, "Fragment2")
  expect_equal(length(GetFragmentData(s, "cells")), 10)
  expect_setequal(names(GetFragmentData(s, "cells")), keep)
})

test_that("subset.Fragment2 returns NULL when no cells match", {
  frags <- make_fragments()
  expect_null(subset(frags, cells = "NotACell"))
})

test_that("subset.Motif works", {
  mat <- matrix(
    sample(c(0, 1), 30, replace = TRUE), nrow = 10,
    dimnames = list(paste0("p", 1:10), c("M1", "M2", "M3"))
  )
  m <- CreateMotifObject(data = mat)
  s <- subset(m, features = c("p1", "p2"), motifs = c("M1"))
  expect_s4_class(s, "Motif")
  expect_equal(dim(GetMotifData(s, slot = "data")), c(2, 1))
})

test_that("subset.RegionAggregation keeps requested cells and rows", {
  mat <- matrix(seq_len(5 * 11), 5, 11)
  rownames(mat) <- paste0("c", 1:5)
  agg <- CreateRegionAggregationObject(
    mat = mat, regions = GRanges("chr1", IRanges(1, 1)),
    upstream = 5, downstream = 5, name = "TF"
  )
  s <- subset(agg, cells = c("c1", "c2"))
  expect_equal(slot(s, "cells"), c("c1", "c2"))
  expect_equal(nrow(slot(s, "matrix")), 2)
  expect_equal(unname(slot(s, "matrix")[, 1]), unname(mat[1:2, 1]))
  # cells = NULL returns the same object
  expect_identical(subset(agg, cells = NULL), agg)
  # No matching cells returns NULL
  expect_null(subset(agg, cells = "x"))
})

test_that("merge.GRangesAssay merges two GRangesAssays", {
  obj <- atac_small[["peaks"]]
  obj2 <- RenameCells(atac_small[["peaks"]],
                      new.names = paste0("x_", colnames(atac_small)))
  m <- suppressWarnings(merge(obj, obj2))
  expect_s4_class(m, "GRangesAssay")
  expect_equal(ncol(m), 2 * ncol(obj))
  expect_equal(nrow(m), nrow(obj))
  expect_equal(length(granges(m)), nrow(obj))
  # cells from both objects present
  expect_true(all(colnames(obj) %in% colnames(m)))
  expect_true(all(paste0("x_", colnames(obj)) %in% colnames(m)))
})

test_that("merge.ChromatinAssay5 merges two Seurat objects", {
  obj1 <- atac_small
  obj2 <- RenameCells(atac_small,
                      new.names = paste0("x_", colnames(atac_small)))
  m <- suppressWarnings(merge(obj1, obj2))
  expect_s4_class(m, "Seurat")
  expect_equal(ncol(m), 2 * ncol(atac_small))
  expect_true(all(colnames(atac_small) %in% colnames(m)))
})

test_that("merge.ChromatinAssay5 with non-Chromatin object falls back to Assay5", {
  obj1 <- as(atac_small[["peaks"]], "ChromatinAssay5")
  m <- matrix(0, 10, 10, dimnames = list(paste0("p", 1:10), paste0("c", 1:10)))
  obj2 <- SeuratObject::CreateAssay5Object(counts = m)
  res <- suppressWarnings(merge(obj1, obj2))
  # Without a second ChromatinAssay5 the merge returns a plain Assay5
  expect_s4_class(res, "Assay5")
  expect_false(inherits(res, "ChromatinAssay5"))
})

# Conversion / class transformation --------------------------------------------

test_that("as.Fragment2 carries path/hash/cells from old Fragment object", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  old <- new(
    Class = "Fragment",
    path = fpath, hash = c("a", "b"), cells = c(x = "BC1", y = "BC2")
  )
  new_obj <- as.Fragment2(old)
  expect_s4_class(new_obj, "Fragment2")
  expect_equal(GetFragmentData(new_obj, "file.path"), fpath)
  expect_equal(GetFragmentData(new_obj, "cells"), c(x = "BC1", y = "BC2"))
})

test_that("as.GRangesAssay attaches provided ranges to a plain Assay5", {
  starts <- seq(1, 10000, length.out = 10)
  ends <- starts + 500
  feat <- paste0("chr1-", starts, "-", ends)
  m <- matrix(
    sample(0:3, 100, replace = TRUE), 10, 10,
    dimnames = list(feat, paste0("c", 1:10))
  )
  a <- SeuratObject::CreateAssay5Object(counts = m)
  gr <- GRanges("chr1", IRanges(starts, ends))
  names(gr) <- feat
  res <- as.GRangesAssay(a, ranges = gr)
  expect_s4_class(res, "GRangesAssay")
  expect_equal(length(granges(res)), 10)
  expect_equal(start(granges(res)), starts)
})

test_that("as.GRangesAssay.ChromatinAssay5 errors on length mismatch", {
  obj <- atac_small[["peaks"]]
  gr <- granges(obj)[1:5]
  expect_error(
    as.GRangesAssay(as(obj, "ChromatinAssay5"), ranges = gr),
    regexp = "match number"
  )
})

test_that("as.ChromatinAssay5.Assay5 converts and preserves the counts", {
  starts <- seq(1, 10000, length.out = 10)
  ends <- starts + 500
  m <- matrix(
    sample(0:3, 100, replace = TRUE), 10, 10,
    dimnames = list(
      paste0("chr1-", starts, "-", ends),
      paste0("c", 1:10)
    )
  )
  a <- SeuratObject::CreateAssay5Object(counts = m)
  res <- as.ChromatinAssay5(a)
  expect_s4_class(res, "ChromatinAssay5")
  # conversion preserves dimensions, names and the count values
  expect_equal(dim(res), dim(a))
  expect_equal(rownames(res), rownames(a))
  expect_equal(
    as.matrix(GetAssayData(object = res, layer = "counts")),
    as.matrix(GetAssayData(object = a, layer = "counts"))
  )
})

# RenameCells methods ----------------------------------------------------------

test_that("RenameCells renames every cell on ChromatinAssay5 and RegionAggregation", {
  obj <- atac_small[["peaks"]]
  new_names <- paste0("renamed_", colnames(obj))
  r <- RenameCells(obj, new.names = new_names)
  expect_equal(colnames(r), new_names)

  agg <- CreateRegionAggregationObject(
    mat = matrix(0, 5, 11, dimnames = list(paste0("c", 1:5), NULL)),
    regions = GRanges("chr1", IRanges(1, 1)),
    upstream = 5L, downstream = 5L, name = "TF"
  )
  agg_new <- RenameCells(agg, new.names = paste0("new_c", 1:5))
  expect_equal(slot(agg_new, "cells"), paste0("new_c", 1:5))
})

# Validation helpers -----------------------------------------------------------

test_that("ValidateHash and ValidateFragments accept a freshly-built Fragment2", {
  frags <- make_fragments()
  expect_true(ValidateHash(frags, verbose = FALSE))
  expect_true(ValidateFragments(frags, verbose = FALSE))
})

# UpdateChromatinObject --------------------------------------------------------

test_that("UpdateChromatinObject errors on missing assay", {
  expect_error(
    UpdateChromatinObject(atac_small, chromatin.assay = "nonexistent"),
    regexp = "not in the object"
  )
})

test_that("UpdateChromatinObject with expression assay carries cells through", {
  res <- UpdateChromatinObject(
    object = atac_small,
    chromatin.assay = "peaks", expression.assay = "RNA",
    features = head(rownames(atac_small[["RNA"]]), 3)
  )
  expect_s4_class(res, "Seurat")
  # both assays preserved
  expect_true(all(c("peaks", "RNA") %in% Assays(res)))
  expect_equal(ncol(res), ncol(atac_small))
})

# AddFragments helper ----------------------------------------------------------

test_that("AddFragments adds fragments to assay", {
  frags <- make_fragments()
  obj <- atac_small[["peaks"]]
  res <- Signac:::AddFragments(obj, fragments = frags)
  expect_equal(length(Fragments(res)), 1)
})

library(GenomicRanges)
library(SeuratObject)
library(ggplot2)

# Shared setup: attach fragments and set cluster idents
setup_obj <- function() {
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  frags <- CreateFragmentObject(
    path = fpath, cells = colnames(atac_small),
    tolerance = 0.5, verbose = FALSE
  )
  obj <- atac_small
  Fragments(obj) <- frags
  Idents(obj) <- "cluster"
  obj
}

# DepthCor ---------------------------------------------------------------------

test_that("DepthCor plots n components when n is set", {
  p <- DepthCor(object = atac_small, n = 5)
  expect_s3_class(p, "ggplot")
  # one row (point) per component plotted
  expect_equal(nrow(p$data), 5)
  expect_equal(p$data$Component, seq_len(5))
})

test_that("DepthCor with n = NULL uses all reduction components", {
  n_components <- ncol(Embeddings(object = atac_small[["lsi"]]))
  expect_equal(n_components, 50)
  p <- DepthCor(object = atac_small, n = NULL)
  expect_s3_class(p, "ggplot")
  # NULL means plot every component in the reduction
  expect_equal(nrow(p$data), n_components)
  expect_equal(p$data$Component, seq_len(n_components))
})

test_that("DepthCor uses the default n when no n argument is supplied", {
  p <- DepthCor(object = atac_small)
  expect_s3_class(p, "ggplot")
  # default n = 10
  expect_equal(nrow(p$data), 10)
  expect_equal(p$data$Component, seq_len(10))
})

# DensityScatter ---------------------------------------------------------------

test_that("DensityScatter returns a ggplot", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("fields")
  atac_small$x_test <- runif(ncol(atac_small))
  atac_small$y_test <- runif(ncol(atac_small))
  p <- DensityScatter(object = atac_small, x = "x_test", y = "y_test")
  expect_s3_class(p, "ggplot")
})

test_that("DensityScatter validates input columns", {
  expect_error(
    DensityScatter(object = atac_small, x = "missing", y = "nCount_peaks"),
    regexp = "not found"
  )
  expect_error(
    DensityScatter(object = atac_small, x = "nCount_peaks", y = "missing"),
    regexp = "not found"
  )
})

test_that("DensityScatter handles log axes and quantiles", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("fields")
  p <- DensityScatter(
    object = atac_small,
    x = "nCount_peaks", y = "nFeature_peaks",
    log_x = TRUE, log_y = TRUE,
    quantiles = c(10, 90)
  )
  expect_s3_class(p, "ggplot")
})

test_that("DensityScatter handles TRUE quantiles", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("fields")
  p <- DensityScatter(
    object = atac_small,
    x = "nCount_peaks", y = "nFeature_peaks",
    quantiles = TRUE
  )
  expect_s3_class(p, "ggplot")
})

test_that("DensityScatter with raster = FALSE", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("fields")
  atac_small$xtest <- runif(ncol(atac_small))
  atac_small$ytest <- runif(ncol(atac_small))
  p <- DensityScatter(
    object = atac_small, x = "xtest", y = "ytest", raster = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("DensityScatter warns on bad quantile values", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("fields")
  atac_small$xtest <- runif(ncol(atac_small))
  atac_small$ytest <- runif(ncol(atac_small))
  expect_warning(
    DensityScatter(
      object = atac_small, x = "xtest", y = "ytest",
      quantiles = c(110, 200)
    ),
    regexp = "between 0 and 100"
  )
})

# PeakPlot ---------------------------------------------------------------------

test_that("PeakPlot returns a ggplot with and without peaks in region", {
  p <- PeakPlot(atac_small, region = "chr1:710000-780000")
  expect_s3_class(p, "ggplot")
  p2 <- PeakPlot(atac_small, region = "chr1:1-100")
  expect_s3_class(p2, "ggplot")
})

test_that("PeakPlot accepts a custom peaks GRanges", {
  pk <- GRanges("chr1", IRanges(start = c(750000, 760000), end = c(750500, 760500)))
  p <- PeakPlot(atac_small, region = "chr1:710000-780000", peaks = pk)
  expect_s3_class(p, "ggplot")
})

test_that("PeakPlot warns and continues with bad group.by", {
  expect_warning(
    PeakPlot(
      atac_small,
      region = "chr1:710000-780000",
      group.by = "nonexistent_col"
    ),
    regexp = "grouping"
  )
})

test_that("PeakPlot with custom color", {
  obj <- setup_obj()
  p <- PeakPlot(object = obj, region = "chr1:710000-780000", color = "red")
  expect_s3_class(p, "ggplot")
})

test_that("PeakPlot returns plot when no peaks in region", {
  p <- PeakPlot(atac_small, region = "chr2:1-100")
  expect_s3_class(p, "ggplot")
})

test_that("PeakPlot group.by works with custom peaks", {
  pk <- GRanges("chr1", IRanges(c(750000, 760000), c(750500, 760500)))
  pk$category <- c("A", "B")
  p <- PeakPlot(
    atac_small,
    region = "chr1:710000-780000",
    peaks = pk,
    group.by = "category"
  )
  expect_s3_class(p, "ggplot")
})

test_that("PeakPlot extends region", {
  p <- PeakPlot(
    atac_small,
    region = "chr1:780000-781000",
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  expect_s3_class(p, "ggplot")
})

# AnnotationPlot ---------------------------------------------------------------

test_that("AnnotationPlot returns a ggplot", {
  p <- AnnotationPlot(object = atac_small, region = "chr1:780000-790000")
  expect_s3_class(p, "ggplot")
})

test_that("AnnotationPlot transcript mode works", {
  p <- AnnotationPlot(
    object = atac_small,
    region = "chr1:780000-790000",
    mode = "transcript"
  )
  expect_s3_class(p, "ggplot")
})

test_that("AnnotationPlot rejects bad mode", {
  expect_error(
    AnnotationPlot(
      object = atac_small,
      region = "chr1:780000-790000",
      mode = "wrong"
    ),
    regexp = "Unknown mode"
  )
})

test_that("AnnotationPlot with empty annotation returns empty plot", {
  obj <- atac_small
  empty <- GRanges()
  mcols(empty)$tx_id <- character(0)
  mcols(empty)$gene_name <- character(0)
  mcols(empty)$gene_biotype <- character(0)
  mcols(empty)$gene_id <- character(0)
  mcols(empty)$type <- character(0)
  Annotation(obj) <- empty
  expect_s3_class(
    AnnotationPlot(object = obj, region = "chr1:780000-790000"),
    "ggplot"
  )
})

test_that("AnnotationPlot with empty region returns empty plot", {
  p <- AnnotationPlot(object = atac_small, region = "chr1:1-100")
  expect_s3_class(p, "ggplot")
})

test_that("AnnotationPlot uses extend args", {
  p <- AnnotationPlot(
    object = atac_small,
    region = "chr1:780000-790000",
    extend.upstream = 1000,
    extend.downstream = 1000
  )
  expect_s3_class(p, "ggplot")
})

test_that("AnnotationPlot extends region with large flanks", {
  p <- AnnotationPlot(
    object = atac_small,
    region = "chr1:780000-781000",
    extend.upstream = 50000,
    extend.downstream = 50000
  )
  expect_s3_class(p, "ggplot")
})

test_that("AnnotationPlot for transcript mode with empty region", {
  p <- AnnotationPlot(
    object = atac_small,
    region = "chr1:1-100",
    mode = "transcript"
  )
  expect_s3_class(p, "ggplot")
})

# theme_browser ----------------------------------------------------------------

test_that("theme_browser returns a theme", {
  th <- theme_browser()
  expect_s3_class(th, "theme")
  th2 <- theme_browser(legend = FALSE, axis.text.y = TRUE)
  expect_s3_class(th2, "theme")
})

# CombineTracks ----------------------------------------------------------------

test_that("CombineTracks combines plot list", {
  p1 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  p2 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_line()
  p3 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  combined <- CombineTracks(plotlist = list(p1, p2, p3))
  expect_s3_class(combined, c("patchwork", "ggplot"), exact = FALSE)
})

test_that("CombineTracks returns a single plot if only one given", {
  p1 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  combined <- CombineTracks(plotlist = list(p1))
  expect_s3_class(combined, "ggplot")
})

test_that("CombineTracks errors on wrong-length heights", {
  p1 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  p2 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  expect_error(
    CombineTracks(plotlist = list(p1, p2), heights = c(1, 2, 3)),
    regexp = "Relative height"
  )
})

test_that("CombineTracks removes NULL plots", {
  p1 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  p2 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_line()
  combined <- CombineTracks(plotlist = list(p1, NULL, p2))
  expect_s3_class(combined, c("patchwork", "ggplot"), exact = FALSE)
})

test_that("CombineTracks with expression plot", {
  p1 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  p2 <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_line()
  ep <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  combined <- CombineTracks(
    plotlist = list(p1, p2),
    expression.plot = ep,
    heights = c(8, 1),
    widths = c(3, 1)
  )
  expect_true(inherits(combined, "ggplot") || inherits(combined, "patchwork"))
})

# MotifPlot --------------------------------------------------------------------

test_that("MotifPlot validates motif argument type", {
  expect_error(
    MotifPlot(object = atac_small, motifs = 1:3),
    regexp = "motif names"
  )
})

# FragmentHistogram ------------------------------------------------------------

test_that("FragmentHistogram works with fragments attached", {
  obj <- setup_obj()
  p <- FragmentHistogram(
    object = obj, region = "chr1:1-2000000", group.by = NULL
  )
  expect_s3_class(p, "ggplot")
})

test_that("FragmentHistogram with group.by and log.scale", {
  obj <- setup_obj()
  p <- FragmentHistogram(
    object = obj, region = "chr1:1-2000000", group.by = "cluster"
  )
  expect_s3_class(p, "ggplot")
  p_log <- FragmentHistogram(
    object = obj, region = "chr1:1-2000000", log.scale = TRUE
  )
  expect_s3_class(p_log, "ggplot")
})

# CoveragePlot -----------------------------------------------------------------

test_that("CoveragePlot works with fragments", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE
  )
  expect_s3_class(p, c("patchwork", "ggplot"), exact = FALSE)
})

test_that("CoveragePlot with annotation and peaks", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = TRUE, peaks = TRUE
  )
  expect_s3_class(p, c("patchwork", "ggplot"), exact = FALSE)
})

test_that("CoveragePlot with extension and window", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE,
    extend.upstream = 1000, extend.downstream = 1000, window = 200
  )
  expect_s3_class(p, c("patchwork", "ggplot"), exact = FALSE)
})

test_that("CoveragePlot with show.bulk", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, show.bulk = TRUE
  )
  expect_s3_class(p, c("patchwork", "ggplot"), exact = FALSE)
})

test_that("CoveragePlot with multiple regions", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj,
    region = c("chr1:780000-790000", "chr1:850000-860000"),
    annotation = FALSE, peaks = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with assay.scale = separate", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, assay.scale = "separate"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot errors on invalid assay.scale", {
  obj <- setup_obj()
  expect_error(
    CoveragePlot(
      object = obj, region = "chr1:780000-790000",
      annotation = FALSE, peaks = FALSE, assay.scale = "wrong"
    ),
    regexp = "Unknown assay.scale"
  )
})

test_that("CoveragePlot with split.by", {
  obj <- setup_obj()
  obj$split_test <- sample(c("a", "b"), ncol(obj), replace = TRUE)
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, split.by = "split_test"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with ymax quantile", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, ymax = "q90"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with ranges", {
  obj <- setup_obj()
  rg <- GRanges("chr1", IRanges(c(782000, 786000), c(783000, 787000)))
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, ranges = rg
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with idents", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, idents = "1"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with region.highlight", {
  obj <- setup_obj()
  rh <- GRanges("chr1", IRanges(782000, 783000))
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, region.highlight = rh
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot tile=TRUE", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, tile = TRUE, tile.cells = 10
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with features and expression.assay", {
  obj <- setup_obj()
  features <- head(rownames(obj[["RNA"]]), 2)
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    features = features, expression.assay = "RNA",
    annotation = FALSE, peaks = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with bulk and group.by", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE,
    show.bulk = TRUE, group.by = "cluster"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with scale.factor", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, scale.factor = 1e6
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with peaks.group.by", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:710000-780000",
    annotation = FALSE, peaks = TRUE, peaks.group.by = "count"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with ranges.group.by", {
  obj <- setup_obj()
  rg <- GRanges("chr1", IRanges(c(782000, 786000), c(783000, 787000)))
  rg$category <- c("A", "B")
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE,
    ranges = rg, ranges.group.by = "category", ranges.title = "MyRanges"
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with links present", {
  obj <- setup_obj()
  gi <- InteractionSet::GInteractions(
    GRanges("chr1", IRanges(c(782000, 786000), c(782500, 786500))),
    GRanges("chr1", IRanges(c(784000, 788000), c(784500, 788500)))
  )
  gi$score <- c(0.6, 0.8)
  Links(obj) <- list(linkpeaks = gi)
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with annotation mode 'gene'", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = "gene", peaks = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with annotation mode 'transcript'", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = "transcript", peaks = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with gwas argument", {
  obj <- setup_obj()
  gwas_df <- data.frame(
    chromosome = "chr1",
    base_pair_location = c(782100, 785000),
    p_value = c(1e-8, 1e-4)
  )
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, gwas = gwas_df
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with variants argument", {
  obj <- setup_obj()
  variants <- data.frame(
    position = c(782100, 785000),
    rsid = c("rs1", "rs2"),
    color = c("steelblue", "darkred")
  )
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, variants = variants
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with multiple assays as list", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, assay = list("peaks")
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("CoveragePlot with split.assays = TRUE", {
  obj <- setup_obj()
  p <- CoveragePlot(
    object = obj, region = "chr1:780000-790000",
    annotation = FALSE, peaks = FALSE, split.assays = TRUE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

# MultiCoveragePlot ------------------------------------------------------------

test_that("MultiCoveragePlot returns a plot", {
  obj <- setup_obj()
  regions <- list("chr1:780000-790000", "chr1:850000-860000")
  p <- MultiCoveragePlot(
    object = obj, regions = regions,
    region_names = c("r1", "r2"),
    annotation = FALSE, peaks = FALSE
  )
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("MultiCoveragePlot errors on bad assay.scale", {
  obj <- setup_obj()
  expect_error(
    MultiCoveragePlot(
      object = obj, regions = list("chr1:780000-790000"),
      assay.scale = "bad", annotation = FALSE, peaks = FALSE
    ),
    regexp = "assay.scale"
  )
})

test_that("MultiCoveragePlot errors on bad bigwig.type", {
  obj <- setup_obj()
  expect_error(
    MultiCoveragePlot(
      object = obj, regions = list("chr1:780000-790000"),
      bigwig.type = "heatmap", annotation = FALSE, peaks = FALSE
    ),
    regexp = "bigwig.type"
  )
})

test_that("MultiCoveragePlot errors on region.highlight length mismatch", {
  obj <- setup_obj()
  expect_error(
    MultiCoveragePlot(
      object = obj, regions = list("chr1:780000-790000"),
      region.highlight = list(
        GRanges("chr1", IRanges(1, 100)),
        GRanges("chr1", IRanges(1, 100))
      ),
      annotation = FALSE, peaks = FALSE
    ),
    regexp = "region.highlight"
  )
})

# TilePlot ---------------------------------------------------------------------

test_that("TilePlot returns a ggplot", {
  obj <- setup_obj()
  p <- TilePlot(object = obj, region = "chr1:713500-714500", tile.cells = 10)
  expect_s3_class(p, "ggplot")
})

# ExpressionPlot ---------------------------------------------------------------

test_that("ExpressionPlot returns a ggplot", {
  Idents(atac_small) <- "cluster"
  features <- head(rownames(atac_small[["RNA"]]), 2)
  p <- ExpressionPlot(object = atac_small, features = features, assay = "RNA")
  expect_s3_class(p, "ggplot")
})

test_that("ExpressionPlot with multiple features and group.by", {
  obj <- setup_obj()
  features <- head(rownames(obj[["RNA"]]), 3)
  p <- ExpressionPlot(
    object = obj, features = features, assay = "RNA", group.by = "cluster"
  )
  expect_s3_class(p, "ggplot")
})

test_that("ExpressionPlot with idents and missing levels", {
  obj <- setup_obj()
  features <- rownames(obj[["RNA"]])[1]
  p <- ExpressionPlot(
    object = obj, features = features, assay = "RNA", idents = "1"
  )
  expect_s3_class(p, "ggplot")
})

test_that("ExpressionPlot errors when no features found", {
  expect_error(
    ExpressionPlot(
      object = atac_small, features = "ZZZZNotAGene", assay = "RNA"
    ),
    regexp = "None of the requested"
  )
})

test_that("ExpressionPlot warns when some features missing", {
  features <- c(rownames(atac_small[["RNA"]])[1], "ZZZZNotAGene")
  expect_warning(
    ExpressionPlot(object = atac_small, features = features, assay = "RNA"),
    regexp = "Some features"
  )
})

# LinkPlot ---------------------------------------------------------------------

test_that("LinkPlot returns NULL when links are empty", {
  links <- InteractionSet::GInteractions()
  Links(atac_small)[["linkpeaks"]] <- links
  p <- LinkPlot(
    object = atac_small,
    region = "chr1:780000-790000",
    key = "linkpeaks"
  )
  expect_null(p)
})

# VariantPlot/VariantTrack -----------------------------------------------------

test_that("VariantTrack carries variant rows into plot data", {
  variants <- data.frame(
    position = c(780100, 780500),
    rsid = c("rs1", "rs2"),
    color = c("steelblue", "darkred")
  )
  p <- VariantTrack(variants = variants, region = "chr1:780000-781000")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 2)
  expect_equal(p$data$position, c(780100, 780500))
  expect_equal(p$data$rsid, c("rs1", "rs2"))
  # region sets the x-axis limits (variants themselves are not filtered)
  expect_equal(p$scales$get_scales("x")$limits, c(780000, 781000))
  # a GRanges region produces the same x-axis limits as the string form
  gr <- GRanges("chr1", IRanges(780000, 781000))
  p2 <- VariantTrack(variants = variants, region = gr)
  expect_equal(p2$scales$get_scales("x")$limits, c(780000, 781000))
})

test_that("VariantPlot filters by n_cells_conf_detected >= min.cells", {
  variants <- data.frame(
    n_cells_conf_detected = c(5, 10, 1),
    vmr = c(0.05, 0.005, 0.5),
    strand_correlation = c(0.7, 0.3, 0.9)
  )
  p <- VariantPlot(variants = variants, min.cells = 2)
  expect_s3_class(p, "ggplot")
  # min.cells = 2 drops the row with n_cells_conf_detected = 1
  expect_equal(nrow(p$data), 2)
  expect_setequal(p$data$n_cells_conf_detected, c(5, 10))
  # 'pos' column flags variants above both vmr and concordance thresholds
  expect_true("pos" %in% colnames(p$data))
})

# BigwigTrack ------------------------------------------------------------------

test_that("BigwigTrack errors on unknown type", {
  expect_error(
    BigwigTrack(region = "chr1:1000-2000", bigwig = "fake.bw", type = "invalid"),
    regexp = "Invalid type"
  )
})

# GWASTrack --------------------------------------------------------------------

test_that("GWASTrack basic with valid data", {
  gwas <- data.frame(
    chromosome = "chr1",
    base_pair_location = c(782100, 785000, 788000),
    p_value = c(1e-8, 1e-4, 1e-6)
  )
  p <- GWASTrack(gwas = gwas, region = "chr1:780000-790000")
  expect_s3_class(p, "ggplot")
})

test_that("GWASTrack with show.axis FALSE", {
  gwas <- data.frame(
    chromosome = "chr1",
    base_pair_location = c(782100, 785000),
    p_value = c(1e-8, 1e-4)
  )
  p <- GWASTrack(
    gwas = gwas, region = "chr1:780000-790000", show.axis = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("GWASTrack with r2 column", {
  gwas <- data.frame(
    chromosome = "chr1",
    base_pair_location = c(782100, 785000),
    p_value = c(1e-8, 1e-4),
    r2 = c(0.9, 0.3)
  )
  p <- GWASTrack(gwas = gwas, region = "chr1:780000-790000")
  expect_s3_class(p, "ggplot")
})

test_that("GWASTrack with credset column", {
  gwas <- data.frame(
    chromosome = "chr1",
    base_pair_location = c(782100, 785000),
    p_value = c(1e-8, 1e-4),
    in_credset = c(TRUE, FALSE)
  )
  p <- GWASTrack(gwas = gwas, region = "chr1:780000-790000")
  expect_s3_class(p, "ggplot")
})

test_that("GWASTrack with LD file", {
  gwas_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tbase_pair_location\tp_value",
    "chr1\t782100\t1e-8",
    "chr1\t785000\t1e-4"
  ), gwas_tf)
  ld_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tr2",
    "chr1\t782100\t0.9",
    "chr1\t785000\t0.3"
  ), ld_tf)
  p <- GWASTrack(
    gwas = gwas_tf, region = "chr1:780000-790000",
    ld.file = ld_tf, ld.lead.snp = "rs1"
  )
  expect_s3_class(p, "ggplot")
})

test_that("GWASTrack errors when ld.file without lead.snp", {
  gwas <- data.frame(
    chromosome = "chr1", base_pair_location = 782100, p_value = 1e-8
  )
  expect_error(
    GWASTrack(
      gwas = gwas, region = "chr1:780000-790000", ld.file = "fake.txt"
    ),
    regexp = "ld.lead.snp"
  )
})

test_that("GWASTrack errors when no data in region", {
  gwas <- data.frame(
    chromosome = "chr1", base_pair_location = 99999999, p_value = 1e-8
  )
  expect_error(
    GWASTrack(gwas = gwas, region = "chr1:780000-790000"),
    regexp = "No GWAS data"
  )
})

test_that("GWASTrack with credset file", {
  gwas_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tbase_pair_location\tp_value",
    "chr1\t782100\t1e-8",
    "chr1\t785000\t1e-4"
  ), gwas_tf)
  cs_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tpip\tcs",
    "chr1\t782100\t0.5\t1"
  ), cs_tf)
  p <- GWASTrack(
    gwas = gwas_tf, region = "chr1:780000-790000", credset.file = cs_tf
  )
  expect_s3_class(p, "ggplot")
})

test_that("GWASTrack with both LD and credset", {
  gwas_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tbase_pair_location\tp_value",
    "chr1\t782100\t1e-8",
    "chr1\t785000\t1e-4"
  ), gwas_tf)
  ld_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tr2",
    "chr1\t782100\t0.9",
    "chr1\t785000\t0.3"
  ), ld_tf)
  cs_tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "chromosome\tposition\tpip\tcs",
    "chr1\t782100\t0.5\t1"
  ), cs_tf)
  p <- GWASTrack(
    gwas = gwas_tf, region = "chr1:780000-790000",
    ld.file = ld_tf, ld.lead.snp = "rs1", credset.file = cs_tf
  )
  expect_s3_class(p, "ggplot")
})

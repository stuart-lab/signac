library(GenomicRanges)

# Simulated AVI SHAP feature attributions: three alt alleles per position
make_avi_df <- function(chrom = "chr1", positions = c(782100, 785000, 788000)) {
  set.seed(1)
  n <- length(positions) * 3
  data.frame(
    chromosome = chrom,
    position = rep(positions, each = 3),
    ref = "A",
    alt = rep(c("C", "G", "T"), times = length(positions)),
    MAX_ABS_ATAC = round(runif(n, 0, 1), 3),
    MAX_ABS_DNASE = round(runif(n, 0, 1), 3),
    MERGED_SPLICING = round(runif(n, 0, 2), 3),
    ALPHAMISSENSE = 0,
    CACTUS_241_WAY = round(runif(n, -0.5, 0.5), 3),
    avi_raw = 0,
    stringsAsFactors = FALSE
  )
}

# Column handling --------------------------------------------------------------

test_that("CheckAVIColumns requires the variant columns", {
  df <- data.frame(chromosome = "chr1", position = 10, ref = "A", alt = "T", x = 1)
  checked <- CheckAVIColumns(df)
  expect_type(checked$position, "integer")
  expect_error(
    CheckAVIColumns(data.frame(chrom = "chr1", ref = "A", alt = "T")),
    regexp = "Missing: chromosome, position"
  )
  expect_error(
    AVITrack(data.frame(chrom = "chr1", pos = 1, ref = "A", alt = "T", MAX_ABS_ATAC = 1),
             region = "chr1:1-10"),
    regexp = "Missing"
  )
})

test_that("AVIFeatureColumns keeps only columns mapping to a known group", {
  avi <- make_avi_df()
  avi$avi_phred <- 10
  avi$label <- "x"
  avi$rank <- 7
  avi$n_tracks <- 3L
  avi$AVI_SCORE <- 1.5
  feats <- AVIFeatureColumns(avi)
  expect_setequal(
    feats,
    c("MAX_ABS_ATAC", "MAX_ABS_DNASE", "MERGED_SPLICING", "ALPHAMISSENSE",
      "CACTUS_241_WAY")
  )
  # tables without feature columns
  scores <- avi[, c("chromosome", "position", "ref", "alt", "AVI_SCORE")]
  expect_error(AVIFeatureColumns(scores), regexp = "No AVI SHAP feature")
  expect_equal(AVIFeatureColumns(scores, allow.empty = TRUE), character())
})

test_that("AVIFeatureGroups maps the 18 Atlas features to modality groups", {
  feats <- c(
    "MAX_ABS_ATAC", "MAX_ABS_DNASE", "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE",
    "MAX_ABS_CAGE", "MAX_ABS_PROCAP", "MAX_ABS_RNA_SEQ",
    "MAX_ABS_POLYADENYLATION", "MERGED_SPLICING", "MAX_ABS_CONTACT_MAPS",
    "ALPHAMISSENSE", "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST",
    "PHASTCONS_470_WAY", "CACTUS_241_WAY", "IS_INSERTION", "IS_DELETION"
  )
  groups <- AVIFeatureGroups(feats)
  expect_equal(names(groups), feats)
  expect_equal(unname(groups[c("MAX_ABS_ATAC", "MAX_ABS_DNASE")]),
               rep("Accessibility", 2))
  expect_equal(unname(groups["MAX_ABS_CHIP_TF"]), "TF binding")
  expect_equal(unname(groups["MAX_ABS_CHIP_HISTONE"]), "Histone")
  expect_equal(unname(groups[c("MAX_ABS_CAGE", "MAX_ABS_PROCAP", "MAX_ABS_RNA_SEQ")]),
               rep("Expression", 3))
  expect_equal(unname(groups["MAX_ABS_POLYADENYLATION"]), "Polyadenylation")
  expect_equal(unname(groups["MERGED_SPLICING"]), "Splicing")
  expect_equal(unname(groups["MAX_ABS_CONTACT_MAPS"]), "Contact maps")
  expect_equal(unname(groups[c("ALPHAMISSENSE", "PROTEIN_TERMINATION",
                               "START_LOST", "STOP_LOST")]),
               rep("Protein", 4))
  expect_equal(unname(groups[c("PHASTCONS_470_WAY", "CACTUS_241_WAY")]),
               rep("Conservation", 2))
  expect_equal(unname(groups[c("IS_INSERTION", "IS_DELETION")]), rep("Indel", 2))
  # lower-case and unknown names; matching is exact apart from case
  expect_equal(unname(AVIFeatureGroups("max_abs_atac")), "Accessibility")
  expect_equal(unname(AVIFeatureGroups("custom")), "custom")
  expect_equal(unname(AVIFeatureGroups("MAX_ABS_ATAC_EXTRA")), "MAX_ABS_ATAC_EXTRA")
  expect_equal(unname(AVIFeatureGroups("protein_score")), "protein_score")
  # the feature table lists each feature once with a color per group
  tbl <- AVIFeatureTable()
  expect_equal(nrow(tbl), 18)
  expect_equal(anyDuplicated(tbl$feature), 0)
  expect_equal(nrow(unique(tbl[, c("group", "color")])), 10)
})

# Collapsing, binning, feature resolution --------------------------------------

test_that("CollapseAVIAlleles keeps the variant with the largest |total| per position", {
  avi <- make_avi_df()
  feats <- c("MAX_ABS_ATAC", "MERGED_SPLICING", "CACTUS_241_WAY")
  # without a Total column the requested features are summed for ranking
  collapsed <- CollapseAVIAlleles(avi, features = feats)
  expect_equal(nrow(collapsed), 3)
  expect_equal(collapsed$position, c(782100, 785000, 788000))
  totals <- rowSums(avi[, feats])
  for (pos in unique(avi$position)) {
    idx <- which(avi$position == pos)
    best <- idx[which.max(abs(totals[idx]))]
    row <- collapsed[collapsed$position == pos, ]
    expect_equal(row$MERGED_SPLICING, avi$MERGED_SPLICING[best])
    expect_equal(row$alt, avi$alt[best])
    expect_equal(row$Total, totals[best])
  }
  # with a Total column (all features) the ranking is independent of the
  # displayed features
  with.total <- avi
  with.total$Total <- rowSums(avi[, AVIFeatureColumns(avi)])
  collapsed <- CollapseAVIAlleles(with.total, features = "MAX_ABS_ATAC")
  for (pos in unique(avi$position)) {
    idx <- which(avi$position == pos)
    best <- idx[which.max(abs(with.total$Total[idx]))]
    row <- collapsed[collapsed$position == pos, ]
    expect_equal(row$alt, avi$alt[best])
    expect_equal(row$Total, with.total$Total[best])
  }
})

test_that("BinAVI fills positions and keeps the largest-magnitude value per bin", {
  values <- data.frame(position = c(1005L, 1010L), a = c(1, -3), b = c(2, 0))
  region <- GRanges("chr1:1001-1020")
  # one bin per position when the region is small
  long <- BinAVI(values, region, tracks = c("a", "b"), bins = 3000)
  expect_equal(attr(long, "bin.size"), 1)
  expect_equal(nrow(long), 20 * 2)
  expect_equal(long$value[long$group == "a" & long$position == 1005], 1)
  expect_equal(long$value[long$group == "a" & long$position == 1010], -3)
  expect_equal(long$value[long$group == "b" & long$position == 1010], 0)
  # binning keeps the signed value of largest magnitude in each bin
  binned <- BinAVI(values, region, tracks = c("a", "b"), bins = 2)
  expect_equal(attr(binned, "bin.size"), 10)
  expect_equal(nrow(binned), 2 * 2)
  a <- binned[binned$group == "a", ]
  expect_equal(a$position, c(1005.5, 1015.5))
  expect_equal(a$value, c(-3, 0))
  b <- binned[binned$group == "b", ]
  expect_equal(b$value, c(2, 0))
  expect_equal(a$position, b$position)
  # rank.by selects one common position per bin for all tracks
  ranked <- BinAVI(values, region, tracks = c("a", "b"), bins = 2, rank.by = "b")
  expect_equal(ranked$value[ranked$group == "a"], c(1, 0))
  expect_equal(ranked$value[ranked$group == "b"], c(2, 0))
  # large regions are reduced to exactly `bins` bins of equal width
  wide <- BinAVI(values, GRanges("chr1:1-100000"), tracks = "a", bins = 200)
  expect_equal(attr(wide, "bin.size"), 500)
  expect_equal(nrow(wide), 200)
  expect_equal(wide$position[1:2], c(250.5, 750.5))
  odd <- BinAVI(values, GRanges("chr1:1-1001"), tracks = "a", bins = 200)
  expect_equal(nrow(odd), 200)
})

test_that("ResolveAVIFeatures accepts column names and group labels", {
  avi <- make_avi_df()
  expect_setequal(
    ResolveAVIFeatures(avi, "Accessibility"), c("MAX_ABS_ATAC", "MAX_ABS_DNASE")
  )
  expect_setequal(
    ResolveAVIFeatures(avi, c("accessibility", "Splicing")),
    c("MAX_ABS_ATAC", "MAX_ABS_DNASE", "MERGED_SPLICING")
  )
  expect_equal(ResolveAVIFeatures(avi, "MAX_ABS_ATAC"), "MAX_ABS_ATAC")
  expect_equal(ResolveAVIFeatures(avi, "max_abs_atac"), "MAX_ABS_ATAC")
  # aggregate score columns can be requested explicitly
  expect_equal(ResolveAVIFeatures(avi, "avi_raw"), "avi_raw")
  # "AVI" is a plain column name, not an alias for the total
  avi$AVI <- 1
  expect_equal(ResolveAVIFeatures(avi, "AVI"), "AVI")
  expect_equal(length(ResolveAVIFeatures(avi, NULL)), 5)
  # known groups without features in the data are skipped with a message
  expect_message(
    feats <- ResolveAVIFeatures(avi, c("Histone", "Accessibility")),
    regexp = "Histone"
  )
  expect_setequal(feats, c("MAX_ABS_ATAC", "MAX_ABS_DNASE"))
  expect_error(ResolveAVIFeatures(avi, "Histone"), regexp = "None of the")
  # unknown names still error
  expect_error(ResolveAVIFeatures(avi, "nope"), regexp = "not found")
})

# AVITrack ---------------------------------------------------------------------

test_that("AVITrack groups all features and sums within groups", {
  avi <- make_avi_df()
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000", features = NULL, bins = Inf
  )
  expect_setequal(
    levels(p$data$group),
    c("Accessibility", "Splicing", "Protein", "Conservation")
  )
  expect_equal(nrow(p$data), 10001 * 4)
  # accessibility = ATAC + DNase of the max-|total| variant at each position
  with.total <- avi
  with.total$Total <- rowSums(avi[, AVIFeatureColumns(avi)])
  collapsed <- CollapseAVIAlleles(with.total, AVIFeatureColumns(avi))
  acc <- p$data[p$data$group == "Accessibility", ]
  expect_equal(
    acc$value[match(collapsed$position, acc$position)],
    collapsed$MAX_ABS_ATAC + collapsed$MAX_ABS_DNASE
  )
  # positions without variants are zero
  expect_equal(sum(acc$value != 0), 3)
})

test_that("AVITrack plots the total AVI score as the sum of attributions", {
  avi <- make_avi_df()
  feats <- AVIFeatureColumns(avi)
  # total alone
  p <- AVITrack(avi = avi, region = "chr1:780000-790000", features = "total",
                bins = Inf)
  expect_equal(levels(p$data$group), "Total")
  expect_equal(nrow(p$data), 10001)
  totals <- rowSums(avi[, feats])
  best <- tapply(seq_len(nrow(avi)), avi$position, function(i) i[which.max(abs(totals[i]))])
  expect_equal(
    p$data$value[match(as.integer(names(best)), p$data$position)],
    unname(totals[best])
  )
  # total with components: components are stacked, total is a separate layer
  p <- AVITrack(avi = avi, region = "chr1:780000-790000",
                features = c("Accessibility", "Splicing", "total impact"))
  expect_setequal(levels(p$data$group), c("Accessibility", "Splicing", "Total"))
  expect_false("Total" %in% p$data$group)
  expect_equal(length(p$layers), 3)
  expect_s3_class(p$layers[[3]]$geom, "GeomPoint")
  # the total appears in the legend as a diamond
  expect_true("shape" %in% names(p$layers[[3]]$mapping))
  expect_equal(unname(p$scales$get_scales("shape")$palette(1)), 18)
})

test_that("AVITrack keeps the requested feature order in the factor levels", {
  avi <- make_avi_df()
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000",
    features = c("Splicing", "Conservation", "Accessibility")
  )
  expect_equal(levels(p$data$group), c("Splicing", "Conservation", "Accessibility"))
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000",
    features = c("CACTUS_241_WAY", "MAX_ABS_DNASE"), group.features = FALSE
  )
  expect_equal(levels(p$data$group), c("CACTUS_241_WAY", "MAX_ABS_DNASE"))
  # features = NULL follows the column order of the data
  p <- AVITrack(avi = avi, region = "chr1:780000-790000", features = NULL)
  expect_equal(
    levels(p$data$group),
    c("Accessibility", "Splicing", "Protein", "Conservation")
  )
  # mixed groups and total, in the given order (reversed for the heatmap rows)
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000",
    features = c("Total", "Splicing", "Accessibility"), type = "heatmap"
  )
  expect_equal(levels(p$data$group), rev(c("Total", "Splicing", "Accessibility")))
})

test_that("AVITrack supports the heatmap type", {
  avi <- make_avi_df()
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000", type = "heatmap",
    features = NULL, bins = 100
  )
  expect_s3_class(p, "ggplot")
  # one tile per bin per group
  expect_equal(nrow(p$data), 100 * 4)
  expect_error(
    AVITrack(avi = avi, region = "chr1:780000-790000", type = "line"),
    regexp = "Invalid type"
  )
})

test_that("AVITrack respects features, group.features, and colors", {
  avi <- make_avi_df()
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000",
    features = c("MAX_ABS_ATAC", "CACTUS_241_WAY"), group.features = FALSE,
    colors = c(MAX_ABS_ATAC = "red", CACTUS_241_WAY = "blue"),
    ymax = 2, y_label = "SHAP", show.axis = FALSE
  )
  expect_s3_class(p, "ggplot")
  expect_setequal(levels(p$data$group), c("MAX_ABS_ATAC", "CACTUS_241_WAY"))
  expect_equal(p$labels$y, "SHAP")
  expect_equal(
    unname(p$scales$get_scales("fill")$palette(2)),
    c("red", "blue")
  )
  expect_error(
    AVITrack(avi = avi, region = "chr1:780000-790000", features = "Histone"),
    regexp = "None of the"
  )
})

test_that("AVITrack draws an empty track when no variants are in the region", {
  expect_message(
    p <- AVITrack(avi = make_avi_df(), region = "chr1:1000-2000",
                  features = "Accessibility"),
    regexp = "No AVI data"
  )
  expect_s3_class(p, "ggplot")
  expect_true(all(p$data$value == 0))
  expect_error(AVITrack(avi = 1, region = "chr1:1000-2000"), regexp = "avi must")
  expect_error(AVITrack(avi = "file.tsv.gz", region = "chr1:1000-2000"),
               regexp = "avi must")
})

test_that("AVITrack fills missing user colors from the default palette", {
  avi <- make_avi_df()
  p <- AVITrack(
    avi = avi, region = "chr1:780000-790000",
    features = c("Accessibility", "Total"), colors = c(Accessibility = "gold")
  )
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 3)
  expect_equal(p$layers[[3]]$aes_params$colour, "#000000")
  expect_equal(unname(p$scales$get_scales("fill")$palette(2))[1], "gold")
})

# AlphaGenome Atlas API -------------------------------------------------------

# A fake Python Atlas client with the same interface as
# alphagenome.atlas.atlas.AtlasClient, returning simulated scores
make_fake_atlas_client <- function(interval.available = TRUE, fail.first = 0L) {
  skip_if_not_installed("reticulate")
  skip_if_not(
    reticulate::py_module_available("alphagenome"),
    "alphagenome Python package not available"
  )
  reticulate::py_run_string('
import numpy as np, pandas as pd, anndata
from alphagenome.data import genome
from alphagenome.atlas.atlas import ScorerMetadata

FEATURES = ["MAX_ABS_ATAC", "MERGED_SPLICING", "CACTUS_241_WAY"]

class FakeAtlasClient:
    def __init__(self, interval_available=True, fail_first=0):
        self.calls = []
        self.metadata_calls = []
        self.interval_available = interval_available
        # number of data requests that fail with UNAVAILABLE before succeeding
        self.fail_first = fail_first

    def scorer_metadata(self):
        self.metadata_calls.append(1)
        return {
            "AVI": ScorerMetadata(name="AVI", is_signed=False,
                track_metadata=pd.DataFrame({"name": ["score"]})),
            "AVI_SCORE_FEATURE_IMPORTANCE": ScorerMetadata(
                name="AVI_SCORE_FEATURE_IMPORTANCE", is_signed=False,
                track_metadata=pd.DataFrame({"name": FEATURES})),
        }

    def query_interval(self, interval, requested_scorers, progress_bar=True, **kw):
        self.calls.append((interval.chromosome, interval.start, interval.end,
                           list(requested_scorers)))
        if not self.interval_available or self.fail_first > 0:
            self.fail_first -= 1
            raise RuntimeError(
                "grpc._channel._InactiveRpcError: status = StatusCode.UNAVAILABLE "
                "details = The service is currently unavailable.")
        variants = [genome.Variant(interval.chromosome, pos, "A", alt)
                    for pos in range(interval.start + 1, interval.end + 1)
                    for alt in "CGT"]
        idx = [str(i) for i in range(len(variants))]
        out = {}
        for scorer in requested_scorers:
            if scorer == "AVI_SCORE_FEATURE_IMPORTANCE":
                X = np.arange(len(variants) * len(FEATURES), dtype=np.float32)
                X = X.reshape(len(variants), len(FEATURES)) / 10.0
                var = pd.DataFrame({"name": FEATURES, "strand": "."},
                                   index=[str(i) for i in range(len(FEATURES))])
            elif scorer == "AVI":
                X = np.arange(len(variants), dtype=np.float32).reshape(-1, 1)
                var = pd.DataFrame({"name": ["score"], "strand": "."}, index=["0"])
            else:
                continue
            obs = pd.DataFrame({"variant": variants}, index=idx)
            out[scorer] = anndata.AnnData(X=X, obs=obs, var=var)
        return out
')
  reticulate::py$FakeAtlasClient(
    interval_available = interval.available, fail_first = fail.first
  )
}

test_that("AtlasScorerTable describes the fake client", {
  client <- make_fake_atlas_client()
  scorers <- AtlasScorerTable(client)
  expect_equal(scorers$name, c("AVI", "AVI_SCORE_FEATURE_IMPORTANCE"))
  expect_equal(scorers$is_signed, c(FALSE, FALSE))
  expect_equal(scorers$n_scores, c(1L, 3L))
})

test_that("ListAtlasScorers returns the scorer table for the client", {
  client <- make_fake_atlas_client()
  local_mocked_bindings(AtlasClient = function(api.key = NULL) client)
  scorers <- ListAtlasScorers(api.key = "x")
  expect_equal(scorers, AtlasScorerTable(client))
})

test_that("QueryAtlas converts Atlas results to a variant data.frame", {
  client <- make_fake_atlas_client()
  avi <- QueryAtlas(client, region = GRanges("chr1:1001-1003"), verbose = FALSE)
  expect_s3_class(avi, "data.frame")
  # 3 positions x 3 alt alleles
  expect_equal(nrow(avi), 9)
  expect_equal(
    colnames(avi),
    c("chromosome", "position", "ref", "alt",
      "MAX_ABS_ATAC", "MERGED_SPLICING", "CACTUS_241_WAY")
  )
  expect_equal(unique(avi$chromosome), "chr1")
  expect_equal(avi$position, rep(1001:1003, each = 3))
  expect_equal(avi$alt, rep(c("C", "G", "T"), times = 3))
  expect_equal(avi$ref, rep("A", 9))
  # scores are float32 in the Atlas
  expect_equal(avi$MAX_ABS_ATAC, seq(0, 2.4, by = 0.3), tolerance = 1e-6)
  expect_equal(avi$CACTUS_241_WAY, seq(0.2, 2.6, by = 0.3), tolerance = 1e-6)
  # the interval passed to Python is 0-based half-open
  calls <- reticulate::py_to_r(client$calls)
  expect_equal(calls[[1]][[1]], "chr1")
  expect_equal(calls[[1]][[2]], 1000L)
  expect_equal(calls[[1]][[3]], 1003L)
  expect_equal(unlist(calls[[1]][[4]]), "AVI_SCORE_FEATURE_IMPORTANCE")
  # no scorer metadata is fetched for a query
  expect_equal(length(reticulate::py_to_r(client$metadata_calls)), 0)
})

test_that("QueryAtlas handles single-score scorers and unknown scorers", {
  client <- make_fake_atlas_client()
  avi <- QueryAtlas(
    client, region = GRanges("chr1:1001-1002"), scorer = "AVI", verbose = FALSE
  )
  expect_equal(colnames(avi), c("chromosome", "position", "ref", "alt", "AVI"))
  expect_equal(avi$AVI, as.numeric(0:5))
  expect_error(
    QueryAtlas(client, region = GRanges("chr1:1001-1002"), scorer = "nope",
               verbose = FALSE),
    regexp = "no scores.*ListAtlasScorers"
  )
})

test_that("QueryAtlas splits the interval into 32 bp chunks", {
  client <- make_fake_atlas_client()
  avi <- QueryAtlas(client, region = GRanges("chr1:1001-1070"), verbose = FALSE)
  expect_equal(nrow(avi), 70 * 3)
  expect_equal(avi$position, rep(1001:1070, each = 3))
  calls <- reticulate::py_to_r(client$calls)
  expect_equal(length(calls), 3)
  expect_equal(sapply(calls, `[[`, 2), c(1000L, 1032L, 1064L))
  expect_equal(sapply(calls, `[[`, 3), c(1032L, 1064L, 1070L))
})

test_that("QueryAtlasCached only fetches uncovered parts of a region", {
  ClearAtlasCache()
  client <- make_fake_atlas_client()
  avi1 <- QueryAtlasCached(client, region = GRanges("chr1:1001-1010"), verbose = FALSE)
  expect_equal(nrow(avi1), 30)
  expect_equal(length(reticulate::py_to_r(client$calls)), 1)
  # overlapping region: only 1011-1020 is requested
  avi2 <- QueryAtlasCached(client, region = GRanges("chr1:1005-1020"), verbose = FALSE)
  expect_equal(nrow(avi2), 16 * 3)
  expect_equal(avi2$position, rep(1005:1020, each = 3))
  calls <- reticulate::py_to_r(client$calls)
  expect_equal(length(calls), 2)
  expect_equal(calls[[2]][[2]], 1010L)
  expect_equal(calls[[2]][[3]], 1020L)
  # fully covered region: no request
  avi3 <- QueryAtlasCached(client, region = GRanges("chr1:1002-1003"), verbose = FALSE)
  expect_equal(nrow(avi3), 6)
  expect_equal(length(reticulate::py_to_r(client$calls)), 2)
  # a region with a gap in the middle fetches both flanks
  avi4 <- QueryAtlasCached(client, region = GRanges("chr1:995-1025"), verbose = FALSE)
  expect_equal(nrow(avi4), 31 * 3)
  expect_equal(length(reticulate::py_to_r(client$calls)), 4)
  # different scorers are cached separately
  QueryAtlasCached(client, region = GRanges("chr1:1001-1002"), scorer = "AVI",
                   verbose = FALSE)
  expect_equal(length(reticulate::py_to_r(client$calls)), 5)
  # cache = FALSE always queries
  QueryAtlasCached(client, region = GRanges("chr1:1001-1002"), cache = FALSE,
                   verbose = FALSE)
  expect_equal(length(reticulate::py_to_r(client$calls)), 6)
  # clearing the cache forces a new request
  ClearAtlasCache()
  QueryAtlasCached(client, region = GRanges("chr1:1001-1002"), verbose = FALSE)
  expect_equal(length(reticulate::py_to_r(client$calls)), 7)
  ClearAtlasCache()
})

test_that("Python Atlas helper unit tests pass", {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  skip_if_not(
    reticulate::py_module_available("alphagenome"),
    "alphagenome Python package not available"
  )
  AtlasHelper()
  tests <- reticulate::import_from_path(
    "test_signac_atlas", path = normalizePath(file.path("..", "python"))
  )
  expect_true(tests$run_suite())
})

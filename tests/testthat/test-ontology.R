test_that("RankFeatures caps infinite scores at the correct sign", {
  # avg_log2FC sign * -log10(p_val); p_val == 0 gives an infinite score
  markers <- data.frame(
    avg_log2FC = c(2, -2, 1, -1, 0.5),
    p_val = c(0, 0, 1e-10, 1e-5, 0.5),
    row.names = c("UP_INF", "DOWN_INF", "UP", "DOWN", "NS")
  )
  rl <- RankFeatures(markers = markers)

  # no infinite (or NaN) values remain
  expect_false(any(is.infinite(x = rl)))
  expect_false(any(is.nan(x = rl)))

  # finite scores are unchanged
  expect_equal(rl[["UP"]], 10)
  expect_equal(rl[["DOWN"]], -5)

  # +Inf capped at the maximum finite score, -Inf at the minimum finite score
  expect_equal(rl[["UP_INF"]], 10)
  expect_equal(rl[["DOWN_INF"]], -5)

  # the strongly down-regulated (p = 0) feature stays at the BOTTOM, not the
  # top: the previous implementation moved it to the maximum
  expect_equal(unname(rl[["DOWN_INF"]]), min(rl))
  expect_equal(unname(rl[["UP_INF"]]), max(rl))
  expect_lt(rl[["DOWN_INF"]], rl[["NS"]])
  expect_lt(rl[["DOWN_INF"]], rl[["UP"]])

  # names are preserved
  expect_setequal(names(x = rl), rownames(x = markers))
})

test_that("RankFeatures leaves finite scores untouched", {
  markers <- data.frame(
    avg_log2FC = c(1, -1),
    p_val = c(0.01, 0.01),
    row.names = c("a", "b")
  )
  rl <- RankFeatures(markers = markers)
  expect_equal(unname(rl), c(2, -2))
})

test_that("RankFeatures handles an all-infinite ranking without NaN", {
  markers <- data.frame(
    avg_log2FC = c(1, -1),
    p_val = c(0, 0),
    row.names = c("u", "d")
  )
  rl <- RankFeatures(markers = markers)
  expect_false(any(is.infinite(x = rl)))
  expect_false(any(is.nan(x = rl)))
  expect_gt(rl[["u"]], rl[["d"]])
})

# Minimal Seurat object with clear differential structure for EnrichedTerms.
make_ontology_object <- function(variable = paste0("g", 1:50)) {
  set.seed(1)
  ng <- 200
  nc <- 80
  counts <- matrix(
    data = rpois(n = ng * nc, lambda = 3), nrow = ng,
    dimnames = list(paste0("g", 1:ng), paste0("c", 1:nc))
  )
  grp <- rep(x = c("A", "B"), each = nc / 2)
  # g1-20 strongly up in group A
  counts[1:20, grp == "A"] <- counts[1:20, grp == "A"] + 40L
  obj <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  obj$grp <- grp
  obj <- Seurat::NormalizeData(object = obj, verbose = FALSE)
  SeuratObject::VariableFeatures(object = obj) <- variable
  obj
}

test_that("EnrichedTerms errors clearly when no variable features are set", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("fgsea")
  obj <- make_ontology_object(variable = character(0))
  terms <- list(T_up = paste0("g", 1:15))
  expect_error(
    EnrichedTerms(
      object = obj, terms = terms, group.by = "grp",
      var.features = TRUE, verbose = FALSE
    ),
    "no variable features"
  )
})

test_that("EnrichedTerms detects an enriched term within the variable set", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("fgsea")
  obj <- make_ontology_object()
  # T_up: up-regulated genes (all variable); T_bg: background genes
  terms <- list(T_up = paste0("g", 1:15), T_bg = paste0("g", 100:140))
  res <- suppressWarnings(EnrichedTerms(
    object = obj, terms = terms, group.by = "grp", direction = "up",
    var.features = TRUE, verbose = FALSE,
    logfc.threshold = 0, min.pct = 0
  ))
  expect_type(res, "list")
  expect_true("A" %in% names(res))
  # the up-regulated term is recovered as enriched in group A
  expect_true("T_up" %in% res[["A"]]$pathway)
  expect_true(all(res[["A"]]$NES > 0))
})

test_that("EnrichedTerms merges a user-supplied features argument", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("fgsea")
  obj <- make_ontology_object()
  terms <- list(T_up = paste0("g", 1:15), T_bg = paste0("g", 100:140))
  # passing `features` via ... alongside var.features = TRUE must intersect
  # rather than raise an argument-conflict error in FindMarkers
  res <- suppressWarnings(EnrichedTerms(
    object = obj, terms = terms, group.by = "grp", direction = "up",
    var.features = TRUE, verbose = FALSE,
    features = paste0("g", 1:30),
    logfc.threshold = 0, min.pct = 0
  ))
  expect_type(res, "list")
  expect_true("T_up" %in% res[["A"]]$pathway)
})

test_that("EnrichedTerms runs with var.features = FALSE (full universe)", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("fgsea")
  obj <- make_ontology_object(variable = character(0))
  terms <- list(T_up = paste0("g", 1:15), T_bg = paste0("g", 100:140))
  res <- suppressWarnings(EnrichedTerms(
    object = obj, terms = terms, group.by = "grp", direction = "up",
    var.features = FALSE, verbose = FALSE,
    logfc.threshold = 0, min.pct = 0
  ))
  expect_type(res, "list")
  expect_true("A" %in% names(res))
})

test_that("EnrichedTerms padj.cutoff controls which terms are retained", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("fgsea")
  obj <- make_ontology_object()
  terms <- list(T_up = paste0("g", 1:15), T_bg = paste0("g", 100:140))
  n.terms <- function(x) sum(sapply(X = x, FUN = nrow))
  loose <- EnrichedTerms(
    object = obj, terms = terms, group.by = "grp",
    padj.cutoff = 1, verbose = FALSE
  )
  strict <- EnrichedTerms(
    object = obj, terms = terms, group.by = "grp",
    padj.cutoff = 1e-12, verbose = FALSE
  )
  expect_gt(n.terms(loose), n.terms(strict))
  expect_equal(n.terms(strict), 0)
})

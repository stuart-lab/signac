library(GenomicRanges)
library(SeuratObject)
library(Matrix)

# atac_small already carries GC.percent and sequence.length feature metadata,
# so LinkPeaks can be invoked directly without injecting fake RegionStats info.

# ConnectionsToLinks -----------------------------------------------------------

test_that("ConnectionsToLinks works", {
  conns <- data.frame(
    Peak1 = c("chr1:1-100", "chr1:200-300", "chr1:400-500"),
    Peak2 = c("chr1:150-200", "chr1:350-400", "chr1:550-600"),
    coaccess = c(0.1, 0.5, NA)
  )
  gi <- ConnectionsToLinks(conns = conns, threshold = 0)
  expect_s4_class(gi, "GInteractions")
  expect_equal(length(gi), 2)
  expect_equal(gi$score, c(0.1, 0.5))
  expect_true(all(is.na(gi$group)))
})

test_that("ConnectionsToLinks with threshold filters", {
  conns <- data.frame(
    Peak1 = c("chr1:1-100", "chr1:200-300"),
    Peak2 = c("chr1:150-200", "chr1:350-400"),
    coaccess = c(0.1, 0.5)
  )
  gi <- ConnectionsToLinks(conns = conns, threshold = 0.2)
  expect_equal(length(gi), 1)
  expect_equal(gi$score, 0.5)
})

test_that("ConnectionsToLinks with ccans assigns groups", {
  conns <- data.frame(
    Peak1 = c("chr1:1-100", "chr1:200-300"),
    Peak2 = c("chr1:150-200", "chr1:350-400"),
    coaccess = c(0.1, 0.5)
  )
  ccans <- data.frame(
    Peak = c("chr1:1-100", "chr1:200-300", "chr1:150-200", "chr1:350-400"),
    CCAN = c("A", "B", "A", "B")
  )
  gi <- ConnectionsToLinks(conns = conns, ccans = ccans, threshold = 0)
  expect_s4_class(gi, "GInteractions")
  expect_equal(length(gi), 2)
  expect_true(all(!is.na(gi$group)))
})

test_that("ConnectionsToLinks with all NA coaccess returns empty", {
  conns <- data.frame(
    Peak1 = c("chr1:1-100", "chr1:200-300"),
    Peak2 = c("chr1:150-200", "chr1:350-400"),
    coaccess = c(NA, NA)
  )
  gi <- ConnectionsToLinks(conns = conns, threshold = 0)
  expect_equal(length(gi), 0)
})

# GetLinkedPeaks / GetLinkedGenes ----------------------------------------------

test_that("GetLinkedPeaks errors when no key", {
  expect_error(
    GetLinkedPeaks(object = atac_small, features = "gene1"),
    regexp = "key"
  )
})

test_that("GetLinkedPeaks errors when no links present", {
  expect_error(
    GetLinkedPeaks(
      object = atac_small, features = "gene1", key = "linkpeaks"
    ),
    regexp = "No links present|Run LinkPeaks"
  )
})

test_that("GetLinkedGenes errors when no links present", {
  expect_error(
    GetLinkedGenes(
      object = atac_small, features = "chr1:1-100", key = "linkpeaks"
    ),
    regexp = "No links present|Run LinkPeaks"
  )
})

test_that("GetLinkedPeaks works on assay with links present", {
  gi <- InteractionSet::GInteractions(
    GRanges("chr1", IRanges(c(100, 200), c(150, 250))),
    GRanges("chr1", IRanges(c(500, 600), c(550, 650)))
  )
  gi$score <- c(0.6, 0.8)
  gi$anchor2.gene_name <- c("geneA", "geneB")
  gi$peak <- c("chr1:100-150", "chr1:200-250")
  obj <- atac_small
  Links(obj) <- list(linkpeaks = gi)
  pks <- GetLinkedPeaks(object = obj, features = "geneA", key = "linkpeaks")
  expect_equal(pks, "chr1:100-150")
})

test_that("GetLinkedPeaks with min.abs.score filter", {
  gi <- InteractionSet::GInteractions(
    GRanges("chr1", IRanges(c(100, 200), c(150, 250))),
    GRanges("chr1", IRanges(c(500, 600), c(550, 650)))
  )
  gi$score <- c(0.1, 0.9)
  gi$anchor2.gene_name <- c("geneA", "geneA")
  gi$peak <- c("chr1:100-150", "chr1:200-250")
  obj <- atac_small
  Links(obj) <- list(linkpeaks = gi)
  res <- GetLinkedPeaks(
    object = obj, features = "geneA", key = "linkpeaks", min.abs.score = 0.5
  )
  expect_equal(res, "chr1:200-250")
})

test_that("GetLinkedGenes works on assay with links present", {
  gi <- InteractionSet::GInteractions(
    GRanges("chr1", IRanges(c(100, 200), c(150, 250))),
    GRanges("chr1", IRanges(c(500, 600), c(550, 650)))
  )
  gi$score <- c(0.6, 0.8)
  gi$anchor2.gene_name <- c("geneA", "geneB")
  obj <- atac_small
  Links(obj) <- list(linkpeaks = gi)
  genes <- GetLinkedGenes(
    object = obj, features = "chr1:100-150", key = "linkpeaks"
  )
  expect_equal(genes, "geneA")
})

test_that("GetLinkedPeaks/Genes Assay5 method errors", {
  a <- SeuratObject::CreateAssay5Object(counts = matrix(0, 2, 2))
  rownames(a) <- c("a", "b")
  colnames(a) <- c("x", "y")
  expect_error(
    GetLinkedPeaks(object = a, features = "a", key = "k"),
    regexp = "ChromatinAssay5"
  )
  expect_error(
    GetLinkedGenes(object = a, features = "a", key = "k"),
    regexp = "ChromatinAssay5"
  )
})

# LinkPeaks --------------------------------------------------------------------

test_that("LinkPeaks errors on non-GRangesAssay", {
  obj <- atac_small
  obj[["normal"]] <- SeuratObject::CreateAssay5Object(
    counts = matrix(0, 3, ncol(atac_small),
      dimnames = list(c("a", "b", "c"), colnames(atac_small)))
  )
  expect_error(
    LinkPeaks(
      object = obj, peak.assay = "normal", expression.assay = "RNA",
      verbose = FALSE
    ),
    regexp = "GRangesAssay"
  )
})

test_that("LinkPeaks errors on bad min.distance", {
  expect_error(
    LinkPeaks(
      object = atac_small, peak.assay = "peaks", expression.assay = "RNA",
      min.distance = "bad", verbose = FALSE
    ),
    regexp = "numeric"
  )
})

test_that("LinkPeaks warns and sets to zero on negative min.distance", {
  expect_warning(
    expect_error(
      LinkPeaks(
        object = atac_small,
        peak.assay = "peaks", expression.assay = "RNA",
        min.distance = -10, verbose = FALSE
      ),
      regexp = "Could not find gene coordinates"
    ),
    regexp = "negative min.distance"
  )
})

test_that("LinkPeaks errors on unknown method", {
  expect_error(
    LinkPeaks(
      object = atac_small, peak.assay = "peaks", expression.assay = "RNA",
      method = "kendall", verbose = FALSE
    ),
    regexp = "pearson|spearman"
  )
})

test_that("LinkPeaks errors when expression layer missing", {
  expect_error(
    LinkPeaks(
      object = atac_small, peak.assay = "peaks", expression.assay = "RNA",
      expression.layer = "wronglayer", verbose = FALSE
    ),
    regexp = "expression layer not found"
  )
})

test_that("LinkPeaks errors with peak.layer not found", {
  # 'count' column must be populated so LinkPeaks skips FindTopFeatures and
  # reaches the explicit peak.layer-not-found check.
  obj <- atac_small
  md <- obj[["peaks"]][[]]
  md$count <- runif(nrow(obj[["peaks"]]), 5, 100)
  obj[["peaks"]][[]] <- md
  expect_error(
    LinkPeaks(
      object = obj, peak.assay = "peaks", expression.assay = "RNA",
      peak.layer = "noexist", verbose = FALSE
    ),
    regexp = "peak layer not found"
  )
})


test_that("LinkPeaks errors when no gene names match expression features", {
  # atac_small annotation contains genes that don't overlap the RNA assay's
  # feature set, so coordinate matching fails after validation passes.
  obj <- atac_small
  expect_error(
    suppressWarnings(LinkPeaks(
      object = obj, peak.assay = "peaks", expression.assay = "RNA",
      min.cells = 0, n_sample = 10, verbose = FALSE
    )),
    regexp = "Could not find gene coordinates"
  )
})

test_that("LinkPeaks errors when no gene coordinates match", {
  obj <- atac_small
  empty_gc <- GRanges()
  mcols(empty_gc)$gene_name <- character(0)
  expect_error(
    suppressWarnings(LinkPeaks(
      object = obj, peak.assay = "peaks", expression.assay = "RNA",
      gene.coords = empty_gc, verbose = FALSE
    )),
    regexp = "Could not find gene coordinates"
  )
})

# DistanceToTSS / LinksToGRanges / LinksToGInteractions -----------------------

test_that("DistanceToTSS returns sparse matrix", {
  peaks <- GRanges("chr1", IRanges(c(100, 5000, 20000), c(200, 5100, 20100)))
  genes <- GRanges(
    "chr1", IRanges(c(1000, 10000), c(2000, 11000)),
    gene_name = c("g1", "g2")
  )
  res <- Signac:::DistanceToTSS(peaks = peaks, genes = genes, distance = 5000)
  expect_s4_class(res, "CsparseMatrix")
  expect_equal(dim(res), c(3, 2))
  # peak 1 (100-200) is within 5kb of g1's TSS only; peak 2 (5000-5100) is
  # within 5kb of both TSSs; peak 3 (20000) is within range of neither
  expect_equal(as.numeric(res[1, ]), c(1, 0))
  expect_equal(as.numeric(res[2, ]), c(1, 1))
  expect_equal(as.numeric(res[3, ]), c(0, 0))
})

test_that("LinksToGRanges converts a link matrix", {
  set.seed(1)
  gene.coords <- GRanges(
    "chr1", IRanges(c(1000, 5000), c(2000, 6000)),
    strand = c("+", "+"), gene_name = c("g1", "g2")
  )
  linkmat <- sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(0.5, 0.8),
    dimnames = list(c("g1", "g2"), c("chr1:1500-1600", "chr1:5500-5600"))
  )
  gr <- Signac:::LinksToGRanges(linkmat = linkmat, gene.coords = gene.coords)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(gr$score, c(0.5, 0.8))
  expect_equal(gr$gene, c("g1", "g2"))
  expect_equal(gr$peak, c("chr1:1500-1600", "chr1:5500-5600"))
})

test_that("LinksToGInteractions converts link matrix", {
  gene.coords <- GRanges(
    "chr1", IRanges(c(1000, 5000), c(2000, 6000)),
    gene_name = c("g1", "g2")
  )
  linkmat <- sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(0.5, 0.8),
    dimnames = list(c("g1", "g2"), c("chr1:1500-1600", "chr1:5500-5600"))
  )
  gi <- Signac:::LinksToGInteractions(
    linkmat = linkmat, gene.coords = gene.coords
  )
  expect_s4_class(gi, "GInteractions")
  expect_equal(gi$score, c(0.5, 0.8))
})

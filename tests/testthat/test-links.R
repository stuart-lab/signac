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
  # no `peak` metadata column: this is the shape LinkPeaks actually stores, and
  # the peak is recovered from the first anchor
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
  obj <- atac_small
  Links(obj) <- list(linkpeaks = gi)
  res <- GetLinkedPeaks(
    object = obj, features = "geneA", key = "linkpeaks", min.abs.score = 0.5
  )
  expect_equal(res, "chr1:200-250")
})

test_that("GetLinkedPeaks reads peaks from links built by LinksToGInteractions", {
  linkmat <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(1, 2), x = c(0.6, 0.8), dims = c(2, 2),
    dimnames = list(c("geneA", "geneB"), c("chr1:100-150", "chr1:200-250"))
  )
  gene.coords <- GRanges(
    "chr1", IRanges(c(500, 600), c(550, 650)), gene_name = c("geneA", "geneB")
  )
  gi <- Signac:::LinksToGInteractions(
    linkmat = linkmat, gene.coords = gene.coords
  )
  expect_null(gi$peak)

  obj <- atac_small
  Links(obj) <- list(linkpeaks = gi)
  pks <- GetLinkedPeaks(object = obj, features = "geneA", key = "linkpeaks")
  expect_equal(pks, "chr1:100-150")
  expect_false(is.null(x = pks))
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
      regexp = "No transcript/gene coordinates found"
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
    regexp = "No transcript/gene coordinates found"
  )
})

test_that("LinkPeaks errors when no gene coordinates match", {
  obj <- atac_small
  empty_gc <- GRanges()
  mcols(empty_gc)$gene_id <- character(0)
  mcols(empty_gc)$gene_name <- character(0)
  expect_error(
    suppressWarnings(LinkPeaks(
      object = obj, peak.assay = "peaks", expression.assay = "RNA",
      gene.coords = empty_gc, verbose = FALSE
    )),
    regexp = "No transcript/gene coordinates found"
  )
})

# FindCandidateLinks / LinkPeaks results ---------------------------------------

# Build an exon-level annotation of the shape GetGRangesFromEnsDb() returns:
# one row per exon, with the transcript recorded in tx_id
MakeExonAnnotation <- function(
  seqnames,
  start,
  end,
  tx_id,
  gene_id,
  strand = "+",
  gene_name = gene_id
) {
  annot <- GRanges(
    seqnames = seqnames,
    ranges = IRanges(start = start, end = end),
    strand = strand
  )
  annot$tx_id <- tx_id
  annot$gene_id <- gene_id
  annot$gene_name <- gene_name
  annot$gene_biotype <- "protein_coding"
  annot$type <- factor(x = "exon", levels = c("cds", "exon", "gap", "utr"))
  return(annot)
}

# Build a peak + RNA object for peak-gene linking tests. The first gene is made
# to correlate with the first peak, so there is a signal for LinkPeaks to find.
MakeLinkObject <- function(annotation, peaks, genes, ncell = 60, seed = 1) {
  set.seed(seed = seed)
  cells <- paste0("cell", seq_len(length.out = ncell))
  counts <- matrix(
    data = rpois(n = length(x = peaks) * ncell, lambda = 3),
    nrow = length(x = peaks),
    dimnames = list(peaks, cells)
  )
  rna <- matrix(
    data = rpois(n = length(x = genes) * ncell, lambda = 5),
    nrow = length(x = genes),
    dimnames = list(genes, cells)
  )
  rna[1, ] <- counts[1, ] * 2 + rpois(n = ncell, lambda = 1)

  assay <- suppressWarnings(
    CreateGRangesAssay(counts = counts, annotation = annotation)
  )
  obj <- SeuratObject::CreateSeuratObject(counts = assay, assay = "peaks")
  obj[["RNA"]] <- suppressWarnings(
    SeuratObject::CreateAssay5Object(counts = rna, data = log1p(x = rna))
  )
  return(obj)
}

MakeLinkTestObject <- function(ncell = 60, npeak = 20, ngene = 4, seed = 42) {
  peak.chr <- rep(x = c("chr1", "chr2"), each = npeak / 2)
  peak.start <- rep(
    x = seq(from = 10000, by = 5000, length.out = npeak / 2), times = 2
  )
  peak.names <- paste0(peak.chr, ":", peak.start, "-", peak.start + 500)

  annot <- GRanges(
    seqnames = rep(x = c("chr1", "chr2"), each = ngene / 2),
    ranges = IRanges(
      start = rep(
        x = seq(from = 12000, by = 20000, length.out = ngene / 2), times = 2
      ),
      width = 3000
    ),
    strand = "+"
  )
  annot$gene_id <- paste0("gene", seq_len(length.out = ngene))
  annot$gene_name <- annot$gene_id
  annot$gene_biotype <- "protein_coding"
  annot$type <- "transcript"
  annot$transcript_id <- paste0("tx", seq_len(length.out = ngene))

  return(MakeLinkObject(
    annotation = annot,
    peaks = peak.names,
    genes = paste0("gene", seq_len(length.out = ngene)),
    ncell = ncell,
    seed = seed
  ))
}

test_that("FindCandidateLinks returns the expected structure", {
  obj <- MakeLinkTestObject()
  cand <- Signac:::FindCandidateLinks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, verbose = FALSE
  )
  expect_named(
    object = cand,
    expected = c(
      "candidate.matrix", "candidate.links", "peaks", "gene.coords",
      "peak.data", "expression.data"
    )
  )
  expect_s4_class(cand$candidate.links, "GInteractions")
  expect_gt(length(cand$candidate.links), 0)
  expect_true(all(cand$candidate.links$tss_distance >= 0))
  # gene.coords must line up with the candidate matrix columns
  expect_equal(
    cand$gene.coords$gene_name, colnames(cand$candidate.matrix)
  )
})

tss.cases <- list(
  list(
    name = "plus strand transcript",
    # TSS at 10000. The final exon starts 45kb away and must not be treated as
    # a TSS, so the peak sitting on it is not a candidate at distance = 10kb
    annotation = MakeExonAnnotation(
      seqnames = "chr1",
      start = c(10000, 30000, 55000),
      end = c(10500, 30500, 60000),
      strand = "+",
      tx_id = "ENST1", gene_id = "ENSG1", gene_name = "geneA"
    ),
    peaks = c("chr1:9000-9500", "chr1:54800-55300",
              "chr2:10000-10500", "chr2:20000-20500"),
    distance = 1e4,
    peak = "chr1:9000-9500",
    tss_distance = 500,
    tss = 10000
  ),
  list(
    name = "minus strand transcript",
    # spans 10000-60000 on the minus strand, so the TSS is at 60000, not 10000
    annotation = MakeExonAnnotation(
      seqnames = "chr1",
      start = c(10000, 55000),
      end = c(10500, 60000),
      strand = "-",
      tx_id = "ENST1", gene_id = "ENSG1", gene_name = "geneA"
    ),
    peaks = c("chr1:9000-9500", "chr1:60500-61000",
              "chr2:10000-10500", "chr2:20000-20500"),
    distance = 1e4,
    peak = "chr1:60500-61000",
    tss_distance = 500,
    tss = 60000
  ),
  list(
    name = "nearest of two transcript TSSs",
    # one gene, two transcripts with TSSs at 10000 and 100000. The peak is
    # 1.5kb from the second, so the pair must resolve to that transcript
    annotation = MakeExonAnnotation(
      seqnames = "chr1",
      start = c(10000, 20000, 100000, 110000),
      end = c(10500, 20500, 100500, 110500),
      strand = "+",
      tx_id = c("ENST1", "ENST1", "ENST2", "ENST2"),
      gene_id = "ENSG1", gene_name = "geneA"
    ),
    peaks = c("chr1:101500-102000", "chr2:10000-10500", "chr2:20000-20500"),
    distance = 5e5,
    peak = "chr1:101500-102000",
    tss_distance = 1500,
    tss = 100000
  )
)

test_that("FindCandidateLinks measures distance to each transcript TSS", {
  for (case in tss.cases) {
    obj <- MakeLinkObject(
      annotation = case$annotation, peaks = case$peaks, genes = "geneA"
    )
    cand <- Signac:::FindCandidateLinks(
      object = obj, peak.assay = "peaks", expression.assay = "RNA",
      min.cells = 0, distance = case$distance, verbose = FALSE
    )
    # a single candidate peak-gene pair, resolved to the nearest TSS
    expect_equal(
      object = cand$candidate.links$peak, expected = case$peak,
      label = paste0(case$name, ": candidate peak")
    )
    expect_equal(
      object = cand$candidate.links$tss_distance, expected = case$tss_distance,
      label = paste0(case$name, ": tss_distance")
    )
    expect_equal(
      object = start(InteractionSet::anchors(cand$candidate.links)$second),
      expected = case$tss,
      label = paste0(case$name, ": TSS used")
    )
  }
})

test_that("FindCandidateLinks does not mix genes sharing a gene symbol", {
  # two gene_ids share the symbol SHARED. GetGeneRanges() uniquifies the names,
  # so the second becomes SHARED.1 and is absent from the expression assay. Its
  # transcripts must not be attributed to SHARED.
  obj <- MakeLinkObject(
    annotation = MakeExonAnnotation(
      seqnames = c("chr1", "chr1", "chr2", "chr2"),
      start = c(10000, 20000, 500000, 510000),
      end = c(10500, 20500, 500500, 510500),
      tx_id = c("txA", "txA", "txB", "txB"),
      gene_id = c("ENSG_A", "ENSG_A", "ENSG_B", "ENSG_B"),
      gene_name = "SHARED"
    ),
    peaks = c("chr1:9000-9500", "chr2:499000-499500", "chr2:600000-600500"),
    genes = c("SHARED", "OTHER")
  )
  cand <- Signac:::FindCandidateLinks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, distance = 1e4, verbose = FALSE
  )
  # only the chr1 peak, next to ENSG_A's TSS
  expect_equal(cand$candidate.links$peak, "chr1:9000-9500")
  # the TSS used must be on the same chromosome as the gene body
  expect_equal(
    as.character(seqnames(InteractionSet::anchors(
      cand$candidate.links)$second)),
    "chr1"
  )
  expect_true(all(as.character(seqnames(cand$gene.coords)) == "chr1"))
})

test_that("FindCandidateLinks excludes multi-chromosome genes entirely", {
  annot <- MakeExonAnnotation(
    seqnames = c("chr1", "chr1", "chrX", "chrX", "chrY", "chrY"),
    start = c(10000, 20000, 1000, 2000, 50000, 51000),
    end = c(10500, 20500, 1200, 2200, 50200, 51200),
    tx_id = c("txA", "txA", "txX", "txX", "txY", "txY"),
    gene_id = c("ENSGA", "ENSGA", "PAR1", "PAR1", "PAR1", "PAR1"),
    gene_name = c("geneA", "geneA", "PAR1", "PAR1", "PAR1", "PAR1")
  )
  peaks <- c("chr1:9000-9500", "chrX:500-900", "chrY:49500-49900",
             "chr1:100000-100500")

  # both naming modes must exclude the gene, and neither may leave an NA in
  # gene.coords, which LinkPeaks indexes into by gene
  for (use.id in c(FALSE, TRUE)) {
    obj <- MakeLinkObject(
      annotation = annot, peaks = peaks,
      genes = if (use.id) c("ENSGA", "PAR1") else c("geneA", "PAR1")
    )
    expect_warning(
      cand <- Signac:::FindCandidateLinks(
        object = obj, peak.assay = "peaks", expression.assay = "RNA",
        min.cells = 0, distance = 1e4, gene.id = use.id, verbose = FALSE
      ),
      regexp = "more than one chromosome"
    )
    expect_equal(
      object = cand$candidate.links$peak, expected = "chr1:9000-9500",
      label = paste0("gene.id = ", use.id, ": candidate peak")
    )
    expect_false(anyNA(as.character(seqnames(cand$gene.coords))))
    expect_true(all(as.character(seqnames(cand$gene.coords)) == "chr1"))
  }
})

test_that("FindCandidateLinks respects min.distance", {
  obj <- MakeLinkTestObject()
  all.cand <- Signac:::FindCandidateLinks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, verbose = FALSE
  )
  filtered <- Signac:::FindCandidateLinks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, min.distance = 10000, verbose = FALSE
  )
  expect_lt(length(filtered$candidate.links), length(all.cand$candidate.links))
  expect_true(all(filtered$candidate.links$tss_distance >= 10000))
})

test_that("LinkPeaks stores links with the expected metadata", {
  obj <- MakeLinkTestObject()
  obj <- LinkPeaks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, verbose = FALSE
  )
  lnk <- Links(object = obj[["peaks"]])[["linkpeaks"]]
  expect_s4_class(lnk, "GInteractions")
  expect_gt(length(lnk), 0)
  expect_equal(
    colnames(S4Vectors::mcols(lnk)),
    c(
      "peak", "gene", "tss_distance", "frac_peak_in_gene", "closest_gene",
      "score", "zscore", "pvalue"
    )
  )
  expect_type(lnk$closest_gene, "logical")
  # zscore is not computed by default
  expect_true(all(is.na(lnk$zscore)))
  expect_true(all(is.na(lnk$pvalue)))
  expect_true(all(abs(lnk$score) > 0.05))
  # gene1 was constructed to correlate with the first peak
  expect_true("chr1:10000-10500" %in% lnk$peak[lnk$gene == "gene1"])
})

test_that("LinkPeaks computes z-scores when zscore = TRUE", {
  obj <- MakeLinkTestObject()
  md <- obj[["peaks"]][[]]
  md$GC.percent <- runif(n = nrow(x = md), min = 30, max = 60)
  md$sequence.length <- seq_len(length.out = nrow(x = md)) + 500
  obj[["peaks"]][[]] <- md
  obj <- LinkPeaks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, zscore = TRUE, n_sample = 10,
    pvalue_cutoff = 0.9, verbose = FALSE
  )
  lnk <- Links(object = obj[["peaks"]])[["linkpeaks"]]
  expect_gt(length(lnk), 0)
  expect_false(any(is.na(lnk$zscore)))
  expect_true(all(lnk$pvalue < 0.9))
})

test_that("LinkPeaks handles zero-variance peaks and genes", {
  obj <- MakeLinkTestObject()
  cnt <- SeuratObject::LayerData(obj[["peaks"]], layer = "counts")
  # make the first peak identical in every cell
  cnt[1, ] <- 5
  SeuratObject::LayerData(obj[["peaks"]], layer = "counts") <- cnt

  expect_no_error(
    obj <- suppressWarnings(LinkPeaks(
      object = obj, peak.assay = "peaks", expression.assay = "RNA",
      min.cells = 0, verbose = FALSE
    ))
  )
  lnk <- Links(object = obj[["peaks"]])[["linkpeaks"]]
  expect_false(any(is.na(lnk$peak)))
  expect_false(any(is.na(lnk$score)))
  expect_false(rownames(cnt)[[1]] %in% lnk$peak)
})

test_that("LinkPeaks stores empty links when nothing passes cor.cutoff", {
  obj <- MakeLinkTestObject()
  obj <- LinkPeaks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, cor.cutoff = 1, verbose = FALSE
  )
  lnk <- Links(object = obj[["peaks"]])[["linkpeaks"]]
  expect_s4_class(lnk, "GInteractions")
  expect_equal(length(lnk), 0)
  # empty links keep the same schema as populated links
  expect_true("gene" %in% colnames(S4Vectors::mcols(lnk)))
})

test_that("GetLinkedPeaks/Genes work on LinkPeaks output", {
  obj <- MakeLinkTestObject()
  obj <- LinkPeaks(
    object = obj, peak.assay = "peaks", expression.assay = "RNA",
    min.cells = 0, verbose = FALSE
  )
  lnk <- Links(object = obj[["peaks"]])[["linkpeaks"]]
  expected.peaks <- unique(x = lnk$peak[lnk$gene == "gene1"])
  pks <- GetLinkedPeaks(
    object = obj[["peaks"]], features = "gene1", key = "linkpeaks"
  )
  expect_setequal(pks, expected.peaks)
  gns <- GetLinkedGenes(
    object = obj[["peaks"]], features = pks[[1]], key = "linkpeaks"
  )
  expect_true("gene1" %in% gns)
})

# LinksToGInteractions ---------------------------------------------------------

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

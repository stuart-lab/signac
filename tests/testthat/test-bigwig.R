test_that("CreateBWGroup works with single tile", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "createBW")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 249250621),
    tiles = GRanges(
      seqnames = "chr1", ranges = IRanges(start = 1, end = 249250621)
    ),
    normBy = NULL,
    tileSize = 249250621,
    normMethod = "RC",
    cutoff = NULL,
    outdir = outdir
  )
  expect_equal(object = length(list.files(outdir)), expected = 2)
  expect(
    file.exists(file.path(outdir, "0-TileSize-249250621-normMethod-rc.bw")),
    "File does not exist."
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "0-TileSize-249250621-normMethod-rc.bw")
  )
  expect_equal(object = bw$score, 20000)
})

test_that("CreateBWGroup ncells uses the supplied per-group cell count", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "createBW_ncells")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  # the bed has a single unique barcode, but nCells says the group has 5 cells;
  # ncells normalization must divide by 5, not by 1
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 249250621),
    tiles = GRanges(
      seqnames = "chr1", ranges = IRanges(start = 1, end = 249250621)
    ),
    normBy = NULL,
    nCells = c("0" = 5),
    tileSize = 249250621,
    normMethod = "ncells",
    cutoff = NULL,
    outdir = outdir
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "0-TileSize-249250621-normMethod-ncells.bw")
  )
  # 5 fragments -> 10 insertion events in the single tile, divided by 5 cells
  expect_equal(object = bw$score, 2)
})

test_that("CreateBWGroup works with 100bp tile", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "createBW2")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 249250621),
    tiles = GRanges(
      seqnames = "chr1",
      ranges = IRanges(
        start = seq(1, 249250621, 100),
        end = c(seq(100, 249250621, 100), 249250621)
      )
    ),
    normBy = NULL,
    tileSize = 100,
    normMethod = "RC",
    cutoff = NULL,
    outdir = outdir
  )
  expect_equal(object = length(list.files(outdir)), expected = 2)
  expect(
    file.exists(file.path(outdir, "0-TileSize-100-normMethod-rc.bw")),
    "File does not exist."
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "0-TileSize-100-normMethod-rc.bw")
  )
  expect_equal(object = bw$score, c(6000, 8000, 2000, 0))
})

test_that("CreateBWGroup works with seqlength equal to final pos", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "createBW3")
  dir.create(outdir, showWarnings = FALSE)
  fake.bed.data <- data.frame(
    seqnames = rep("chr1", 5),
    start = c(0, 10, 100, 110, 300),
    end = c(100, 150, 200, 250, 500),
    cell_name = rep("fake_cell", 5),
    nb = 1:5
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  CreateBWGroup(
    groupNamei = "0",
    availableChr = "chr1",
    chromLengths = c("chr1" = 500),
    tiles = GRanges(
      seqnames = "chr1",
      ranges = IRanges(start = seq(1, 499, 100), end = c(seq(100, 500, 100)))
    ),
    normBy = NULL,
    tileSize = 100,
    normMethod = "RC",
    cutoff = NULL,
    outdir = outdir
  )
  expect_equal(object = length(list.files(outdir)), expected = 2)
  expect(
    file.exists(file.path(outdir, "0-TileSize-100-normMethod-rc.bw")),
    "File does not exist."
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "0-TileSize-100-normMethod-rc.bw")
  )
  expect_equal(object = bw$score, c(6000, 8000, 2000))
})

test_that("ExportBigwig works", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "ExportBigwig")
  if (dir.exists(paths = outdir)) {
    unlink(x = outdir, recursive = TRUE)
  }
  dir.create(outdir, showWarnings = FALSE)
  bdir <- file.path(tempdir(), "ExportBigwig_bed")
  if (dir.exists(paths = bdir)) {
    unlink(x = bdir, recursive = TRUE)
  }
  dir.create(bdir, showWarnings = FALSE)
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- cells
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate.fragments = FALSE
  )
  Fragments(atac_small) <- frags
  # split into two groups, each above minCells
  groups <- rep(x = c("g1", "g2"), length.out = ncol(x = atac_small))
  atac_small <- AddMetaData(
    object = atac_small, metadata = groups, col.name = "bw_group"
  )
  # v2 objects do not store chromosome lengths, so supply them explicitly
  covfiles <- ExportBigwig(
    object = atac_small,
    group.by = "bw_group",
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    seqlengths = c("chr1" = 1e6),
    outdir = outdir,
    temp.dir = bdir,
    verbose = FALSE
  )
  # one bigwig file per group written to outdir
  expect_length(object = covfiles, n = 2)
  expect(
    file.exists(file.path(outdir, "g1-TileSize-100-normMethod-rc.bw")),
    "File does not exist."
  )
  expect(
    file.exists(file.path(outdir, "g2-TileSize-100-normMethod-rc.bw")),
    "File does not exist."
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "g1-TileSize-100-normMethod-rc.bw")
  )
  expect_equal(object = length(seqlengths(bw)), expected = 1)
  expect_equal(object = unname(seqlengths(bw)), expected = 1e6)
  # intermediate bed files go to temp.dir and are cleaned up by default
  expect_false(object = file.exists(file.path(bdir, "g1.bed")))
  expect_false(object = file.exists(file.path(bdir, "g2.bed")))
  # no bed files left behind in the output directory
  expect_length(object = list.files(outdir, pattern = "[.]bed$"), n = 0)
})

test_that("ExportBigwig handles group names containing spaces", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "ExportBigwig_space")
  if (dir.exists(paths = outdir)) {
    unlink(x = outdir, recursive = TRUE)
  }
  dir.create(outdir, showWarnings = FALSE)
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- cells
  frags <- CreateFragmentObject(
    path = fpath, cells = cells, verbose = FALSE, validate.fragments = FALSE
  )
  Fragments(atac_small) <- frags
  # group labels with spaces must be sanitized consistently so the bed files
  # written by SplitFragments are found again
  groups <- rep(x = c("group a", "group b"), length.out = ncol(x = atac_small))
  atac_small <- AddMetaData(
    object = atac_small, metadata = groups, col.name = "bw_group"
  )
  expect_no_error(
    ExportBigwig(
      object = atac_small,
      group.by = "bw_group",
      seqlengths = c("chr1" = 1e6),
      outdir = outdir,
      verbose = FALSE
    )
  )
  expect(
    file.exists(file.path(outdir, "group_a-TileSize-100-normMethod-rc.bw")),
    "File does not exist."
  )
})

test_that("CreateBWGroup clamps fragments beyond the chromosome end", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "createBW_clamp")
  dir.create(outdir, showWarnings = FALSE)
  # a fragment whose end (250) extends past the chromosome length (200)
  fake.bed.data <- data.frame(
    seqnames = "chr1",
    start = 0,
    end = 250,
    cell_name = "fake_cell",
    nb = 1
  )
  write.table(
    fake.bed.data, file.path(outdir, "0.bed"),
    col.names = FALSE, quote = FALSE, sep = "\t",
    row.names = FALSE
  )
  # should not error despite the out-of-bounds fragment end
  expect_no_error(
    CreateBWGroup(
      groupNamei = "0",
      availableChr = "chr1",
      chromLengths = c("chr1" = 200),
      tiles = GRanges(
        seqnames = "chr1", ranges = IRanges(start = c(1, 101), end = c(100, 200))
      ),
      normBy = NULL,
      tileSize = 100,
      normMethod = "RC",
      cutoff = NULL,
      outdir = outdir
    )
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "0-TileSize-100-normMethod-rc.bw")
  )
  # both ends counted, clamped within [1, 200]: 1 event in tile1, 1 in tile2,
  # RC-normalized by the single fragment -> 10000 across the whole chromosome
  expect_equal(object = bw$score, 10000)
})

test_that("ExportBigwig supports NULL and metadata-column normalization", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "ExportBigwig_norm")
  if (dir.exists(paths = outdir)) {
    unlink(x = outdir, recursive = TRUE)
  }
  dir.create(outdir, showWarnings = FALSE)
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- cells
  frags <- CreateFragmentObject(
    path = fpath, cells = cells, verbose = FALSE, validate.fragments = FALSE
  )
  Fragments(atac_small) <- frags
  atac_small$depth <- seq_len(length.out = ncol(x = atac_small))

  # NULL normMethod must not error and is treated as "none"
  expect_no_error(
    ExportBigwig(
      object = atac_small,
      normMethod = NULL,
      seqlengths = c("chr1" = 1e6),
      outdir = outdir,
      verbose = FALSE
    )
  )
  expect(
    file.exists(
      file.path(outdir, "SeuratProject-TileSize-100-normMethod-none.bw")
    ),
    "File does not exist."
  )

  # metadata-column normalization must produce finite (non-NaN) scores
  ExportBigwig(
    object = atac_small,
    normMethod = "depth",
    seqlengths = c("chr1" = 1e6),
    outdir = outdir,
    verbose = FALSE
  )
  bw <- rtracklayer::import.bw(
    file.path(outdir, "SeuratProject-TileSize-100-normMethod-depth.bw")
  )
  expect_true(all(is.finite(bw$score)))
  expect_true(any(bw$score > 0))

  # a non-numeric normalization column gives a clear error
  atac_small$celltype <- rep(x = c("a", "b"), length.out = ncol(x = atac_small))
  expect_error(
    object = ExportBigwig(
      object = atac_small,
      normMethod = "celltype",
      seqlengths = c("chr1" = 1e6),
      outdir = outdir,
      verbose = FALSE
    ),
    regexp = "must be numeric"
  )
})

test_that("ExportBigwig warns when no group meets minCells", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "ExportBigwig_mincells")
  dir.create(outdir, showWarnings = FALSE)
  fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- cells
  frags <- CreateFragmentObject(
    path = fpath, cells = cells, verbose = FALSE, validate.fragments = FALSE
  )
  Fragments(atac_small) <- frags
  expect_warning(
    object = res <- ExportBigwig(
      object = atac_small,
      minCells = 1e6,
      seqlengths = c("chr1" = 1e6),
      outdir = outdir,
      verbose = FALSE
    ),
    regexp = "minCells"
  )
  expect_length(object = res, n = 0)
})

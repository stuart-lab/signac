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

test_that("ExportGroupBW works", {
  skip_if_not_installed("rtracklayer")
  outdir <- file.path(tempdir(), "ExportGroupBW")
  if (dir.exists(paths = outdir)) {
    unlink(x = outdir, recursive = TRUE)
  }
  dir.create(outdir, showWarnings = FALSE)
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
  covfiles <- ExportGroupBW(
    object = atac_small,
    group.by = "bw_group",
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    seqlengths = c("chr1" = 1e6),
    outdir = outdir,
    verbose = FALSE
  )
  # one bed file + one bigwig file per group
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
})

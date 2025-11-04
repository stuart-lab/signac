suppressWarnings(RNGversion(vstr = "3.5.3"))

count_vec <- c(3, 3, 4, 3, 5, 5, 3, 3, 3, 2, 1, 3, 4, 1, 1, 2, 2, 2, 3, 5, 
  5, 6, 2, 1, 3, 5, 2, 3, 5, 8, 4, 5, 1, 3, 3, 3, 1, 4, 2, 2, 4, 
  1, 2, 1, 2, 1, 2, 5, 1, 2, 1, 4, 2, 5, 3, 1, 2, 2, 2, 3, 3, 3, 
  1, 1, 2, 3, 1, 2, 3, 1, 2, 1, 1, 1, 1, 1, 1)

mononucleosome <- c(0, 1, 0, 1, 4, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 
                    1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 
                    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 1, 0, 1, 
                    0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0)

nucleosome_free <- c(3, 2, 4, 2, 1, 5, 2, 2, 2, 2, 1, 3, 2, 1, 1, 1, 2, 1, 1, 4, 
                     4, 5, 0, 0, 2, 3, 1, 2, 3, 5, 3, 5, 1, 3, 2, 2, 1, 4, 1, 1, 3, 
                     0, 2, 1, 1, 1, 2, 5, 1, 2, 1, 3, 1, 3, 3, 1, 2, 2, 2, 2, 2, 2, 
                     1, 1, 0, 2, 1, 2, 3, 0, 1, 1, 1, 1, 1, 1, 1)

read_vec <- c(5, 4, 6, 4, 10, 8, 4, 7, 12, 2, 5, 5, 5, 5, 2, 6, 2, 5, 6, 
              12, 10, 16, 6, 2, 4, 5, 5, 8, 10, 18, 5, 7, 2, 5, 3, 4, 2, 7, 
              7, 2, 4, 1, 3, 4, 8, 1, 2, 10, 1, 4, 3, 6, 4, 11, 7, 2, 8, 3, 
              7, 10, 3, 5, 2, 3, 4, 4, 1, 7, 8, 2, 2, 1, 2, 1, 4, 2, 1)


fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
cells <- colnames(x = atac_small)
names(x = cells) <- paste0("test_", cells)
frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
frags_header <- CreateFragmentObject(path = fpath_headered, cells = cells, verbose = FALSE, tolerance = 0.5)

test_that("Fragment head method works", {
  expect_equal(
    object = head(x = frags),
    expected = structure(list(chrom = c("chr1", "chr1", "chr1", "chr1", "chr1", 
                                        "chr1"),
                              start = c(10168L, 77156L, 181068L, 181348L, 190824L, 267990L),
                              end = c(10216L, 77221L, 181131L, 181514L, 191033L, 268023L),
                              barcode = c("AAACTGCTCTTCCACG-1", "AAAGATGCACGCTGTG-1",
                                          "AAACTGCCAGTAGGCA-1", "AAACTGCAGGCAAGCT-1",
                                          "AAAGATGAGACTAGGC-1", "AAACTCGAGAGGCCTA-1"),
                              readCount = c(2L, 1L, 1L, 1L, 3L, 1L)),
                         class = "data.frame", row.names = c(NA, -6L))
  )
  expect_equal(
    object = head(x = frags_header),
    expected = structure(list(chrom = c("chr1", "chr1", "chr1", "chr1", "chr1", 
                                        "chr1"),
                              start = c(10168L, 77156L, 181068L, 181348L, 190824L, 267990L),
                              end = c(10216L, 77221L, 181131L, 181514L, 191033L, 268023L),
                              barcode = c("AAACTGCTCTTCCACG-1", "AAAGATGCACGCTGTG-1",
                                          "AAACTGCCAGTAGGCA-1", "AAACTGCAGGCAAGCT-1",
                                          "AAAGATGAGACTAGGC-1", "AAACTCGAGAGGCCTA-1"),
                              readCount = c(2L, 1L, 1L, 1L, 3L, 1L)),
                         class = "data.frame", row.names = c(NA, -6L))
  )
})

test_that("Subset fragment object works", {
  cells.use <- c("test_AAAGATGAGCGAGCTA-1", "test_AAACTGCTCCAGGTCG-1")
  subs <- subset(frags, cells = cells.use)
  expect_equal(
    object = Cells(subs),
    expected = c("test_AAACTGCTCCAGGTCG-1", "test_AAAGATGAGCGAGCTA-1")
  )
})

test_that("UpdatePath works", {
  expect_error(
    object = UpdatePath(frags, new.path = "x")
  )
})

test_that("CountFragments works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
  counts <- CountFragments(fragments = fpath)
  counts_headered <- CountFragments(fragments = fpath_headered)

  expect_equal(
    object = counts_headered,
    expected = counts
  )
  expect_equal(
    object = counts$frequency_count,
    expected = count_vec
  )
  expect_equal(
    object = counts$mononucleosomal,
    expected = mononucleosome
  )

  expect_equal(
    object = counts$nucleosome_free,
    expected = nucleosome_free
  )

  expect_equal(
    object = counts$reads_count,
    expected = read_vec
  )
})

test_that("ExtractFragments works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- paste0("test_", cells)
  frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
  counts <- ExtractFragments(fragments = frags, verbose = FALSE)
  frags_headered <- CreateFragmentObject(path = fpath_headered, cells = cells, verbose = FALSE, tolerance = 0.5)
  counts_headered <- ExtractFragments(fragments = frags_headered, verbose = FALSE)

  expect_equal(
    object = counts_headered,
    expected = counts
  )

  expect_equal(
    object = counts$CB,
    expected = names(x = cells)
  )

  expect_equal(
    object = head(x = counts$frequency_count),
    expected = c(1, 0, 2, 3, 3, 5)
  )
  expect_equal(
    object = head(x = counts$mononucleosomal),
    expected = c(0, 0, 0, 1, 1, 0)
  )

  expect_equal(
    object = head(x = counts$nucleosome_free),
    expected = c(1, 0, 2, 2, 2, 5)
  )

  expect_equal(
    object = head(x = counts$reads_count),
    expected = c(1, 0, 3, 10, 8, 10)
  )
})

test_that("ValidateCells works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- paste0("test_", cells)
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  frags_headered <- CreateFragmentObject(
    path = fpath_headered,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  valid <- Signac:::ValidateCells(
    object = frags,
    verbose = FALSE,
    tolerance = 0.5
  )
  invalid <- Signac:::ValidateCells(
    object = frags,
    verbose = FALSE,
    tolerance = 0
  )
  expect_true(object = valid)
  expect_false(object = invalid)

  valid_h <- Signac:::ValidateCells(
    object = frags_headered,
    verbose = FALSE,
    tolerance = 0.5
  )
  invalid_h <- Signac:::ValidateCells(
    object = frags_headered,
    verbose = FALSE,
    tolerance = 0
  )
  expect_true(object = valid_h)
  expect_false(object = invalid_h)
})

test_that("ValidateHash works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- paste0("test_", cells)
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  frags_headered <- CreateFragmentObject(
    path = fpath_headered,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )

  valid <- Signac:::ValidateHash(
    object = frags,
    verbose = FALSE
  )
  expect_true(object = valid)

  valid_h <- Signac:::ValidateHash(
    object = frags_headered,
    verbose = FALSE
  )
  expect_true(object = valid_h)
})

test_that("FilterCells works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
  tmpf <- tempfile(fileext = ".gz")
  FilterCells(
    fragments = fpath,
    cells = head(colnames(atac_small)),
    outfile = tmpf
  )
  output <- read.table(file = tmpf, stringsAsFactors = FALSE)

  file.remove(tmpf)
  FilterCells(
    fragments = fpath_headered,
    cells = head(colnames(atac_small)),
    outfile = tmpf
  )
  output_headered <- read.table(file = tmpf, stringsAsFactors = FALSE)

  expect_equal(
    object = output,
    expected = output_headered
  )

  expected <- structure(list(V1 = c("chr1", "chr1", "chr1", "chr1", "chr1", 
                                    "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                    "chr1"),
                             V2 = c(778433L, 778618L, 778663L, 778719L, 778737L, 778741L, 778808L, 778853L, 778869L, 813095L, 821919L, 827477L, 827489L, 827501L),
                             V3 = c(778658L, 778663L, 778796L, 778810L, 778781L, 778808L, 778844L, 778932L, 778937L, 813174L, 822085L, 827596L, 827614L, 827560L),
                             V4 = c("AAACGAACAAGCACTT-1", "AAACGAACAAGCGGTA-1", "AAACGAACAAGCGGTA-1", "AAACGAAAGGAAGAAC-1", "AAACGAAAGTCGACCC-1", "AAACGAACAAGCGGTA-1", "AAACGAACAAGCACTT-1", "AAACGAAAGGAAGAAC-1", "AAACGAACAAGCGGTA-1", "AAACGAACAAGCGGTA-1", "AAACGAAAGTCGACCC-1", "AAACGAAAGTCGACCC-1", "AAACGAACAAGCACTT-1", "AAACGAAAGAGAGGTA-1"),
                             V5 = c(1L, 3L, 2L, 1L, 3L, 1L, 5L, 2L, 1L, 3L, 1L, 6L, 2L, 1L)),
                        class = "data.frame", row.names = c(NA, -14L))
  expect_equal(object = output, expected = expected)
})

test_that("SplitFragments works", {
  fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
  fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
  cells <- colnames(x = atac_small)
  names(x = cells) <- cells
  frags <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  frags_headered <- CreateFragmentObject(
    path = fpath,
    cells = cells,
    verbose = FALSE,
    validate = FALSE
  )
  Fragments(atac_small) <- frags
  SplitFragments(
    object = atac_small,
    assay = "peaks",
    group.by = "cluster",
    outdir = tempdir()
  )
  of1 <- paste0(tempdir(), .Platform$file.sep, "1.bed")
  of2 <- paste0(tempdir(), .Platform$file.sep, "2.bed")

  bed1 <- read.table(file = of1, sep = "\t", stringsAsFactors = FALSE)
  bed2 <- read.table(file = of2, sep = "\t", stringsAsFactors = FALSE)

  expect_equal(
    object = bed1,
    expected = structure(list(V1 = c("chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1"), V2 = c(267990L, 586175L, 633984L, 633994L, 634004L, 
                                                             634005L, 634018L, 775171L, 777689L, 778403L, 778419L, 778421L, 
                                                             778433L, 778436L, 778558L, 778564L, 778575L, 778587L, 778587L, 
                                                             778592L, 778599L, 778612L, 778617L, 778618L, 778630L, 778637L, 
                                                             778637L, 778637L, 778658L, 778660L, 778660L, 778660L, 778660L, 
                                                             778661L, 778662L, 778663L, 778663L, 778680L, 778719L, 778731L, 
                                                             778736L, 778741L, 778745L, 778745L, 778747L, 778751L, 778762L, 
                                                             778768L, 778781L, 778781L, 778802L, 778808L, 778813L, 778827L, 
                                                             778836L, 778853L, 778869L, 778869L, 778874L, 778902L, 778919L, 
                                                             778927L, 779230L, 813095L, 817147L, 817247L, 817278L, 817298L, 
                                                             817308L, 817312L, 817323L, 817335L, 818777L, 821059L, 822568L, 
                                                             823344L, 825857L, 826702L, 827080L, 827282L, 827293L, 827403L, 
                                                             827408L, 827413L, 827418L, 827431L, 827484L, 827489L, 827501L, 
                                                             827511L, 827543L, 827556L, 827560L, 827566L, 827571L), V3 = c(268023L, 
                                                                                                                           586210L, 634051L, 634068L, 634051L, 634038L, 634063L, 775241L, 
                                                                                                                           777881L, 778612L, 778613L, 778658L, 778658L, 778557L, 778618L, 
                                                                                                                           778636L, 778755L, 778612L, 778660L, 778635L, 778660L, 778684L, 
                                                                                                                           778737L, 778663L, 778735L, 778668L, 778735L, 778755L, 778745L, 
                                                                                                                           778728L, 778741L, 778761L, 778762L, 778762L, 778771L, 778796L, 
                                                                                                                           778844L, 778755L, 778810L, 778781L, 778765L, 778808L, 778772L, 
                                                                                                                           778779L, 778781L, 778808L, 778808L, 778910L, 778829L, 778848L, 
                                                                                                                           778831L, 778844L, 778963L, 778927L, 778938L, 778932L, 778937L, 
                                                                                                                           779099L, 778910L, 778943L, 778949L, 778960L, 779384L, 813174L, 
                                                                                                                           817209L, 817379L, 817418L, 817418L, 817508L, 817440L, 817368L, 
                                                                                                                           817502L, 819129L, 821302L, 822796L, 823502L, 825905L, 826789L, 
                                                                                                                           827575L, 827318L, 827348L, 827439L, 827613L, 827652L, 827505L, 
                                                                                                                           827547L, 827540L, 827614L, 827560L, 827558L, 827569L, 827607L, 
                                                                                                                           827602L, 827602L, 827801L), V4 = c("AAACTCGAGAGGCCTA-1", "AAAGGATCAATACTGC-1", 
                                                                                                                                                              "AAAGATGAGATAGGTT-1", "AAACTCGTCTTCCACG-1", "AAACTGCGTCCGCTTT-1", 
                                                                                                                                                              "AAAGATGCACCGTTGG-1", "AAACTGCAGATCTAAG-1", "AAACGAAGTAAACCCT-1", 
                                                                                                                                                              "AAAGATGTCGTTCCGT-1", "AAACTCGGTAGAAGCC-1", "AAAGATGCACCGAAAG-1", 
                                                                                                                                                              "AAACTGCAGATCTAAG-1", "AAACGAACAAGCACTT-1", "AAACTCGTCCTCAAGA-1", 
                                                                                                                                                              "AAACTCGTCACTAGGT-1", "AAACTCGGTAGAAGCC-1", "AAACTGCAGAACGTCG-1", 
                                                                                                                                                              "AAACGAACACCGTTGG-1", "AAACTCGTCAGGTCTA-1", "AAAGGATCAACGTCCG-1", 
                                                                                                                                                              "AAAGGATCAACGTCCG-1", "AAACTGCGTAGGTGAC-1", "AAAGATGTCTGGCGCA-1", 
                                                                                                                                                              "AAACGAACAAGCGGTA-1", "AAACTCGTCTTCCACG-1", "AAAGATGGTTTGCATG-1", 
                                                                                                                                                              "AAACTCGGTAGAAGCC-1", "AAAGGATAGTCGGGAT-1", "AAACTCGCAAGCACTT-1", 
                                                                                                                                                              "AAACTCGGTAGAAGCC-1", "AAACTGCTCCTTTGAT-1", "AAACTGCGTTCTCGAA-1", 
                                                                                                                                                              "AAACTGCAGAACGTCG-1", "AAAGATGTCAAGAGGC-1", "AAAGATGGTCCGTCGA-1", 
                                                                                                                                                              "AAACGAACAAGCGGTA-1", "AAACTGCGTTCTCGAA-1", "AAACTGCGTCGGCTGT-1", 
                                                                                                                                                              "AAACGAAAGGAAGAAC-1", "AAAGATGGTAGACGCA-1", "AAACTCGAGAGGCCTA-1", 
                                                                                                                                                              "AAACGAACAAGCGGTA-1", "AAACTCGTCCTCAAGA-1", "AAACGAACACCGTTGG-1", 
                                                                                                                                                              "AAAGATGTCAAGAGGC-1", "AAACTCGCAAGTGGCA-1", "AAACTGCAGAACGTCG-1", 
                                                                                                                                                              "AAACTCGTCTTCCACG-1", "AAACTCGGTAGAAGCC-1", "AAAGGATGTACATGGG-1", 
                                                                                                                                                              "AAAGGATCAACGTCCG-1", "AAACGAACAAGCACTT-1", "AAAGATGTCAAGAGGC-1", 
                                                                                                                                                              "AAACTGCTCCTTTGAT-1", "AAACTGCGTTAGAGAT-1", "AAACGAAAGGAAGAAC-1", 
                                                                                                                                                              "AAACGAACAAGCGGTA-1", "AAAGGATCAACTCGTA-1", "AAAGATGTCTGGCGCA-1", 
                                                                                                                                                              "AAAGGATCAATACTGC-1", "AAAGATGGTCCGTCGA-1", "AAACTGCTCCTTTGAT-1", 
                                                                                                                                                              "AAACTCGTCCTCAAGA-1", "AAACGAACAAGCGGTA-1", "AAACTGCTCAAGGCAG-1", 
                                                                                                                                                              "AAACTGCTCAAGGCAG-1", "AAACTCGAGAGGCCTA-1", "AAAGATGTCAAGAGGC-1", 
                                                                                                                                                              "AAAGGATCAATACTGC-1", "AAAGATGTCCACGCTT-1", "AAAGATGTCCACGCTT-1", 
                                                                                                                                                              "AAAGATGTCAAGAGGC-1", "AAACGAACACCGTTGG-1", "AAACTGCGTAACCCAT-1", 
                                                                                                                                                              "AAACTCGTCCTCAAGA-1", "AAAGGATCAACTCGTA-1", "AAAGATGGTCCGTCGA-1", 
                                                                                                                                                              "AAAGGATCAACGTCCG-1", "AAACTCGCAAGTGGCA-1", "AAACTCGTCCTCAAGA-1", 
                                                                                                                                                              "AAACGAATCGTGCTGG-1", "AAAGATGTCGTTCCGT-1", "AAACTGCTCCTTTGAT-1", 
                                                                                                                                                              "AAACGAATCGTGCTGG-1", "AAACTCGCAAGTGGCA-1", "AAAGATGAGAACTCCT-1", 
                                                                                                                                                              "AAACTCGAGAGGCCTA-1", "AAACGAACAAGCACTT-1", "AAACGAAAGAGAGGTA-1", 
                                                                                                                                                              "AAACTCGAGAGGCCTA-1", "AAAGATGGTAGACGCA-1", "AAAGGATAGTCGGGAT-1", 
                                                                                                                                                              "AAACTCGGTAGAAGCC-1", "AAACGAAGTAAACCCT-1", "AAACTGCTCAAGGCAG-1"
                                                                                                                           ), V5 = c(1L, 2L, 5L, 1L, 5L, 2L, 3L, 1L, 4L, 3L, 2L, 3L, 1L, 
                                                                                                                                     4L, 2L, 3L, 1L, 1L, 2L, 4L, 1L, 1L, 1L, 3L, 2L, 1L, 5L, 2L, 3L, 
                                                                                                                                     2L, 2L, 3L, 1L, 3L, 2L, 2L, 1L, 2L, 1L, 6L, 1L, 1L, 1L, 1L, 3L, 
                                                                                                                                     1L, 1L, 2L, 2L, 2L, 1L, 5L, 2L, 1L, 3L, 2L, 1L, 2L, 1L, 1L, 1L, 
                                                                                                                                     1L, 3L, 3L, 1L, 2L, 2L, 2L, 1L, 5L, 2L, 1L, 2L, 2L, 1L, 2L, 4L, 
                                                                                                                                     1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 
                                                                                                                                     1L, 1L)), class = "data.frame", row.names = c(NA, -95L))
  )

  expect_equal(
    object = bed2,
    expected = structure(list(V1 = c("chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", 
                                     "chr1", "chr1", "chr1", "chr1"), V2 = c(10168L, 77156L, 181068L, 
                                                                             181348L, 190824L, 586180L, 586197L, 605580L, 629910L, 634004L, 
                                                                             634023L, 778274L, 778309L, 778357L, 778385L, 778400L, 778403L, 
                                                                             778414L, 778419L, 778420L, 778423L, 778424L, 778436L, 778446L, 
                                                                             778446L, 778451L, 778565L, 778582L, 778593L, 778593L, 778602L, 
                                                                             778602L, 778608L, 778612L, 778612L, 778615L, 778624L, 778626L, 
                                                                             778637L, 778637L, 778655L, 778658L, 778659L, 778660L, 778660L, 
                                                                             778660L, 778662L, 778663L, 778668L, 778678L, 778695L, 778696L, 
                                                                             778712L, 778719L, 778737L, 778745L, 778745L, 778745L, 778755L, 
                                                                             778758L, 778768L, 778784L, 778806L, 778853L, 778869L, 778888L, 
                                                                             778889L, 778895L, 778907L, 778910L, 811333L, 817069L, 817264L, 
                                                                             817265L, 817330L, 817332L, 817344L, 820527L, 820681L, 821919L, 
                                                                             826722L, 827070L, 827275L, 827278L, 827280L, 827300L, 827306L, 
                                                                             827306L, 827360L, 827405L, 827428L, 827477L, 827500L, 827500L, 
                                                                             827501L, 827512L, 827516L, 827528L, 827543L, 827568L, 827568L, 
                                                                             827579L, 827581L, 827591L, 827719L), V3 = c(10216L, 77221L, 181131L, 
                                                                                                                         181514L, 191033L, 586231L, 586234L, 605805L, 629978L, 634039L, 
                                                                                                                         634052L, 778637L, 778637L, 778660L, 778582L, 778592L, 778612L, 
                                                                                                                         778593L, 778667L, 778654L, 778613L, 778637L, 778613L, 778609L, 
                                                                                                                         778637L, 778592L, 778602L, 778624L, 778644L, 778719L, 778712L, 
                                                                                                                         778781L, 778668L, 778655L, 778674L, 778655L, 778662L, 778660L, 
                                                                                                                         778661L, 778784L, 778745L, 778726L, 778762L, 778735L, 778745L, 
                                                                                                                         778766L, 778712L, 778755L, 778750L, 778747L, 778764L, 778808L, 
                                                                                                                         778758L, 778868L, 778781L, 778776L, 778781L, 778810L, 778827L, 
                                                                                                                         778853L, 778827L, 778823L, 778850L, 778948L, 778942L, 778948L, 
                                                                                                                         778949L, 779097L, 778948L, 778965L, 811504L, 817349L, 817344L, 
                                                                                                                         817349L, 817377L, 817408L, 817478L, 820684L, 820854L, 822085L, 
                                                                                                                         826961L, 827132L, 827564L, 827556L, 827338L, 827553L, 827348L, 
                                                                                                                         827360L, 827607L, 827501L, 827500L, 827596L, 827540L, 827561L, 
                                                                                                                         827568L, 827602L, 827569L, 827562L, 827617L, 827612L, 827619L, 
                                                                                                                         827648L, 827650L, 827635L, 827768L), V4 = c("AAACTGCTCTTCCACG-1", 
                                                                                                                                                                     "AAAGATGCACGCTGTG-1", "AAACTGCCAGTAGGCA-1", "AAACTGCAGGCAAGCT-1", 
                                                                                                                                                                     "AAAGATGAGACTAGGC-1", "AAACGAATCTCGTAGA-1", "AAACTGCAGGCAAGCT-1", 
                                                                                                                                                                     "AAAGATGGTGCAAGAC-1", "AAACGAACACCCTTTG-1", "AAACTCGAGGTCTTTG-1", 
                                                                                                                                                                     "AAACGAACACCCTTTG-1", "AAACGAATCTCGTAGA-1", "AAACTGCCAATGGTCT-1", 
                                                                                                                                                                     "AAACTCGAGGTCTTTG-1", "AAAGGATGTAGGTGAC-1", "AAACTCGGTCCGAGCT-1", 
                                                                                                                                                                     "AAAGATGCACGCTGTG-1", "AAACTCGGTGTGAGGT-1", "AAACTCGGTGTGAGGT-1", 
                                                                                                                                                                     "AAACTGCCAAAGCTGG-1", "AAACTGCCATCGCCTT-1", "AAACTCGCAAAGAAGG-1", 
                                                                                                                                                                     "AAACTGCCATCGCCTT-1", "AAACTCGAGCTATCCA-1", "AAACTGCCAAACGTTC-1", 
                                                                                                                                                                     "AAACTGCCAGTGCTCG-1", "AAACTCGTCCGTTAGA-1", "AAAGGATGTAGGTGAC-1", 
                                                                                                                                                                     "AAAGATGTCCGCTCTA-1", "AAACTCGTCTTAGTGG-1", "AAAGATGAGTACCACT-1", 
                                                                                                                                                                     "AAACGAATCGAGTTAC-1", "AAACTGCGTCCTTCAC-1", "AAAGATGGTGCTTGAT-1", 
                                                                                                                                                                     "AAACTGCTCGAGGTAG-1", "AAACTCGGTCCGAGCT-1", "AAAGGATGTAGGTGAC-1", 
                                                                                                                                                                     "AAAGATGGTGCAAGAC-1", "AAACTGCCAATGGTCT-1", "AAACTCGAGCTATCCA-1", 
                                                                                                                                                                     "AAACTCGTCCGTTAGA-1", "AAACGAATCTCGTAGA-1", "AAACTGCCAGTAGGCA-1", 
                                                                                                                                                                     "AAACTGCCAAACGTTC-1", "AAACTGCCAGTGCTCG-1", "AAACTGCCATCGCCTT-1", 
                                                                                                                                                                     "AAACTCGAGCTATCCA-1", "AAACTGCCAGTGCTCG-1", "AAACTGCGTCCTTCAC-1", 
                                                                                                                                                                     "AAACTGCCAAAGCTGG-1", "AAACTGCAGGCAAGCT-1", "AAACTGCTCACCCTTG-1", 
                                                                                                                                                                     "AAACTGCTCTTCCACG-1", "AAACTCGTCTTAGTGG-1", "AAACGAAAGTCGACCC-1", 
                                                                                                                                                                     "AAAGATGAGTACCACT-1", "AAACTCGAGGTCTTTG-1", "AAACTGCTCACCCTTG-1", 
                                                                                                                                                                     "AAACTGCCAAACGTTC-1", "AAACTGCTCTTCCACG-1", "AAACGAACACATAAAG-1", 
                                                                                                                                                                     "AAACTCGAGCTATCCA-1", "AAACTGCCAAAGCTGG-1", "AAAGATGAGTACCACT-1", 
                                                                                                                                                                     "AAACTCGGTCCGAGCT-1", "AAACTCGAGCTATCCA-1", "AAACTCGGTCCGAGCT-1", 
                                                                                                                                                                     "AAAGATGTCCGCTCTA-1", "AAACTCGTCCGTTAGA-1", "AAACTGCCAGTGCTCG-1", 
                                                                                                                                                                     "AAAGATGAGACTAGGC-1", "AAAGATGAGACTAGGC-1", "AAAGGATGTAGGTGAC-1", 
                                                                                                                                                                     "AAACGAAGTGCTGGCT-1", "AAAGATGCACGCTGTG-1", "AAAGGATAGACGCCCT-1", 
                                                                                                                                                                     "AAAGGATGTAGGTGAC-1", "AAACTCGAGCTATCCA-1", "AAAGATGAGTACCACT-1", 
                                                                                                                                                                     "AAACGAAAGTCGACCC-1", "AAACTGCCAATGGTCT-1", "AAAGATGGTGCAAGAC-1", 
                                                                                                                                                                     "AAACTCGAGGTCTTTG-1", "AAACTGCTCGAGGTAG-1", "AAACTGCCAAACGTTC-1", 
                                                                                                                                                                     "AAACGAACACATAAAG-1", "AAACTCGTCGGTGATT-1", "AAAGATGAGACTAGGC-1", 
                                                                                                                                                                     "AAAGATGAGACTAGGC-1", "AAACTCGAGCTATCCA-1", "AAAGGATAGACGCCCT-1", 
                                                                                                                                                                     "AAACGAAAGTCGACCC-1", "AAAGGATAGACGCCCT-1", "AAACTGCCATCGCCTT-1", 
                                                                                                                                                                     "AAACTGCCAGTAGGCA-1", "AAACGAACACATAAAG-1", "AAACTGCCAGTGCTCG-1", 
                                                                                                                                                                     "AAACTCGAGCTATCCA-1", "AAACTCGCAAAGAAGG-1", "AAACTCGGTCCGAGCT-1", 
                                                                                                                                                                     "AAACTGCCAGTAGGCA-1", "AAACGAATCCAGAATC-1", "AAAGATGCAAAGAAGG-1", 
                                                                                                                                                                     "AAACTGCCATCGCCTT-1", "AAACGAACACATATCG-1"), V5 = c(2L, 1L, 1L, 
                                                                                                                                                                                                                         1L, 3L, 1L, 2L, 9L, 1L, 1L, 1L, 1L, 2L, 1L, 6L, 1L, 1L, 3L, 3L, 
                                                                                                                                                                                                                         1L, 1L, 4L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 4L, 4L, 
                                                                                                                                                                                                                         1L, 1L, 1L, 3L, 3L, 3L, 5L, 1L, 1L, 1L, 1L, 3L, 1L, 2L, 2L, 1L, 
                                                                                                                                                                                                                         3L, 2L, 1L, 3L, 1L, 2L, 5L, 1L, 1L, 3L, 1L, 1L, 1L, 4L, 4L, 1L, 
                                                                                                                                                                                                                         5L, 1L, 2L, 3L, 1L, 1L, 1L, 2L, 5L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 
                                                                                                                                                                                                                         4L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 6L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 
                                                                                                                                                                                                                         3L, 2L, 4L, 2L, 1L, 1L)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                       -105L))
  )
  file.remove(of1, of2)

  Fragments(atac_small) <- NULL
  Fragments(atac_small) <- frags_headered
  SplitFragments(
    object = atac_small,
    assay = "peaks",
    group.by = "cluster",
    outdir = tempdir()
  )
  of1 <- paste0(tempdir(), .Platform$file.sep, "1.bed")
  of2 <- paste0(tempdir(), .Platform$file.sep, "2.bed")

  bed1_h <- read.table(file = of1, sep = "\t", stringsAsFactors = FALSE)
  bed2_h <- read.table(file = of2, sep = "\t", stringsAsFactors = FALSE)

  expect_equal(
    object = bed1_h,
    expected = bed1
  )

  expect_equal(
    object = bed2_h,
    expected = bed2
  )

})

suppressWarnings(RNGversion(vstr = "3.5.3"))

count_vec <- c(4, 3, 2, 3, 2, 2, 2, 8, 1, 1, 1, 4, 3, 1, 5, 2, 2, 2, 1, 1, 2,
               1, 1, 3, 1, 1, 2, 1, 3, 2, 1, 2, 1, 2, 2, 1, 3, 1, 2, 1, 1, 2,
               1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1)

mononucleosome <- c(1, 1, 0, 2, 0, 0, 1, 5, 0, 0, 1, 4, 1, 1, 4, 1, 1, 1, 0, 1, 1,
                    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1,
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)

nucleosome_free <- c(3, 1, 2, 0, 2, 2, 0, 2, 0, 0, 0, 0, 2, 0, 1, 1, 1, 1, 1, 0, 1,
                     0, 1, 2, 1, 1, 2, 1, 3, 2, 1, 2, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1,
                     0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1)

read_vec <- c(5, 4, 2, 5, 3, 2, 3, 13, 1, 2, 6, 8, 5, 1, 8, 4, 6, 3, 1, 1, 2,
              4, 1, 4, 1, 8, 12, 1, 4, 4, 2, 19, 4, 3, 3, 1, 4, 3, 4, 1, 1,
              7, 1, 1, 1, 1, 1, 3, 2, 1, 3, 5, 2, 1)

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
    expected = c(0, 0, 3, 8, 0, 0)
  )
  expect_equal(
    object = head(x = counts$mononucleosomal),
    expected = c(0, 0, 0, 5, 0, 0)
  )

  expect_equal(
    object = head(x = counts$nucleosome_free),
    expected = c(0, 0, 3, 2, 0, 0)
  )

  expect_equal(
    object = head(x = counts$reads_count),
    expected = c(0, 0, 4, 13, 0, 0)
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

  expected <- structure(list(V1 = c("chr1", "chr1", "chr1",
                                    "chr1", "chr1", "chr1",
                                    "chr1", "chr1", "chr1",
                                    "chr1", "chr1"),
                             V2 = c(712868L, 713783L, 713944L,
                                    714140L, 714144L, 714263L,
                                    757378L, 762811L, 762874L,
                                    762951L, 773225L),
                             V3 = c(713146L, 714045L, 713997L,
                                    714174L, 714209L, 714732L,
                                    757538L, 762953L, 762953L,
                                    763224L, 773453L),
                             V4 = c("AAACGAAAGGCTTCGC-1", "AAACGAAAGGCTTCGC-1",
                                    "AAACGAAAGGCTTCGC-1", "AAACGAAAGCGAGCTA-1",
                                    "AAACGAAAGGCTTCGC-1", "AAACGAAAGGCTTCGC-1",
                                    "AAACGAAAGGCTTCGC-1", "AAACGAAAGCGAGCTA-1",
                                    "AAACGAAAGCGAGCTA-1", "AAACGAAAGGCTTCGC-1",
                                    "AAACGAAAGGCTTCGC-1"),
                             V5 = c(1L, 2L, 1L, 1L, 3L, 1L,
                                    2L, 1L, 2L, 1L, 2L)),
                        class = "data.frame", row.names = c(NA, -11L))
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
    group.by = "seurat_clusters",
    outdir = tempdir()
  )
  of1 <- paste0(tempdir(), .Platform$file.sep, "0.bed")
  of2 <- paste0(tempdir(), .Platform$file.sep, "1.bed")

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
                                     "chr1", "chr1", "chr1", "chr1"),
                              V2 = c(10245L, 55313L, 56455L,
                                     237741L, 241022L, 713591L, 713734L, 713751L, 713832L, 713934L,
                                     714017L, 714017L, 714017L, 714028L, 714043L, 714043L, 714131L,
                                     714140L, 714140L, 714182L, 714183L, 714241L, 714289L, 714454L,
                                     714464L, 714597L, 752259L, 752682L, 752684L, 752699L, 752702L,
                                     755549L, 762180L, 762531L, 762629L, 762652L, 762668L, 762718L,
                                     762738L, 762811L, 762811L, 762845L, 762874L, 762927L, 770711L,
                                     777349L, 778251L, 779103L, 779777L),
                              V3 = c(10302L, 55699L, 56658L,
                                     237772L, 241062L, 714040L, 714120L, 714035L, 714003L, 714035L,
                                     714043L, 714092L, 714213L, 714145L, 714148L, 714178L, 714186L,
                                     714174L, 714241L, 714284L, 714221L, 714345L, 714339L, 714485L,
                                     714777L, 715147L, 752335L, 752837L, 752825L, 752877L, 752831L,
                                     755587L, 762568L, 762907L, 762881L, 762938L, 762726L, 762904L,
                                     762881L, 762864L, 762953L, 762881L, 762953L, 762972L, 771123L,
                                     777590L, 778277L, 779346L, 780007L),
                              V4 = c("AAAGATGAGGCTAAAT-1",
                                     "AAACTCGTCTGGCACG-1", "AAACTCGTCTGGCACG-1", "AAAGGATTCCTTACGC-1",
                                     "AAAGATGAGAAGGGCG-1", "AAACTCGCAGCGTCGT-1", "AAACGAATCCTTACGC-1",
                                     "AAAGGATCATGGAGGT-1", "AAACTGCAGAATCAAC-1", "AAAGATGCAAGTCTGT-1",
                                     "AAACTGCCATGTATCG-1", "AAAGGATTCCTTACGC-1", "AAACTGCCAAGCCAGA-1",
                                     "AAAGGATTCTATGAGC-1", "AAACTGCAGCTCCATA-1", "AAAGATGAGTCCAGAG-1",
                                     "AAACTGCTCATTCATC-1", "AAACGAAAGCGAGCTA-1", "AAACTCGGTGGATTCT-1",
                                     "AAACTCGTCTGGCACG-1", "AAAGATGTCTAGCAGT-1", "AAACTGCGTACTAGAA-1",
                                     "AAAGATGAGGCTAAAT-1", "AAAGATGAGGCTAAAT-1", "AAAGATGTCCTGAAAC-1",
                                     "AAAGGATCAAGGGAGG-1", "AAACTCGGTGGATTCT-1", "AAAGGATTCTATGAGC-1",
                                     "AAAGGATGTCCCTTTG-1", "AAAGATGCATGACTGT-1", "AAAGATGAGAAGGGCG-1",
                                     "AAAGATGCAGCGTCGT-1", "AAACGAAGTTGTATCG-1", "AAACTCGCATGTGGGA-1",
                                     "AAACTGCGTGCTGAAG-1", "AAAGATGCAAGTCTGT-1", "AAAGATGCATGACTGT-1",
                                     "AAACTGCTCTGAGTAC-1", "AAACTCGAGTCACGCC-1", "AAAGGATTCCAACCTC-1",
                                     "AAACGAAAGCGAGCTA-1", "AAAGGATTCTATGAGC-1", "AAACGAAAGCGAGCTA-1",
                                     "AAAGATGAGTCCAGAG-1", "AAACTCGAGTCACGCC-1", "AAAGGATCAAGGGAGG-1",
                                     "AAAGATGTCCTGAAAC-1", "AAACTCGCATGTGGGA-1", "AAAGATGAGGCTAAAT-1"
                                                                                                                         ),
                              V5 = c(1L, 2L, 1L, 1L, 1L, 1L, 2L, 6L, 1L, 3L, 1L, 2L, 1L,
                                     1L, 8L, 9L, 1L, 1L, 3L, 1L, 2L, 4L, 2L, 1L, 2L, 1L, 1L, 1L, 1L,
                                     3L, 1L, 1L, 1L, 1L, 2L, 1L, 4L, 3L, 4L, 2L, 1L, 2L, 2L, 3L, 1L,
                                     2L, 1L, 2L, 1L)), class = "data.frame", row.names = c(NA, -49L
                                     ))
  )

  expect_equal(
    object = bed2,
    expected = structure(list(V1 = c("chr1", "chr1", "chr1", "chr1", "chr1",
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1",
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1",
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1",
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1",
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1",
                                     "chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
                              V2 = c(60687L,
                                     235723L, 526949L, 565296L, 712868L, 713684L, 713752L, 713780L,
                                     713783L, 713783L, 713856L, 713934L, 713944L, 713989L, 713990L,
                                     714022L, 714027L, 714028L, 714040L, 714114L, 714125L, 714144L,
                                     714190L, 714263L, 714464L, 722942L, 738117L, 752435L, 752651L,
                                     752660L, 752688L, 753227L, 754851L, 757378L, 762083L, 762450L,
                                     762629L, 762657L, 762778L, 762778L, 762801L, 762859L, 762892L,
                                     762895L, 762923L, 762951L, 764656L, 773187L, 773225L, 778213L,
                                     779523L),
                              V3 = c(60726L, 235936L, 527161L, 565351L, 713146L,
                                    714108L, 713990L, 713992L, 714025L, 714045L, 714048L, 714176L,
                                    713997L, 714247L, 714322L, 714182L, 714214L, 714144L, 714108L,
                                    714207L, 714193L, 714209L, 714335L, 714732L, 714573L, 723275L,
                                    738167L, 752819L, 752827L, 752875L, 752772L, 753684L, 755017L,
                                    757538L, 762426L, 762501L, 762982L, 762859L, 762836L, 763004L,
                                    762830L, 762971L, 762994L, 762938L, 763111L, 763224L, 764836L,
                                    773410L, 773453L, 778402L, 779769L),
                              V4 = c("AAACTGCAGTCTGTGT-1",
                                    "AAACTGCTCCTATCCG-1", "AAACGAAGTGCCCGAT-1", "AAACTGCAGTCTGTGT-1",
                                    "AAACGAAAGGCTTCGC-1", "AAACGAAGTGCCCGAT-1", "AAACTGCTCCTATCCG-1",
                                    "AAAGATGCACGTTACA-1", "AAAGATGTCCACACCT-1", "AAACGAAAGGCTTCGC-1",
                                    "AAACTCGGTTTGATCG-1", "AAACTCGTCAGGTCTA-1", "AAACGAAAGGCTTCGC-1",
                                    "AAACTGCGTGCATTCA-1", "AAACTGCTCCTATCCG-1", "AAACTCGAGTACTCTG-1",
                                    "AAACTGCTCGTTCAGA-1", "AAAGGATGTCGTAATC-1", "AAACGAAGTCAGGCTC-1",
                                    "AAAGATGTCCACACCT-1", "AAAGATGTCCACACCT-1", "AAACGAAAGGCTTCGC-1",
                                    "AAAGGATCAGATGGCA-1", "AAACGAAAGGCTTCGC-1", "AAACTCGAGTACTCTG-1",
                                    "AAAGATGCAGCACATT-1", "AAACTCGGTTTGATCG-1", "AAACTCGCAAAGAGAG-1",
                                    "AAAGATGAGTCGACCC-1", "AAACTCGCATCACAGT-1", "AAACTGCCAACAAACA-1",
                                    "AAACTCGTCCACTAGA-1", "AAAGATGCACGTTACA-1", "AAACGAAAGGCTTCGC-1",
                                    "AAAGGATGTCTAAGAA-1", "AAAGGATAGAACCATA-1", "AAACTCGCATCACAGT-1",
                                    "AAAGGATGTCTTAGCA-1", "AAACTCGCAAAGAGAG-1", "AAAGATGCACGTTACA-1",
                                    "AAACTCGTCAGGTCTA-1", "AAAGGATCAGATGGCA-1", "AAACTCGAGCGCATTT-1",
                                    "AAACTGCGTGCATTCA-1", "AAACTCGCAAAGAGAG-1", "AAACGAAAGGCTTCGC-1",
                                    "AAACTCGGTTTGATCG-1", "AAACTCGGTTTGATCG-1", "AAACGAAAGGCTTCGC-1",
                                    "AAAGATGCACGTTACA-1", "AAACTCGGTTTGATCG-1"),
                              V5 = c(1L, 1L, 2L,
                                    1L, 1L, 1L, 3L, 3L, 2L, 2L, 1L, 3L, 1L, 2L, 1L, 1L, 4L, 1L, 1L,
                                    1L, 2L, 3L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 1L, 1L, 1L, 2L, 1L,
                                    1L, 1L, 1L, 2L, 2L, 3L, 18L, 1L, 1L, 1L, 1L, 1L, 4L, 2L, 2L,
                                                                                                                                      1L)),
                         class = "data.frame", row.names = c(NA, -51L))
  )
  file.remove(of1, of2)

  Fragments(atac_small) <- NULL
  Fragments(atac_small) <- frags_headered
  SplitFragments(
    object = atac_small,
    assay = "peaks",
    group.by = "seurat_clusters",
    outdir = tempdir()
  )
  of1 <- paste0(tempdir(), .Platform$file.sep, "0.bed")
  of2 <- paste0(tempdir(), .Platform$file.sep, "1.bed")

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

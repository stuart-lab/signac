library(GenomicRanges)
library(SeuratObject)

# data set up
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
cells <- colnames(x = atac_small)
atac_small2 <- atac_small

# subset cells & write to file
bc <- colnames(x = atac_small)[1:50]
cells_path <- file.path(tempdir(), "test_barcodes.txt")
writeLines(text = bc, con = cells_path)

frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE)
Fragments(atac_small) <- frags

# rename cells
atac_small <- RenameCells(atac_small, "test")

# make groups
seurat_clusters <- c(rep("cell1",20), rep("cell2",20), rep("cell3",20), rep("cell4",20), rep("cell5",20))
atac_small <- AddMetaData(atac_small, metadata = seurat_clusters, col.name = "seurat_clusters")
Idents(atac_small) <- atac_small$seurat_clusters

# make 2nd object
fpath_headered <- system.file("extdata", "fragments_header.tsv.gz", package="Signac")
frags_headered <- CreateFragmentObject(path = fpath_headered, cells = cells, verbose = FALSE)
Fragments(atac_small2) <- frags_headered
atac_small2 <- RenameCells(atac_small2, "test")

seurat_clusters <- c(rep("cell1",20), rep("cell2",20), rep("cell3",20), rep("cell4",20), rep("cell5",20))
atac_small2 <- AddMetaData(atac_small2, metadata = seurat_clusters, col.name = "seurat_clusters")
Idents(atac_small2) <- atac_small2$seurat_clusters

# combine obj 
combined <- merge(
  x = atac_small,
  y = atac_small2,
  add.cell.ids = c("obj1", "obj2")
)


# Default parameters
test_that("CallPeaks with default parameters works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    output_default_params <- GRanges(
        seqnames = Rle(rep("chr1", 9)),
        ranges = IRanges(
            start = c(1, 605581, 778386, 811334, 817148, 818778, 820528, 821060, 827279),
            end = c(581234, 605805, 779097, 811504, 817502, 819129, 820684, 821302, 827650)
        ),
        strand = Rle(rep("*", 9)),
        name = paste0("macs3_peak_", 1:9),
        score = c(45, 136, 571, 34, 348, 19, 21, 18, 410),
        fold_change = c(3.87432, 9.11316, 17.90100, 3.75631, 18.99140, 2.11016, 2.81354, 2.04786, 19.64450),
        neg_log10pvalue_summit = c(7.34726, 16.71640, 63.09270, 6.15439, 38.33860, 2.03903, 3.02567, 1.92558, 44.78700),
        neg_log10qvalue_summit = c(4.54594, 13.68160, 57.17480, 3.44248, 34.82890, 1.97269, 2.19277, 1.87557, 41.04510),
        relative_summit_position = c(190928, 112, 280, 85, 199, 176, 155, 121, 256)
    )
    
    # Fragment path
    expect_equal(
        object = CallPeaks(object = fpath, verbose = FALSE),
        expected = output_default_params,
        tolerance = 0.5
    )
    # Fragment2 class
    expect_equal(
        object = CallPeaks(object = frags, verbose = FALSE),
        expected = output_default_params,
        tolerance = 0.5
    )
    # ChromatinAssay5 class
    expect_equal(
        object = CallPeaks(object = atac_small[['peaks']], verbose = FALSE),
        expected = output_default_params,
        tolerance = 0.5
    )
    # Seurat class
    expect_equal(
        object = CallPeaks(object = atac_small, name="macs3", verbose = FALSE),
        expected = output_default_params,
        tolerance = 0.5
    )
})

test_that("CallPeaks with cell barcodes works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    output_cells <- GRanges(
        seqnames = Rle(rep("chr1", 8)),
        ranges = IRanges(
            start = c(1, 778358, 817266, 818778, 820528, 821060, 821920, 827276),
            end   = c(629052, 778948, 817418, 819129, 820684, 821302, 822085, 827648)
        ),
        strand = Rle(rep("*", 8)),
        name = paste0("macs31_peak_", 1:8),
        score = c(42, 374, 21, 21, 21, 21, 20, 293),
        fold_change = c(3.87558, 16.62590, 3.51151, 2.68406, 2.68406, 2.68406, 1.77240, 16.63400),
        neg_log10pvalue_summit = c(7.36518, 43.33380, 4.85497, 3.60393, 3.60393, 3.60393, 2.12080, 32.86260),
        neg_log10qvalue_summit = c(4.26018, 37.41590, 2.16093, 2.16093, 2.16093, 2.16093, 2.07912, 29.30640),
        relative_summit_position = c(586214, 303, 48, 176, 78, 121, 101, 259)
    )
    
    # CallPeaks.Default with barcode path
    expect_equal(
        object = CallPeaks(
          object = fpath,
          barcodes = cells_path,
          cleanup = TRUE,
          name="macs31",
          verbose = FALSE
        ),
        expected = output_cells,
        tolerance = 0.5
    )
    # CallPeaks.Fragment2 with cells vector
    expect_equal(
        object = CallPeaks(
          object = frags,
          cells = names(x = frags@cells)[1:50],
          cleanup = TRUE,
          name = "macs31",
          verbose = FALSE
        ),
        expected = output_cells,
        tolerance = 0.5
    )
    # CallPeaks.ChromatinAssay5 with cells vector
    expect_equal(
        object = CallPeaks(
          object = atac_small[['peaks']],
          cells = Cells(x = atac_small[['peaks']])[1:50],
          cleanup = TRUE,
          verbose = FALSE
        ),
        expected = output_cells,
        tolerance = 0.5
    )
    # CallPeaks.Seurat with cells vector
    expect_equal(
        object = CallPeaks(
          object = atac_small,
          cells = Cells(x = atac_small)[1:50],
          cleanup = TRUE,
          name = "macs3",
          verbose = FALSE
        ),
        expected = output_cells,
        tolerance = 0.5
    )
})

test_that("CallPeaks group.by works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    expect_equal(
        object = CallPeaks(object = atac_small, group.by = "seurat_clusters"),
        expected = GRanges(
            seqnames = Rle(rep("chr1", 9)),
            ranges = IRanges(
                start = c(1, 777690, 778275, 779231, 783782, 822569, 823345, 826723, 827081),
                end = c(773637, 777881, 779099, 779384, 822477, 822796, 823502, 826961, 827801)
            ),
            strand = Rle(rep("*", 9)),
            peak_called_in = c("cell5,cell4,cell3,cell2,cell1", "cell5", "cell1,cell3,cell5,cell2,cell4", "cell2", "cell4,cell1,cell2,cell3,cell5", "cell2", "cell5", "cell3", "cell2,cell4,cell1,cell5,cell3")
        ),
        tolerance = 0.5
    )
    expect_length(
        object = CallPeaks(object = atac_small, group.by = "seurat_clusters", combine.peaks = FALSE),
        n = length(unique(atac_small$seurat_clusters))
    )
})

test_that("CallPeaks group.by & idents works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    expect_equal(
        object = CallPeaks(object = atac_small, group.by = "seurat_clusters", idents = "cell5"),
        expected = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = c(1, 605581, 777690, 778386, 817265, 823345, 827404),
                end = c(600580, 605805, 777881, 779099, 817508, 823502, 827607)
            ),
            strand = "*",
            peak_called_in = rep("cell5", 7)
        ),
        tolerance = 0.5
    )
    expect_equal(
        object = CallPeaks(object = atac_small, group.by = "seurat_clusters", idents = c("cell1","cell5")),
        expected = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = c(1, 777690, 778275, 783845, 817265, 818175, 823345, 827294),
                end = c(773637, 777881, 779099, 814129, 817508, 822477, 823502, 827768)
            ),
            strand = "*",
            peak_called_in = c("cell5,cell1", "cell5", "cell1,cell5", "cell1", "cell5", "cell1", "cell5", "cell1,cell5")
        ),
        tolerance = 0.5
    )
    expect_length(
        object = CallPeaks(object = atac_small, group.by = "seurat_clusters", idents = c("cell1","cell5"), combine.peaks = FALSE),
        n = length(c("cell1","cell5"))
    )
    expect_error(
        object = CallPeaks(object = atac_small, group.by = NULL, idents = "cell1")
    )
})

test_that("CallPeaks group.by & idents & cells works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    expect_equal(
        object = CallPeaks(object = atac_small, group.by = "seurat_clusters", idents = c("cell1","cell3"), cells = Cells(atac_small)[1:50]),
        expected = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = c(1, 778275, 783756, 826723, 827294),
                end   = c(773637, 778965, 822501, 826961, 827768)
            ),
            strand = "*",
            peak_called_in = c("cell3,cell1", "cell1,cell3", "cell3,cell1", "cell3", "cell1,cell3")
        ),
        tolerance = 0.5
    )
    expect_error(
        object = CallPeaks(object = atac_small, group.by = NULL, idents = c("cell1","cell3"), cells = Cells(atac_small)[1:50]),
    )
})

test_that("CallPeaks broad option works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    # broad option
    expect_equal(
        object = CallPeaks(object = atac_small, broad = TRUE, name = 'macs3'),
        expected = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = c(1, 778386, 783950, 817070, 818778, 820528, 821920, 823345, 826723),
                end   = c(629018, 779097, 812265, 817508, 819129, 821302, 822085, 823502, 827768)
            ),
            strand = "*",
            name = paste0("macs3_peak_", 1:9),
            score = c(21, 241, 20, 78, 19, 16, 10, 11, 119),
            fold_change = c(1.00103, 9.33481, 1.00555, 5.71701, 2.11016, 1.87752, 1.34606, 1.56325, 7.61677),
            neg_log10pvalue_summit = c(4.46216, 26.96830, 4.07226, 9.23500, 2.03903, 1.73695, 1.06585, 1.18023, 14.00300),
            neg_log10qvalue_summit = c(2.12409, 24.13070, 2.00120, 7.80066, 1.97269, 1.68636, 1.04291, 1.15680, 11.94190)
        ),
        tolerance = 0.5
    )
})

test_that("CallPeaks multiple fragments works", {
    # Check if macs3 is available in the system PATH
    skip_if(Sys.which("macs3") == "", message = "macs3 command line tool not found, skipping test-callpeaks")
    
    # CallPeaks.Seurat with multiple fragments
    expect_equal(
        object = CallPeaks(object = combined, name = "macs3", verbose = FALSE),
        expected = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = c(1, 190825, 195825, 605581, 778386, 811334, 817148, 818778, 820528, 821060, 823345, 827276),
                end   = c(185824, 191033, 581197, 605805, 779097, 811504, 817502, 819129, 820684, 821302, 823502, 827768)
            ),
            strand = "*",
            name = paste0("macs3_peak_", 1:12),
            score = c(62, 92, 33, 276, 1177, 67, 704, 19, 19, 19, 13, 831),
            fold_change = c(4.79271, 6.57353, 2.93650, 15.90450, 20.69980, 6.19605, 28.75140, 2.71239, 3.79735, 2.59082, 1.76171, 26.64010),
            neg_log10pvalue_summit = c(8.91486, 12.04250, 5.78030, 30.67050, 123.65500, 9.95989, 73.91090, 2.75106, 4.53877, 2.56963, 1.40638, 86.91940),
            neg_log10qvalue_summit = c(6.25629, 9.24149, 3.15374, 27.63560, 117.73700, 7.24798, 70.40790, 1.99291, 1.99412, 1.99291, 1.35706, 83.17760),
            relative_summit_position = c(10192, 104, 72182, 112, 280, 85, 199, 176, 155, 121, 79, 259)
        ),
        tolerance = 0.5
    )
    # output types & length
    groupby_combined <- CallPeaks(object = combined, combine.peaks=TRUE, group.by="seurat_clusters")
    groupby_separate <- CallPeaks(object = combined, combine.peaks=FALSE, group.by="seurat_clusters")
    
    expect_s4_class(groupby_combined, 'GRanges')
    expect_type(groupby_separate, 'list')
    expect_length(
        object = groupby_separate,
        n = 5
    )
    
    expect_length(
        object = CallPeaks(object = combined, combine.peaks=FALSE, group.by="seurat_clusters", idents = c("cell1", "cell2")),
        n = 2
    )
})
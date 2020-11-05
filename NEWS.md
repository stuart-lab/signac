# Signac 1.1.0

New functionality:

* Added `CallPeaks()` function to call peaks using MACS2. Peaks can be called
for different groups of cells separately by setting the `group.by` parameter
* Added `LinkPeaks()` function to link peaks to correlated genes.
* Added `AddMotifs()` function to add motif information to a Seurat object or ChromatinAssay.
* Added `AggregateTiles()` function to combine adjacent genome tiles
* Added `ranges` parameter to `CoveragePlot()` to plot addition sets of genomic ranges
* Added `show.bulk` parameter to `CoveragePlot()` to plot accessibility of all cells combined
* Added ability to remove `Fragment` objects and modify the file path for existing
fragment objects ([#206](https://github.com/timoast/signac/issues/206))

Bug fixes: 

* Fixed bugs in `AlleleFreq()` ([#196](https://github.com/timoast/signac/issues/196)
and [#260](https://github.com/timoast/signac/issues/260))
* Fixed bug in `FeatureMatrix()` ([#205](https://github.com/timoast/signac/issues/205), [#291](https://github.com/timoast/signac/issues/291))
* Fixed bug in `CreateChromatinAssay()` when setting `min.features` argument ([#194](https://github.com/timoast/signac/issues/194))
* Fixed bug in `CreateChromatinAssay()` when setting `min.cells` argument ([#292](https://github.com/timoast/signac/issues/292))
* Fixed bug in `TSSEnrichment()` when cell information not set for fragment files ([#203](https://github.com/timoast/signac/issues/203))
* Fixed bug in `TSSEnrichment()` when no fragments present in TSS region ([#244](https://github.com/timoast/signac/issues/244))
* Removed `qvalue` calculation from `FindMotifs()` ([#223](https://github.com/timoast/signac/issues/223))
* Fixed bug in `SetAssayData()` when setting the `scale.data` slot

Other changes:

* Improved feature matching in `MatchRegionStats()` function when matching distribution of multiple features (eg, GC content and overall accessibility)
* Changed parameter names in `MatchRegionStats()`

# Signac 1.0.0

This release includes major updates to the Signac package, including new
functionality, performance improvements, and new data structures.

The entire package has been updated to use the new `ChromatinAssay` class for the
storage of single-cell chromatin data. This is an extension of the standard 
Seurat `Assay` that adds additional slots needed for the analysis of chromatin 
data, including genomic ranges, genome information, fragment file information,
motifs, gene annotations, and genomic links.

In addition, we have defined a new `Fragment` class to store information 
relating to a fragment file. This makes use of the fragment files within Signac
more robust, as checks are now performed to verify that the expected cells are
present in the fragment file, and that the fragment file or index are not
modified on disk.

Key new functionality:

* **Store multiple fragment files**: you can now store as many fragment
files as needed in a single object, and all functions that use the fragment file
will pull data from each of the files. Cell barcodes in the fragment files do
_not_ need to match the cell barcodes in the object.
* **Use remote fragment files**: you can now use all the same functionality with
fragment files hosted on remote servers accessible through `http` or `ftp`.
* **Transcription factor footprinting**: New `Footprint()` and `PlotFootprint()`
functions for TF footprinting analysis.
* **Bioconductor methods**: call `granges()`, `findOverlaps()`, `seqinfo()`, and
other Bioconductor generic functions directly on the `ChromatinAssay` or
`Seurat` object.
* **New multi-modal visualization methods**: Jointly visualize RNA expression
and chromatin accessibility using the `CoveragePlot()` function.
* **New interactive visualizations**: Interactively browse the genome using the
`CoverageBrowser()` function.
* **Mitochondrial lineage tracing**: New functions to identify informative
mitochondrial alleles, find clonotypes, and predict cell lineage relationships
using mitochondrial mutations.

Other changes:

* Updates to `NucleosomeSignal()`: we have greatly improved the scalability of
`NucleosomeSignal()`, and fixed a bug present in previous versions. The score
computed by `NucleosomeSignal()` in 1.0.0 will be different to that computed by
previous versions of Signac.
* New `CountFragments()` function: a fast, memory-efficient function implemented
in C++ that counts the total number of fragments for each cell barcode present
in a fragment file.
* New `fast` option in the `TSSEnrichment()` function. Setting this to `TRUE`
will compute the TSS enrichment score per cell without storing the entire
cell by TSS position matrix. This can significantly reduce memory requirements
for large datasets, but does not allow subsequent plotting of the TSS signal
for different groups of cells.
* New `TilePlot()` function and `tile` parameter for `CoveragePlot()` to plot
Tn5 integration events in a genomic region for individual cells.
* Performance improvements for `FeatureMatrix()`, `CoveragePlot()`, and
`TSSEnrichment()`
* Added the manually curated hg38 genomic blacklist regions curated by Anshul
Kundaje and Anna Shcherbina. These are available as the `blacklist_hg38_unified`
object.
* Updated the `FRiP()` function to use total fragment counts per cell stored
in object metadata.

# Signac 0.2.5

* New `DepthCor` function to compute the correlation between sequencing depth and
reduced dimension components. 
* Performance improvements for `RunTFIDF`. 
* Removed option to use EnsDb object in `ClosestFeatures` and `CoveragePlot`. Use GRanges instead. 
* Removed `ucsc` parameter from `CoveragePlot`. 
* Fixed bug in FeatureMatrix that would cause fragments to be counted multiple
times if `nchunk` was greater than the number of features used. 
* Fixed bug in `CoveragePlot` that would prevent plotting multiple regions when
using `GRanges`. 
* Fixed bug in `CoveragePlot` that would prevent plotting when a different 
assay was active. 
* Removed dependencies: GenomicFeatures
* Moved dependencies to suggests: Biostrings, BSgenome
* Removed from suggests: BSgenome.Hsapiens.UCSC.hg19, EnsDb.Hsapiens.v75, JASPAR2018

# Signac 0.2.4

* First CRAN release.
* New `SubsetMatrix` function to subset a matrix based on number of non-zero elements 
in the rows or columns.
* Removed `seed.use` parameter from `RunSVD`.

# Signac 0.2.3

* New `UnifyPeaks` function to create a merged set of peaks from multiple samples.

# Signac 0.2.2

* Bug fix for `RunSVD`: previously, scaling was applied to each cell rather than each component.
Now, mean centering and SD scaling are applied to the cell embeddings within a component.
* Added `scale.embeddings` option to `RunSVD` to control whether embeddings are scaled
and centered.
* Added `irlba.work` parameter to `RunSVD`.
* Update to allow comment characters in fragment file cell names

# Signac 0.2.1

* Removed `SingleCoveragePlot` from exported functions
* Added executable examples for all functions
* Store raw SVD output in DimReduc misc slot in `RunSVD`
* Fixed strand orientation for gene plot in `CoveragePlot`
* Fix missing x-axis when plotting peaks but not genes in `CoveragePlot`

# Signac 0.2.0

* Removed dependency on TFBSTools, motifmatchr, AnnotationDbi, ggbio, AnnotationFilter
* Renamed `PeriodPlot` to `FragmentHistogram`
* Removed motif dimension reduction functions
* Removed motif clustering functions
* Removed `neighbors` and `reductions` slots from `motif` class
* Added `motif.names` slot to `motif` class
* Added ability to plot peak ranges in `CoveragePlot`
* Added ability to plot gene annotations from `GRanges` object
* Changed gene plot style in `CoveragePlot`
* Allow passing additional arguments to `FilterFragments`
* Add inst/extdata 
* Change DESCRIPTION file so that Bioconductor dependencies are automatically installed

# Signac 0.1.6

* Bug fix for `GetCellsInRegion`
* Improve documentation

# Signac 0.1.5

* New `TSSEnrichment` and `TSSPlot` functions for TSS enrichment scoring
* New `InsertionBias` function
* New options in `CoveragePlot` for scaling tracks 
* Major speed improvements for `CoveragePlot` 
* Improved documentation (added examples)

# Signac 0.1.4

* Updates to `CoveragePlot`: now plots a Tn5 integration score per base, rather than the whole fragment.

# Signac 0.1.3

* New `GetIntersectingFeatures` function to find overlapping peaks between objects  
* New `MergeWithRegions` function to perform region-aware Seurat object merging  

# Signac 0.1.2

* New `RunChromVAR` function to run chromVAR through Signac  
* New `RegionStats` function to add statistics about peak sequences to the feature metadata  
* Improvements to `FindMotifs`: now selects a set of background peaks matching the sequence characteristics of the input

# Signac 0.1.1

* Added `IntersectMatrix`  
* Added unit tests  
* Bug fixes for `ChunkGRanges`

# Signac 0.1.0

* This is the first release of Signac!

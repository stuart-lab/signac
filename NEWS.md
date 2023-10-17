# Signac 1.11.0

Bug fixes:

* Fix error when supplying genome string to `AddMotifs()` function ([#1437](https://github.com/stuart-lab/signac/issues/1437))
* Fix error in quantification with multiple fragment files ([#1474](https://github.com/stuart-lab/signac/issues/1474))

Other changes:

* Added `region_extension` parameter to `TSSEnrichment()` function to enable changing the size of the region used to compute TSS enrichment scores ([#1444](https://github.com/stuart-lab/signac/pull/1402); @twmcart)
* Allow fragment files with a CSI index ([#1494](https://github.com/stuart-lab/signac/issues/1494))

# Signac 1.10.0

New functionality:

* Added `DensityScatter()` function
* Added `scale.linewidth` parameter to `LinkPlot()` ([#1412](https://github.com/stuart-lab/signac/issues/1412))

Bug fixes:

* Fixed error in `GetReadsInRegion()` when no fragments present that overlap the region (@nrockweiler; [#1348](https://github.com/stuart-lab/signac/pull/1348))
* Make motif names unique to avoid errors in motif-related functions ([#1311](https://github.com/stuart-lab/signac/issues/1311))
* Ignore missing seqnames when constructing a Motif object ([#1381](https://github.com/stuart-lab/signac/issues/1381))
* Fixed error in `BigwigTrack()` ([#1389](https://github.com/stuart-lab/signac/issues/1389))
* Fixed bug in `RegionMatrix()` where regions on minus strand were not handled correctly ([#1368](https://github.com/stuart-lab/signac/issues/1368))

Other changes:

* Improved documentation for `AddMotifs()` function ([#1361](https://github.com/stuart-lab/signac/issues/1361))

# Signac 1.9.0

Bug fixes:

* Fixed Motif object subsetting when using only one feature ([#1266](https://github.com/stuart-lab/signac/discussions/1266))
* Avoid error in chromVAR when matrix contains NA values ([#1254](https://github.com/stuart-lab/signac/issues/1254))
* Fixed error in `CoveragePlot()` when setting `show.bulk=TRUE` ([#1294](https://github.com/stuart-lab/signac/issues/1294))

Other changes:

* Improved error message when no fragment file information stored ([#1236](https://github.com/stuart-lab/signac/issues/1236))
* Add `tol` parameter to `RunSVD()` to control the `irlba::irlba()` `tol` parameter
* Added `split.by` parameters to `CoveragePlot()` and `PlotFootprint()` ([#523](https://github.com/stuart-lab/signac/issues/523))

# Signac 1.8.0

Bug fixes:

* Fixed bug in Motif object subsetting ([#1162](https://github.com/stuart-lab/signac/issues/1162))
* Fixed bug in `FeatureMatrix()` where cell names would not be converted correctly when running
with `cells=NULL` ([#1198](https://github.com/stuart-lab/signac/issues/1198))
* Fixed bug causing crash when fragment file line lengths are very large ([#723](https://github.com/stuart-lab/signac/issues/723))

Other changes:

* Improved fold change calculation ([#1174](https://github.com/stuart-lab/signac/issues/1174))
* Added `subset()` method for `Fragment`-class objects
* Set random seed within `CoverageTrack` to enable reproducible coverage plots ([#1206](https://github.com/stuart-lab/signac/issues/1206))
* Add ability to plot region using gene name in `PeakPlot()`, `LinkPlot()`, and 
`AnnotationPlot()` ([#919](https://github.com/stuart-lab/signac/issues/919))
* Change `Motif` class definition to allow any `CsparseMatrix` in the data slot

# Signac 1.7.0

New functionality:

* Added `method` parameter to `LinkPeaks()` ([@saketkc](https://github.com/saketkc); [#1030](https://github.com/stuart-lab/signac/pull/1030))
* Added ability to plot transcript isoforms (set `annotation="transcript"` in `CoveragePlot()`)
* Added multiple testing correction to `FindMotifs()`, and new `p.adjust` column in output dataframe
* Added `p.adjust.method` parameter to `FindMotifs()` to control multiple testing correction method used

Bug fixes:

* Fixed bug in `CallPeaks()` when project name contained whitespace ([#981](https://github.com/stuart-lab/signac/issues/981))
* Fixed bug in `CoveragePlot()` when `tile=TRUE`
* Fixed bug in `InsertionBias()` that set the coordinates beyond the end of some chromosomes ([#986](https://github.com/stuart-lab/signac/issues/986))
* Fixed bug in `BigwigTrack()` when supplying a single bigwig file ([#1053](https://github.com/stuart-lab/signac/issues/1053))
* Fixed bug in `GeneActivity()` when specifying biotypes ([#1058](https://github.com/stuart-lab/signac/issues/1058))
* Fixed bug in `GeneActivity()` when gene name is an empty string ([#1055](https://github.com/stuart-lab/signac/issues/1055))
* Fixed bug in `FeatureMatrix()` when using list of Fragment objects ([#1056](https://github.com/stuart-lab/signac/issues/1056))
* Fixed bug in `RegionMatrix()` when running on objects containing renamed cells ([#1076](https://github.com/stuart-lab/signac/issues/1076))
* Fixed bug in `Footprint()` when using a FASTA file ([#1092](https://github.com/stuart-lab/signac/issues/1092))
* Fixed bug in `Footprint()` when using list of genomic regions ([#1098](https://github.com/stuart-lab/signac/issues/1098))
* Add check for invalid features in `FindMotifs()` ([#1109](https://github.com/stuart-lab/signac/issues/1109))

Other changes:

* Improved error message when no fragment file is supplied in `CallPeaks()` ([#1062](https://github.com/stuart-lab/signac/issues/1062))
* Moved `Seurat`, `ggforce`, `ggrepel`, `ggseqlogo`, `lsa`, `qlcMatrix` to suggested packages
* Added progress bar to `GetGRangesFromEnsDb()`

# Signac 1.6.0

New functionality:

* Added `RegionHeatmap()` function
* Added `RegionPlot()` function
* Added `RegionMatrix()` function
* Added `MotifCounts()` function
* Added ability to plot multiple assays simultaneously in `CoveragePlot()` by providing a list of assay names
* Added `assay.scale`, `bigwig.scale`, and `split.assay` parameters to `CoveragePlot()`
* Added new option to provide quantile `ymax` value in `CoveragePlot()`
* Added `bigwig.scale` and `ymax` parameters to `BigwigTrack()`
* Added `peak.slot` parameter to `LinkPeaks()` ([#932](https://github.com/stuart-lab/signac/issues/932))

Bug fixes:

* Fixed behaviour of `min.features` parameter in `CreateChromatinAssay()` to retain cells with `>= min.features` (#[902](https://github.com/stuart-lab/signac/discussions/902))
* Fixed behaviour of `min.cells` parameter in `LinkPeaks()` ([#932](https://github.com/stuart-lab/signac/issues/932))

Other changes:

* Updated documentation for `CoveragePlot()` ([#992](https://github.com/stuart-lab/signac/issues/992))
* Removed `method` parameter from `LinkPeaks()` ([#932](https://github.com/stuart-lab/signac/issues/932))

# Signac 1.5.0

Bug fixes:

* Fixed bug in `FeatureMatrix()` when cells information not present in Fragment object ([#803](https://github.com/stuart-lab/signac/issues/803))
* Fixed bug in object merging ([#804](https://github.com/stuart-lab/signac/issues/804))
* Add ability to run `LinkPeaks()` using Ensembl IDs ([#858](https://github.com/stuart-lab/signac/issues/858))
* Fix issue in `GeneActivity()` when gene names are `NA` ([#865](https://github.com/stuart-lab/signac/issues/865))
* Fix bug in `FeatureMatrix()` when only one region supplied
* Allow negative values in `ExpressionPlot()` when using scaled data ([#893](https://github.com/stuart-lab/signac/issues/893))

Other changes:

* Added `idf` parameter to `RunTFIDF()` to use precomputed IDF vector
* Added `gene.id` parameter to `GeneActivity()` to allow output genes named using gene ID ([#837](https://github.com/stuart-lab/signac/issues/837))
* Added `sep` parameter to `ConnectionsToLinks()` ([#841](https://github.com/stuart-lab/signac/issues/841))

# Signac 1.4.0

New functionality:

* Added ability to display multiple assays in `CoveragePlot()`. The `assay` parameter can now be a list of assays to plot data from, with signal colored by assay of origin.

Bug fixes:

* Fixed bug in `FindMotifs()` when using only one region as input ([#732](https://github.com/stuart-lab/signac/issues/732))
* Add check for correct number of columns in fragment file ([#748](https://github.com/stuart-lab/signac/issues/748))
* Fixed gene lookup when annotations contain NA values ([#771](https://github.com/stuart-lab/signac/issues/771))
* Fixed error in `ClosestFeature()` when query contained regions on contigs not present in gene annotation ([#758](https://github.com/stuart-lab/signac/issues/758))
* Fixed bug in `TSSEnrichment()` when using multiple fragment files ([#783](https://github.com/stuart-lab/signac/issues/783))
* Fixed bug in `CallPeaks()` when multiple fragment files used as input
* Fixed bug in `CallPeaks()` to account for 0-based starts in called peaks
* Fixed bug in gene name lookup when gene names contain `-` characters ([#759](https://github.com/stuart-lab/signac/issues/759))

Other changes:

* Updated documentation for `genome` parameter in `AddMotifs()` and `RunChromVAR()` ([#712](https://github.com/stuart-lab/signac/issues/712))
* Updated the `FoldChange()` function to use normalized counts rather than raw counts ([#795](https://github.com/stuart-lab/signac/issues/795))
* Improved error checking in `GeneActivity()` ([#797](https://github.com/stuart-lab/signac/issues/797))
* Added `format` parameter to `CallPeaks()` ([#682](https://github.com/stuart-lab/signac/issues/682))

# Signac 1.3.0 

Bug fixes:

* Fixed `LinkPeaks()` function when running on a single gene ([#629](https://github.com/stuart-lab/signac/issues/629))
* Added `fragment.tempdir` parameter to `CallPeaks()` to enable setting directory
that split fragment files are written to during peak calling ([#579](https://github.com/stuart-lab/signac/issues/579))
* Fixed error in `FeatureMatrix()` when setting `sep` parameter ([#626](https://github.com/stuart-lab/signac/discussions/626))
* Fixed peak calling error when group names contain special characters
* Fixed issue with `RenameCells()` when cell information not present in Fragment object ([#704](https://github.com/stuart-lab/signac/issues/704))

Other changes: 

* Improved error checking for `GeneActivity()` ([#625](https://github.com/stuart-lab/signac/issues/625))
* Added `FoldChange()` method for `ChromatinAssay()` object that sets proper parameters for 
chromatin data. This fixes the calculation of fold changes when running `Seurat::FindMarkers()` on
single-cell chromatin data.

# Signac 1.2.1

New functionality:

* Added `head()` method for `Fragment`-class objects.

Bug fixes:

* Fixed bug in `ChromatinAssay` merging ([#596](https://github.com/stuart-lab/signac/pull/596))

Other changes:

* Added support for fragment files containing headers (cellranger-atac v2; [#609](https://github.com/stuart-lab/signac/issues/609))

# Signac 1.2.0

New functionality:

* Added `BigwigTrack()` function to plot data from bigWig files
* Added `bigwig` and `bigwig.type` arguments to `CoveragePlot()` to
include bigWig files in `CoveragePlot()`
* Added `region.highlight` parameter to `CoveragePlot()`
* Added `biotypes` parameter to `GeneActivity()` and `GetTSSPositions()` functions
* Added `max.width` parameter to `GeneActivity()`
* Added `min.distance` parameter to `LinkPeaks()` ([#561](https://github.com/stuart-lab/signac/pull/561))

Bug fixes:

* Fixed fragment file reading when only one fragment found in requested region ([#474](https://github.com/stuart-lab/signac/issues/474))
* Fixed `standard.chromosomes` parameter in `GetGRangesFromEnsDb()` ([#513](https://github.com/stuart-lab/signac/issues/513))
* Fixed `group.by` parameter in `PlotFootprint()` ([#522](https://github.com/stuart-lab/signac/issues/522))
* Fixed bug that would cause some gene coordinates used by `GeneActivity()` to be 
incorrect ([#521](https://github.com/stuart-lab/signac/issues/521))
* Fixed error message in `FindMotifs()` ([#549](https://github.com/stuart-lab/signac/issues/549))
* Fixed bug in `CountsInRegion()` ([#563](https://github.com/stuart-lab/signac/issues/563))

Other changes:

* Improved speed of ChromatinAssay merging
* Improved error message for `TSSEnrichment()` ([#485](https://github.com/stuart-lab/signac/issues/485))
* Improved error messages when trying to run `ChromatinAssay`-specific functions
on non-`ChromatinAssay` assays
* Performance improvements
* Changed default value for `n` in `NucleosomeSignal()`
* Enabled parallization in `TSSEnrichment()` when `fast=TRUE`
* Added early error checking in `LinkPeaks()` ([#550](https://github.com/stuart-lab/signac/pull/550))
* Change to sparse matrix correlation in `LinkPeaks()` ([#550](https://github.com/stuart-lab/signac/pull/550))
* Moved `biovizBase` and `Biostrings` to suggested packages
* Removed `ggbio` dependency
* Re-implemented `AnnotationPlot()`

# Signac 1.1.1 

New functionality:

* Added `group.by` parameter to `PeakPlot()` to allow coloring plotted genomic 
ranges by metadata variables.
* Added `peaks.group.by` and `ranges.group.by` parameters to `CoveragePlot()` to
allow coloring plotted genomic ranges in `CoveragePlot()` to be colored by metadata
variables.

Bug fixes:

* Update meta feature information (overall peak accessibility) when subsetting 
objects to avoid counts becoming inaccurate ([#332](https://github.com/stuart-lab/signac/issues/332))
* Prevent dropping features when creating a merged ChromatinAssay ([#340](https://github.com/stuart-lab/signac/pull/340))
* Fix compilation error when using g++ version <5 ([#326](https://github.com/stuart-lab/signac/issues/326))
* Retain motif positions during subset ([#364](https://github.com/stuart-lab/signac/issues/364))
* Fix `assay` parameter in `CoveragePlot()`
* Fix error when merging ChromatinAssay object ([#355](https://github.com/stuart-lab/signac/issues/355))
* Add more informative error message when all features or cells removed by parameter choices in `CreateChromatinAssay()` ([#387](https://github.com/stuart-lab/signac/issues/387))
* Fix bug in `CreateChromatinAssay()` when setting both `min.cells` and `min.features` arguments ([#390](https://github.com/stuart-lab/signac/issues/390))
* Improved support for remote fragment files
* Fixed bug in `PlotFootprint()` when only one cell in an identity class ([#406](https://github.com/stuart-lab/signac/issues/406))

Other changes:

* Added citation information to the package
* Added `SeuratObject` dependency

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
fragment objects ([#206](https://github.com/stuart-lab/signac/issues/206))

Bug fixes: 

* Fixed bugs in `AlleleFreq()` ([#196](https://github.com/stuart-lab/signac/issues/196)
and [#260](https://github.com/stuart-lab/signac/issues/260))
* Fixed bug in `FeatureMatrix()` ([#205](https://github.com/stuart-lab/signac/issues/205), [#291](https://github.com/stuart-lab/signac/issues/291))
* Fixed bug in `CreateChromatinAssay()` when setting `min.features` argument ([#194](https://github.com/stuart-lab/signac/issues/194))
* Fixed bug in `CreateChromatinAssay()` when setting `min.cells` argument ([#292](https://github.com/stuart-lab/signac/issues/292))
* Fixed bug in `TSSEnrichment()` when cell information not set for fragment files ([#203](https://github.com/stuart-lab/signac/issues/203))
* Fixed bug in `TSSEnrichment()` when no fragments present in TSS region ([#244](https://github.com/stuart-lab/signac/issues/244))
* Removed `qvalue` calculation from `FindMotifs()` ([#223](https://github.com/stuart-lab/signac/issues/223))
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

# Signac 1.0.0

The entire package has been updated to use a new `ChromatinAssay` class for the
storage of single-cell chromatin data. This is an extension of the standard 
Seurat `Assay` that adds additional slots needed for the analysis of chromatin 
data, including genomic ranges, genome information, fragment file information,
motifs, and gene annotations.

In addition, we have defined a new `Fragment` class to store information 
regarding the fragment file. This makes use of the fragment file more robust,
as checks are now performed to verify that the expected cells are present in
the fragment file, and that the fragment file or index are not modified on disk.

Key new functionality:

* Ability to store multiple fragment files: you can now store as many fragment
files as needed in a single object, and all functions that use the fragment file
will pull data from each of the files. Cell barcodes can also be modified in the 
object without needing to update the cell barcodes in the fragment file. The 
`Fragment` object stores an internal mapping of the cell barcode in the file to
the cell barcode in the Signac object.

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

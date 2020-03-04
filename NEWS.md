# Signac 0.2.4

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

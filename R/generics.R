#' Add DNA sequence motif information
#'
#' Construct a [Motif-class] object containing DNA sequence motif
#' information and add it to an existing Seurat object or ChromatinAssay.
#' If running on a Seurat object, `AddMotifs` will also run
#' [RegionStats()] to compute the GC content of each peak and store
#' the results in the feature metadata. PFMs or PWMs are matched to the genome
#' sequence using the [motifmatchr::matchMotifs()] function with
#' default parameters to construct a matrix of motif positions in genomic
#' regions.
#'
#' @param object A Seurat object or ChromatinAssay object
#' @param ... Additional arguments passed to other methods
#' @export AddMotifs
#' @rdname AddMotifs
#' @return When running on a `ChromatinAssay` or `Seurat` object,
#' returns a modified version of the input object. When running on a matrix,
#' returns a `Motif` object.
#' @seealso \pkg{motifmatchr}
AddMotifs <- function(object, ...) {
  UseMethod(generic = "AddMotifs", object = object)
}

#' Quantify aggregated genome tiles
#'
#' Quantifies fragment counts per cell in fixed-size genome bins across the
#' whole genome, then removes bins with less than a desired minimum number of
#' counts in the bin, then merges adjacent tiles into a single region.
#'
#' @param object A Seurat object or ChromatinAssay object
#' @param ... Additional arguments passed to other methods
#' @export AggregateTiles
#' @rdname AggregateTiles
AggregateTiles <- function(object, ...) {
  UseMethod(generic = "AggregateTiles", object = object)
}

#' Compute scATAC-seq QC metrics
#' 
#' Wrapper function to run `fragtk qc` and add the metadata to the Seurat
#' object.
#' 
#' @param object A [SeuratObject::Seurat] object, [ChromatinAssay5-class] object
#' , or path to a fragment file.
#' @param fragtk.path Path to `fragtk` executable. If NULL, try to find `fragtk`
#' automatically.
#' @param annotations [GenomicRanges::GRanges()] object containing
#' gene annotations. If `NULL`, attempt to extract this from the `ChromatinAssay5`
#' object if provided.
#' @param outdir Path for output directory
#' @param cleanup Remove output files created by fragtk
#' @param verbose Display messages
#' @param ... Arguments passed to other methods
#' 
#' @export ATACqc
#' @rdname ATACqc
ATACqc <- function(object, ...) {
  UseMethod(generic = "ATACqc", object = object)
}

#' Convert objects to a [ChromatinAssay5-class] object
#' 
#' @param x An object to convert to class [ChromatinAssay5-class]
#' @param ... Arguments passed to other methods
#' @rdname as.ChromatinAssay5
#' @export as.ChromatinAssay5
as.ChromatinAssay5 <- function(x, ...) {
  UseMethod(generic = "as.ChromatinAssay5", object = x)
}

#' Convert objects to a [GRangesAssay-class] object
#' 
#' @param x An object to convert to class [GRangesAssay-class]
#' @param ... Arguments passed to other methods
#' @rdname as.GRangesAssay
#' @export as.GRangesAssay
as.GRangesAssay <- function(x, ...) {
  UseMethod(generic = "as.GRangesAssay", object = x)
}

#' Convert objects to a [Fragment2-class] object
#' 
#' @param x An object to convert to class [Fragment2-class]
#' @param ... Arguments passed to other methods
#' @rdname as.Fragment2
#' @export as.Fragment2
as.Fragment2 <- function(x, ...) {
  UseMethod(generic = "as.Fragment2", object = x)
}

#' Compute allele frequencies per cell
#'
#' Collapses allele counts for each strand and normalize by the total number of
#' counts at each nucleotide position.
#'
#' @param object A [SeuratObject::Seurat] object, Assay, or matrix
#' @param variants A character vector of informative variants to keep. For
#' example, `c("627G>A","709G>A","1045G>A","1793G>A")`.
#' @param ... Arguments passed to other methods
#'
#' @export
#' @return Returns a [SeuratObject::Seurat] object with a new assay
#' containing the allele frequencies for the informative variants.
AlleleFreq <- function(object, ...) {
  UseMethod(generic = "AlleleFreq", object = object)
}

#' Annotation
#'
#' Get the annotation from a [ChromatinAssay5-class] object
#'
#' @param ... Arguments passed to other methods
#' @return Returns a [GenomicRanges::GRanges()] object
#' if the annotation data is present, otherwise returns NULL
#' @rdname Annotation
#' @export Annotation
Annotation <- function(object, ...) {
  UseMethod(generic = "Annotation", object = object)
}

#' @param value A value to set. Can be `NULL`, to remove the current annotation
#' information, or a [GenomicRanges::GRanges()] object.
#'
#' @rdname Annotation
#' @export Annotation<-
#'
"Annotation<-" <- function(object, ..., value) {
  UseMethod(generic = 'Annotation<-', object = object)
}

#' Bias
#'
#' Get the bias information from a [ChromatinAssay5-class] object
#'
#' @param ... Arguments passed to other methods
#' @return Returns a numeric vector or `NULL`
#' @rdname Bias
#' @export Bias
Bias <- function(object, ...) {
  UseMethod(generic = "Bias", object = object)
}

#' @param value A value to set. Can be `NULL`, to remove the current bias
#' information, or a numeric vector object.
#'
#' @rdname Bias
#' @export Bias<-
#'
"Bias<-" <- function(object, ..., value) {
  UseMethod(generic = 'Bias<-', object = object)
}

#' Binarize counts
#'
#' Set counts >1 to 1 in a count matrix
#'
#' @param object A [SeuratObject::Seurat] object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat] object
#' @rdname BinarizeCounts
#' @export BinarizeCounts
BinarizeCounts <- function(object, ...) {
  UseMethod(generic = "BinarizeCounts", object = object)
}

#' Call peaks
#'
#' Call peaks using `MACS`. Fragment files linked to the specified assay will be
#' used to call peaks. If multiple fragment files are present, all will be used
#' in a single MACS invocation. Returns the `.narrowPeak` MACS output as a
#' `GRanges` object.
#'
#' See <https://macs3-project.github.io/MACS/> for MACS documentation.
#'
#' If you call peaks using MACS2 please cite:
#' \doi{10.1186/gb-2008-9-9-r137}
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a [GenomicRanges::GRanges()] object
#' @rdname CallPeaks
#' @export CallPeaks
CallPeaks <- function(object, ...) {
  UseMethod(generic = "CallPeaks", object = object)
}

#' Set and get cell barcode information for a [Fragment2-class] object
#'
#' @param x A Seurat object
#' @param value A character vector of cell barcodes
#' @param ... Arguments passed to other methods
#' @export Cells<-
#' @concept assay
"Cells<-" <- function(x, ..., value) {
  UseMethod(generic = "Cells<-", object = x)
}

#' Convert between motif name and motif ID
#'
#' Converts from motif name to motif ID or vice versa. To convert common names
#' to IDs, use the `name` parameter. To convert IDs to common names, use
#' the `id` parameter.
#'
#' @param object A [SeuratObject::Seurat], [ChromatinAssay5-class], or
#' [Motif-class] object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a character vector with the same length and order as the
#' input. Any names or IDs that were not found will be stored as `NA`.
#'
#' @rdname ConvertMotifID
#' @export ConvertMotifID
#'
ConvertMotifID <- function(object, ...) {
  UseMethod(generic = "ConvertMotifID", object = object)
}

#' Find variable features fitting LOESS
#'
#' Find highly variable features by fitting a locally polynomial regression
#' model to the log(mean) and log(variance) of downsampled features.
#' 
#' This function is similar to the [Seurat::FindVariableFeatures()]
#' function (with `selection.method="vst"`), but downsamples the features
#' evenly across the range of mean values. This speeds up fitting of the loess 
#' curve when the number of features is large.
#' 
#' The function also provides the ability to combine ranking of features
#' according to their mean count and their residual variance, using a weighted
#' rank sum with weights set by the `weight.mean` parameter. This can help
#' to avoid selecting features with high residual variance but very low mean.
#'
#' @param object A [SeuratObject::Seurat] object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname FitMeanVar
#' @export FitMeanVar
FitMeanVar <- function(object, ...) {
  UseMethod(generic = "FitMeanVar", object = object)
}

#' Find most frequently observed features
#'
#' Find top features for a given assay based on total number of counts for the
#' feature. Can specify a minimum cell count, or a lower percentile
#' bound to determine the set of variable features. Running this function will
#' store the total counts and percentile rank for each feature in the feature
#' metadata for the assay. To only compute the feature metadata, without
#' changing the variable features for the assay, set `min.cutoff=NA`.
#'
#' @param object A [SeuratObject::Seurat] object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname FindTopFeatures
#' @export FindTopFeatures
FindTopFeatures <- function(object, ...) {
  UseMethod(generic = "FindTopFeatures", object = object)
}

#' Transcription factor footprinting analysis
#'
#' Compute the normalized observed/expected Tn5 insertion frequency
#' for each position surrounding a set of motif instances.
#'
#' @param object A [SeuratObject::Seurat] or [ChromatinAssay5-class] object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat] object
#' @rdname Footprint
#' @export Footprint
Footprint <- function(object, ...) {
  UseMethod(generic = "Footprint", object = object)
}

#' Get the Fragment objects
#'
#' @param ... Arguments passed to other methods
#' @return Returns a list of [Fragment2-class] objects. If there are
#' no [Fragment2-class] objects present, returns an empty list.
#' @rdname Fragments
#' @export Fragments
Fragments <- function(object, ...) {
  UseMethod(generic = "Fragments", object = object)
}

#' @param value A [Fragment2-class] object or list of Fragment objects
#'
#' @rdname Fragments
#' @export Fragments<-
#'
"Fragments<-" <- function(object, ..., value) {
  UseMethod(generic = 'Fragments<-', object = object)
}

#' Compute Tn5 insertion bias
#'
#' Counts the Tn5 insertion frequency for each DNA hexamer.
#' @param object A [SeuratObject::Seurat] or [ChromatinAssay5-class] object
#' @param ... Arguments passed to other methods
#' @return Returns a Seurat object
#' @rdname InsertionBias
#' @export InsertionBias
InsertionBias <- function(object, ...) {
  UseMethod(generic = "InsertionBias", object = object)
}


#' Get genes linked to peaks
#' 
#' Retrieve peak-gene links for a given set of genes. Links must be first
#' obtained by running the [LinkPeaks()] function.
#' 
#' This function is designed to obtain the stored results from running the
#' [LinkPeaks()] function. Alternatively, custom peak-gene linkage methods
#' can be used as long as they store the gene name, peak name, and a peak-gene
#' score information as metadata columns named "gene," "peak," and "score"
#' respectively.
#' @param object A [SeuratObject::Seurat] object
#' @param ... Arguments passed to other methods
#' @return Returns a character vector of peaks
#'
#' @export GetLinkedGenes
#' @rdname GetLinkedGenes
#' @seealso GetLinkedPeaks
GetLinkedGenes <- function(object, ...) {
  UseMethod(generic = "GetLinkedGenes", object = object)
}

#' Get peaks linked to genes
#' 
#' Retrieve peak-gene links for a given set of genes. Links must be first
#' obtained by running the [LinkPeaks()] function.
#' 
#' This function is designed to obtain the stored results from running the
#' [LinkPeaks()] function. Alternatively, custom peak-gene linkage methods
#' can be used as long as they store the gene name, peak name, and a peak-gene
#' score information as metadata columns named "gene," "peak," and "score"
#' respectively.
#'
#' @param object A [SeuratObject::Seurat] object
#' @param ... Arguments passed to other methods
#' @return Returns a character vector of genes
#' 
#' @export GetLinkedPeaks
#' @rdname GetLinkedPeaks
#' @seealso GetLinkedGenes
GetLinkedPeaks <- function(object, ...) {
  UseMethod(generic = "GetLinkedPeaks", object = object)
}

#' Retrieve a motif matrix
#'
#' Get motif matrix for given assay
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname GetMotifData
#' @export GetMotifData
GetMotifData <- function(object, ...) {
  UseMethod(generic = "GetMotifData", object = object)
}

#' Get or set a motif information
#'
#' Get or set the [Motif-class] object for a Seurat object or
#' [ChromatinAssay5-class].
#'
#' @param ... Arguments passed to other methods
#' @rdname Motifs
#' @export Motifs
Motifs <- function(object, ...) {
  UseMethod(generic = "Motifs", object = object)
}

#' @param value A [Motif-class] object
#' @rdname Motifs
#' @export Motifs<-
"Motifs<-" <- function(object, ..., value) {
  UseMethod(generic = 'Motifs<-', object = object)
}

#' Get or set links information
#'
#' Get or set the genomic link information for a Seurat object or
#' [ChromatinAssay5-class]
#'
#' @param ... Arguments passed to other methods
#' @rdname Links
#' @export Links
Links <- function(object, ...) {
  UseMethod(generic = "Links", object = object)
}

#' @param value A [GenomicRanges::GRanges()] object
#' @rdname Links
#' @export Links<-
"Links<-" <- function(object, ..., value) {
  UseMethod(generic = "Links<-", object = object)
}

#' Identify mitochondrial variants
#'
#' Identify mitochondrial variants present in single cells.
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @rdname IdentifyVariants
#' @export IdentifyVariants
IdentifyVariants <- function(object, ...) {
  UseMethod(generic = "IdentifyVariants", object = object)
}

#' Compute analytic Pearson residual variance
#' 
#' Find the top features for a given assay based on analytic Pearson residual
#' variance. This function computes the Pearson residual variance for each
#' feature without constructing the entire dense matrix of Pearson residuals to
#' reduce the memory required.
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname PearsonResidualVar
#' @export PearsonResidualVar
#' @references
#' Lause, J., Berens, P. & Kobak, D. Analytic Pearson residuals for
#' normalization of single-cell RNA-seq UMI data. Genome Biol 22, 258 (2021).
#' <https://doi.org/10.1186/s13059-021-02451-7>
PearsonResidualVar <- function(object, ...) {
  UseMethod(generic = "PearsonResidualVar", object = object)
}

#' Region Aggregation
#'
#' Get the region aggregation information from a [ChromatinAssay5-class] object
#' 
#' @param object A Seurat object 
#' @param ... Arguments passed to other methods
#' @return Returns a list of [RegionAggregation-class] objects
#' @rdname RegionAggr
#' @export RegionAggr
RegionAggr <- function(object, ...) {
  UseMethod(generic = "RegionAggr", object = object)
}

#' @param value A [RegionAggregation-class] object or list of [RegionAggregation-class] objects
#'
#' @rdname RegionAggr
#' @export RegionAggr<-
#'
"RegionAggr<-" <- function(object, ..., value) {
  UseMethod(generic = 'RegionAggr<-', object = object)
}

#' Region enrichment analysis
#'
#' Count fragments within a set of regions for different groups of
#' cells.
#'
#' @param object A Seurat or [ChromatinAssay5-class] object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname RegionMatrix
#' @export RegionMatrix
RegionMatrix <- function(object, ...) {
  UseMethod(generic = "RegionMatrix", object = object)
}

#' Compute base composition information for genomic ranges
#'
#' Compute the GC content, region lengths, and dinucleotide base frequencies
#' for regions in the assay and add to the feature metadata.
#'
#' @param object A Seurat object, Assay object, or set of genomic ranges
#' @param ... Arguments passed to other methods
#' @return Returns a dataframe
#' @rdname RegionStats
#' @export RegionStats
RegionStats <- function(object, ...) {
  UseMethod(generic = "RegionStats", object = object)
}

#' Run chromVAR
#'
#' Wrapper to run [chromVAR::chromVAR()] on an assay with a motif
#' object present. Will return a new Seurat assay with the motif activities
#' (the deviations in chromatin accessibility across the set of regions) as
#' a new assay.
#'
#' See the chromVAR documentation for more information:
#' <https://greenleaflab.github.io/chromVAR/index.html>
#'
#' See the chromVAR paper: <https://www.nature.com/articles/nmeth.4401>
#'
#' @param object A Seurat object
#' @param genome A `BSgenome` object or string stating the genome build
#' recognized by `getBSgenome`.
#' @param motif.matrix A peak x motif matrix. If NULL, pull the peak x motif
#' matrix from a Motif object stored in the assay.
#' @param verbose Display messages
#' @param ... Additional arguments passed to
#' [chromVAR::getBackgroundPeaks()]
#' @return Returns a [SeuratObject::Seurat()] object with a new assay
#' @rdname RunChromVAR
#' @export RunChromVAR
RunChromVAR <- function(object, ...) {
  UseMethod(generic = "RunChromVAR", object = object)
}

#' Run singular value decomposition
#'
#' Run partial singular value decomposition using [RSpectra::svds()]
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname RunSVD
#' @export RunSVD
RunSVD <- function(object, ...) {
  UseMethod(generic = "RunSVD", object = object)
}

#' Compute the term-frequency inverse-document-frequency
#'
#' Run term frequency inverse document frequency (TF-IDF) normalization on a
#' matrix.
#'
#' Four different TF-IDF methods are implemented. We recommend using method 1
#' (the default).
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname RunTFIDF
#' @export RunTFIDF
#' @references <https://en.wikipedia.org/wiki/Latent_semantic_analysis#Latent_semantic_indexing>
RunTFIDF <- function(object, ...) {
  UseMethod(generic = "RunTFIDF", object = object)
}

#' Set motif data
#'
#' Set motif matrix for given assay
#'
#' @param object A Seurat object
#' @param ... Arguments passed to other methods
#' @return Returns a [SeuratObject::Seurat()] object
#' @rdname SetMotifData
#' @export SetMotifData
SetMotifData <- function(object, ...) {
  UseMethod(generic = "SetMotifData", object = object)
}

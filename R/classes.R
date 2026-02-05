#' @importFrom methods setClass setClassUnion
#' slotNames
#' @importClassesFrom Matrix CsparseMatrix
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom InteractionSet GInteractions
#' @useDynLib Signac
NULL


setClassUnion(name = "AnyMatrix", members = c("matrix", "CsparseMatrix"))
setClassUnion(name = "MatrixOrNULL", members = c("AnyMatrix", "NULL"))
setClassUnion(name = "GRangesOrNULL", members = c("GRanges", "NULL"))
setClassUnion(name = "NumericOrNULL", members = c("numeric", "NULL"))

#' The Fragment class
#'
#' The Fragment class is designed to hold information needed for working with
#' fragment files.
#'
#' @slot file.path Path to the fragment file on disk.
#' See <https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments>
#' @slot file.index Path to the fragment file index on disk.
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file, and the name of each element is the
#' corresponding cell barcode as stored in the ChromatinAssay5 object.
#' @slot seqlevels A named vector of sequence levels (eg, chromosome name) where
#' each element is the sequence name as it appears in the fragment file, and then
#' name of each element is the corresponding sequence name as stored in the 
#' ChromatinAssay5 object.
#'
#' @name Fragment2-class
#' @rdname Fragment2-class
#' @aliases Fragment
#' @exportClass Fragment2
#' @concept fragments
Fragment2 <- setClass(
  Class = "Fragment2",
  slots = list(
    file.path = "character",
    file.index = "character",
    hash = "character",
    cells = "ANY",
    seqlevels = "ANY"
  )
)

#' The Motif class
#'
#' The Motif class is designed to store DNA sequence motif information,
#' including motif PWMs or PFMs, motif positions, and metadata.
#'
#' @slot data A feature x motif matrix. Columns
#' correspond to motif IDs, rows correspond to genomic features
#' (peaks or bins). Entries in the matrix should be 1 if the
#' genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A named list of position weight matrices
#' @slot motif.names A list containing the name of each motif
#' @slot positions A [GenomicRanges::GRangesList()] object containing
#' exact positions of each motif.
#' @slot meta.data A dataframe for storage of additional
#' information related to each motif. This could include the
#' names of proteins that bind the motif.
#'
#' @name Motif-class
#' @rdname Motif-class
#' @exportClass Motif
#' @concept motifs
Motif <- setClass(
  Class = "Motif",
  slots = list(
    data = "MatrixOrNULL",
    pwm = "list",
    motif.names = "list",
    positions = "ANY",
    meta.data = "data.frame"
  )
)
setValidity(Class = "Motif", function(object) {
  if (length(x = object@positions) > 0 &&
      !all(vapply(
        X = object@positions,
        FUN = function(x) inherits(x = x, what = "GRanges"),
        logical(1)))) {
    return("All elements of 'positions' must be GRanges objects")
  }
  TRUE
})
setClassUnion(name = "MotifOrNULL", members = c("Motif", "NULL"))

#' RegionAggregation class
#' 
#' The RegionAggregation class enables storage of counts centered on a group of
#' genomic regions, aggregated across regions for each cell. The main data
#' object stored is a cell-by-position matrix.
#' 
#' @slot matrix A cell-by-position matrix
#' @slot regions A [GenomicRanges::granges()] object containing the
#' regions aggregated across.
#' @slot upstream Integer denoting number of bases upstream of the centered
#' position that are stored in the matrix
#' @slot downstream Integer denoting number of bases downstream of the centered
#' position that are stored in the matrix
#' @slot name A name for the set of regions aggregated
#' @slot expected  A vector containing expected number of Tn5 insertions per
#' position
#' @slot cells A vector of cell barcodes present in the region aggregation
#' matrix. Each entry in the cells vector corresponds to a row in the region
#' aggregation matrix.
#'
#' @name RegionAggregation-class
#' @rdname RegionAggregation-class
#' @exportClass RegionAggregation
#' @concept footprinting
RegionAggregation <- setClass(
  Class = "RegionAggregation",
  slots = list(
    matrix = "AnyMatrix",
    regions = "GRanges",
    upstream = "integer",
    downstream = "integer",
    name = "character",
    expected = "numeric",
    cells = "character"
  )
)
setValidity(Class = "RegionAggregation", function(object) {
  # make sure cell names do not contain NA
  if (any(is.na(x = object@cells))){
    return("Cells slot must not contain NA values")
  }
  # matrix rows must match cells number 
  if (dim(x = object@matrix)[1] != length(x = object@cells)){
    return("Number of rows in matrix must match length of cells")
  }
  # region width must be identical
  w <- unique(x = width(object@regions))
  if (length(x = w) != 1) {
    return("All regions must have identical width")  
  }
  # region width must match matrix columns 
  if (dim(x = object@matrix)[2] != (object@upstream + object@downstream + w)){
    return("Matrix columns do not match upstream + downstream + region width")
  }
  # expected vector length must match matrix columns 
  if (length(x = object@expected) != dim(x = object@matrix)[2]){
    return("Expected vector length must match number of matrix columns")
  }
  TRUE
})
setClassUnion(name = "ListOrNULL", members = c("list", "NULL"))


#' @slot fragments A list of [Fragment()] objects.
#' @slot annotation A  [GenomicRanges::GRanges()] object containing
#' genomic annotations. This should be a GRanges object with the following 
#' columns:
#' \itemize{
#'   \item{tx_id: Transcript ID}
#'   \item{gene_name: Gene name}
#'   \item{gene_id: Gene ID}
#'   \item{gene_biotype: Gene biotype (e.g. "protein_coding", "lincRNA")}
#'   \item{type: Annotation type (e.g. "exon", "gap")}
#' }
#' @slot bias A vector containing Tn5 integration bias information
#' (frequency of Tn5 integration at different kmers). This must be a named 
#' numeric vector where the name of the element corresponds to the DNA sequence
#' and the value represents the bias value. All DNA hexamers must be present in
#' the vector.
#' @slot region.aggregation A list of [RegionAggregation()] objects
#' @slot motifs A [Motif()] object
#' @slot links A list of [InteractionSet::GInteractions()] objects
#' describing linked genomic positions, such as co-accessible sites, eQTLs, 
#' Hi-C contact, or enhancer-gene regulatory relationships.
#'
#' @name ChromatinAssay5-class
#' @rdname GRangesAssay-class
#' @importClassesFrom SeuratObject Assay5
#' @exportClass ChromatinAssay5
#' @concept assay
ChromatinAssay5 <- setClass(
  Class = "ChromatinAssay5",
  contains = "Assay5",
  slots = list(
    "fragments" = "list",
    "annotation" = "GRangesOrNULL",
    "bias" = "NumericOrNULL",
    "region.aggregation" = "ListOrNULL",
    "links" = "list",
    "motifs" = "MotifOrNULL"
  )
)

setValidity(Class = "ChromatinAssay5", function(object) {
  if (length(x = object@links) > 0 &&
      !all(vapply(
        X = object@links,
        FUN = function(x) inherits(x = x, what = "GInteractions"),
        logical(1)))) {
    return("All elements of 'links' must be GInteractions objects")
  }
  if (length(x = object@fragments) > 0 &&
      !all(vapply(
        X = object@fragments,
        FUN = function(x) inherits(x = x, what = "Fragment2"),
        logical(1)))) {
    return("All elements of 'fragments' must be Fragment2 objects")
  }
  if (length(x = object@region.aggregation) > 0 &&
      !all(vapply(
        X = object@region.aggregation,
        FUN = function(x) inherits(x = x, what = "RegionAggregation"),
        logical(1)))) {
    return("All elements of 'region.aggregation' must be RegionAggregation objects")
  }
  if (!is.null(x = object@bias)) {
    if (!is.numeric(x = object@bias)) {
      return("Bias must be a numeric vector")
    }
    if (is.null(x = names(x = object@bias))) {
      return("Bias must be a named numeric vector")
    }
    bases <- c("A","C","G","T")
    hexamers <- apply(expand.grid(rep(list(bases), 6)), 1, paste0, collapse = "")
    if (!all(hexamers %in% names(x = object@bias))) {
      return("Bias vector must contain each hexamer")
    }
  }
  TRUE
})

#' ChromatinAssay5 and GRangesAssay object classes
#'
#' The `GRangesAssay` and `ChromatinAssay5` classes are extended
#' [SeuratObject::Assay5()] classes for the storage and analysis of
#' single-cell chromatin data. The `GRangesAssay` class requires that
#' features in the assay are genomic ranges.
#'
#' @slot ranges A [GenomicRanges::GRanges()] object describing the
#' genomic location of features in the object
#' 
#' @name GRangesAssay-class
#' @rdname GRangesAssay-class
#' @exportClass GRangesAssay
#' @concept assay
#' @seealso [SeuratObject::Assay5()]
GRangesAssay <- setClass(
  Class = "GRangesAssay",
  contains = "ChromatinAssay5",
  slots = list(
    "ranges" = "GRanges"
  )
)

##### OLD CLASSES #####

#' The Fragment class (old)
#'
#' The Fragment class is designed to hold information needed for working with
#' fragment files.
#'
#' @slot path Path to the fragment file on disk.
#' See <https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments>
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file, and the name of each element is the
#' corresponding cell barcode as stored in the ChromatinAssay object.
#'
#' @rdname oldfragment-class
#' @keywords internal
#' @noRd
Fragment <- setClass(
  Class = "Fragment",
  slots = list(
    path = "character",
    hash = "character",
    cells = "ANY"
  )
)

#' The ChromatinAssay class (old)
#'
#' The ChromatinAssay object is an extended [SeuratObject::Assay()]
#' for the storage and analysis of single-cell chromatin data.
#' 
#' This is an old object class used in Signac v1
#'
#' @slot ranges A [GenomicRanges::GRanges()] object describing the
#' genomic location of features in the object
#' @slot motifs A [Motif()] object
#' @slot fragments A list of [Fragment()] objects.
#' @slot seqinfo A [Seqinfo::Seqinfo()] object containing basic
#' information about the genome sequence used.
#' @slot annotation A  [GenomicRanges::GRanges()] object containing
#' genomic annotations. This should be a GRanges object with the following 
#' columns:
#' \itemize{
#'   \item{tx_id: Transcript ID}
#'   \item{gene_name: Gene name}
#'   \item{gene_id: Gene ID}
#'   \item{gene_biotype: Gene biotype (e.g. "protein_coding", "lincRNA")}
#'   \item{type: Annotation type (e.g. "exon", "gap")}
#' }
#' @slot bias A vector containing Tn5 integration bias information
#' (frequency of Tn5 integration at different kmers)
#' @slot positionEnrichment A named list of matrices containing positional
#' enrichment scores for Tn5 integration (for example, enrichment at the TSS)
#' @slot links A [GenomicRanges::GRanges()] object describing linked
#' genomic positions, such as co-accessible sites or enhancer-gene regulatory
#' relationships. This should be a `GRanges` object, where the start and
#' end coordinates are the two linked genomic positions, and must contain a
#' "score" metadata column.
#'
#' @name chromatinassay-class
#' @noRd
#' @importClassesFrom SeuratObject Assay
ChromatinAssay <- setClass(
  Class = "ChromatinAssay",
  contains = "Assay",
  slots = list(
    "ranges" = "GRanges",
    "motifs" = "ANY",
    "fragments" = "list",
    "seqinfo" = "ANY",
    "annotation" = "ANY",
    "bias" = "ANY",
    "region.aggregation" = "list",
    "links" = "GRanges"
  )
)

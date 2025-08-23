#' @importFrom methods setClass setClassUnion
#' slotNames
#' @importClassesFrom Matrix dgCMatrix
#' @useDynLib Signac
NULL


setClassUnion(name = "AnyMatrix", c("matrix", "dgCMatrix"))

#' The Fragment class
#'
#' The Fragment class is designed to hold information needed for working with
#' fragment files.
#'
#' @slot file.path Path to the fragment file on disk.
#' See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @slot file.index Path to the fragment file index on disk.
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file, and the name of each element is the
#' corresponding cell barcode as stored in the ChromatinAssay5 object.
#' @slot seqnames A named vector of sequence names (eg, chromosome name) where
#' each element is the sequence name as it appears in the fragment file, and the
#' name of each element is the corresponding sequence name as stored in the 
#' ChromatinAssay5 object.
#'
#' @name Fragment-class
#' @rdname Fragment-class
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
    seqnames = "ANY"
  )
)

#' The Motif class
#'
#' The Motif class is designed to store DNA sequence motif information,
#' including motif PWMs or PFMs, motif positions, and metadata.
#'
#' @slot data A sparse, binary, feature x motif matrix. Columns
#' correspond to motif IDs, rows correspond to genomic features
#' (peaks or bins). Entries in the matrix should be 1 if the
#' genomic feature contains the motif, and 0 otherwise.
#' @slot pwm A named list of position weight matrices
#' @slot motif.names A list containing the name of each motif
#' @slot positions A \code{\link[GenomicRanges]{GRangesList}} object containing
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
    data = "CsparseMatrix",
    pwm = "list",
    motif.names = "list",
    positions = "ANY",
    meta.data = "data.frame"
  )
)

#' The ChromatinAssay5 class
#'
#' The ChromatinAssay5 object is an extended \code{\link[SeuratObject]{Assay5}}
#' for the storage and analysis of single-cell chromatin data. This class does
#' not require that the stored features are genomic ranges.
#'
#' @slot fragments A list of \code{\link{Fragment}} objects.
#' @slot seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome sequence used.
#' @slot annotation A  \code{\link[GenomicRanges]{GRanges}} object containing
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
#'
#' @name ChromatinAssay5-class
#' @rdname ChromatinAssay5-class
#' @importClassesFrom SeuratObject Assay5
#' @exportClass ChromatinAssay5
#' @concept assay
ChromatinAssay5 <- setClass(
  Class = "ChromatinAssay5",
  contains = "Assay5",
  slots = list(
    "fragments" = "list",
    "annotation" = "ANY",
    "bias" = "ANY",
    "positionEnrichment" = "list"
  )
)

#' The GRangesAssay class
#'
#' The GRangesAssay object is an extended \code{\link{ChromatinAssay5}}
#' for the storage and analysis of single-cell chromatin data with associated
#' genomic ranges.
#'
#' @slot ranges A \code{\link[GenomicRanges]{GRanges}} object describing the
#' genomic location of features in the object
#' @slot motifs A \code{\link{Motif}} object
#' @slot links A \code{\link[GenomicRanges]{GRanges}} object describing linked
#' genomic positions, such as co-accessible sites or enhancer-gene regulatory
#' relationships. This should be a \code{GRanges} object, where the start and
#' end coordinates are the two linked genomic positions, and must contain a
#' "score" metadata column.
#'
#' @name GRangesAssay-class
#' @rdname GRangesAssay-class
#' @exportClass GRangesAssay
#' @concept assay
GRangesAssay <- setClass(
  Class = "GRangesAssay",
  contains = "ChromatinAssay5",
  slots = list(
    "ranges" = "GRanges",
    "motifs" = "Motif",
    "links" = "GRanges"
  )
)

##### OLD CLASSES #####

#' The Fragment class (old)
#'
#' The Fragment class is designed to hold information needed for working with
#' fragment files.
#'
#' @slot path Path to the fragment file on disk.
#' See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @slot hash A vector of two md5sums: first element is the md5sum of the
#' fragment file, the second element is the md5sum of the index.
#' @slot cells A named vector of cells where each element is the cell barcode
#' as it appears in the fragment file, and the name of each element is the
#' corresponding cell barcode as stored in the ChromatinAssay object.
#'
#' @rdname oldfragment-class
#' @keywords internal
#' @concept v1
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
#' The ChromatinAssay object is an extended \code{\link[SeuratObject]{Assay}}
#' for the storage and analysis of single-cell chromatin data.
#' 
#' This is an old object class used in Signac v1
#'
#' @slot ranges A \code{\link[GenomicRanges]{GRanges}} object describing the
#' genomic location of features in the object
#' @slot motifs A \code{\link{Motif}} object
#' @slot fragments A list of \code{\link{Fragment}} objects.
#' @slot seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome sequence used.
#' @slot annotation A  \code{\link[GenomicRanges]{GRanges}} object containing
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
#' @slot links A \code{\link[GenomicRanges]{GRanges}} object describing linked
#' genomic positions, such as co-accessible sites or enhancer-gene regulatory
#' relationships. This should be a \code{GRanges} object, where the start and
#' end coordinates are the two linked genomic positions, and must contain a
#' "score" metadata column.
#'
#' @name chromatinassay-class
#' @rdname oldchromatinassay-class
#' @importClassesFrom SeuratObject Assay
#' @concept v1
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
    "positionEnrichment" = "list",
    "links" = "GRanges"
  )
)
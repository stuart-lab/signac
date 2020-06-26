#' @include generics.R
#' @importFrom methods callGeneric
#' @importFrom Seurat DefaultAssay
#' @importFrom IRanges precede follow nearest distance distanceToNearest
#' findOverlaps countOverlaps coverage
#' reduce disjoin gaps isDisjoint disjointBins
NULL

setOldClass(Classes = "ChromatinAssay")

## Nearest methods

# precede

#' Find the nearest range neighbors for ChromatinAssay objects
#'
#' The \code{precede, follow, nearest, distance, distanceToNearest} methods
#' are available for \code{\link{ChromatinAssay}} objects.
#'
#' @name nearest-methods
#' @param x A query \code{\link{ChromatinAssay}} object
#' @param subject The subject \code{\link[GenomicRanges]{GRanges}} or
#' \code{\link{ChromatinAssay}} object. If missing, \code{x} is used as the
#' subject.
#' @param select Logic for handling ties.
#' See \code{\link[GenomicRanges]{nearest-methods}} in the \pkg{GenomicRanges}
#' package.
#' @param ignore.strand Logical argument controlling whether strand information
#' should be ignored.
#' @param ... Additional arguments for methods
#'
#' @aliases precede precede,ANY,ChromatinAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{nearest-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{nearest-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod precede
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("ANY", "ChromatinAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ANY
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("ChromatinAssay", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ChromatinAssay
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ANY, Seurat
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("ANY", "Seurat"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, ANY
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("Seurat", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, Seurat
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("Seurat", "Seurat"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay.s <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay.s]])
    assay.x <- DefaultAssay(object = x)
    x <- granges(x = x[[assay.x]])
    callGeneric()
  }
)

# follow

#' @aliases follow
#' @describeIn nearest-methods method for ANY, ChromatinAssay
#' @exportMethod follow
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("ANY", "ChromatinAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ANY
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("ChromatinAssay", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ChromatinAssay
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ANY, Seurat
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("ANY", "Seurat"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, ANY
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("Seurat", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, Seurat
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("Seurat", "Seurat"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay.s <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay.s]])
    assay.x <- DefaultAssay(object = x)
    x <- granges(x = x[[assay.x]])
    callGeneric()
  }
)

#' @aliases nearest
#' @describeIn nearest-methods method for ANY, ChromatinAssay
#' @exportMethod nearest
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("ANY", "ChromatinAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ANY
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("ChromatinAssay", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ChromatinAssay
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ANY, Seurat
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("ANY", "Seurat"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, ANY
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("Seurat", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, Seurat
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("Seurat", "Seurat"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    assay.s <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay.s]])
    assay.x <- DefaultAssay(object = x)
    x <- granges(x = x)
    callGeneric()
  }
)

#' @param y For the \code{distance} method, a
#' \code{\link[GenomicRanges]{GRanges}} object or a \code{\link{ChromatinAssay}}
#' object
#'
#' @aliases distance
#' @exportMethod distance
#' @describeIn nearest-methods method for ANY, ChromatinAssay
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("ANY", "ChromatinAssay"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    y <- granges(x = y)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ANY
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("ChromatinAssay", "ANY"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ChromatinAssay
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    y <- granges(x = y)
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ANY, Seurat
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("ANY", "Seurat"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    assay <- DefaultAssay(object = y)
    y <- granges(x = y[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, ANY
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("Seurat", "ANY"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, Seurat
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("Seurat", "Seurat"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    assay.x <- DefaultAssay(object = x)
    x <- granges(x = x[[assay.x]])
    assay.y <- DefaultAssay(object = y)
    y <- granges(x = y[[assay.y]])
    callGeneric()
  }
)

#' @aliases distanceToNearest
#' @exportMethod distanceToNearest
#' @describeIn nearest-methods method for ANY, ChromatinAssay
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("ANY", "ChromatinAssay"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ANY
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("ChromatinAssay", "ANY"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ChromatinAssay, ChromatinAssay
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    subject <- granges(x = subject)
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for ANY, Seurat
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("ANY", "Seurat"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    assay <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, ANY
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("Seurat", "ANY"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @describeIn nearest-methods method for Seurat, Seurat
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("Seurat", "Seurat"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    assay.s <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay.s]])
    assay.x <- DefaultAssay(object = x)
    x <- granges(x = x[[assay.x]])
    callGeneric()
  }
)

## Find overlaps methods

#' Find overlapping ranges for ChromatinAssay objects
#'
#' The \code{findOverlaps, countOverlaps} methods are available for
#' \code{\link{ChromatinAssay}} objects. This allows finding overlaps between
#' genomic ranges and the ranges stored in the ChromatinAssay.
#'
#' If a ChromatinAssay is set as the default assay in a
#' \code{\link[Seurat]{Seurat}} object, you can also call \code{findOverlaps}
#' directly on the Seurat object.
#'
#' @param query,subject A \code{\link{ChromatinAssay}} object
#' @param maxgap,minoverlap,type,select,ignore.strand See
#' \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} and
#' \pkg{IRanges} packages.
#' @return See \code{\link[GenomicRanges]{findOverlaps}}
#'
#' @name findOverlaps-methods
#' @aliases findOverlaps findOverlaps,Vector,ChromatinAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{findOverlaps-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{findOverlaps-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#'
#' @exportMethod findOverlaps
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("Vector", "ChromatinAssay"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for ChromatinAssay, Vector
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("ChromatinAssay", "Vector"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = FALSE
  ) {
    query <- granges(x = query)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for ChromatinAssay, ChromatinAssay
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = FALSE
  ) {
    query <- granges(x = query)
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for Vector, Seurat
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("Vector", "Seurat"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay]])
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for Seurat, Vector
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("Seurat", "Vector"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = query)
    query <- granges(x = query[[assay]])
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for Seurat, Seurat
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("Seurat", "Seurat"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    select = c("all", "first", "last", "arbitrary"),
    ignore.strand = FALSE
  ) {
    assay.s <- DefaultAssay(object = subject)
    assay.q <- DefaultAssay(object = query)
    query <- granges(x = query[[assay.q]])
    subject <- granges(x = subject[[assay.s]])
    callGeneric()
  }
)

#' @aliases countOverlaps
#' @describeIn findOverlaps-methods method for Vector, ChromatinAssay
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("Vector", "ChromatinAssay"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for ChromatinAssay, Vector
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("ChromatinAssay", "Vector"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    query <- granges(x = query)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for ChromatinAssay, ChromatinAssay
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("ChromatinAssay", "ChromatinAssay"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    query <- granges(x = query)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for Seurat, Vector
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("Seurat", "Vector"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = query)
    query <- granges(x = query[[assay]])
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for Vector, Seurat
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("Vector", "Seurat"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    assay <- DefaultAssay(object = subject)
    subject <- granges(x = subject[[assay]])
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for Seurat, Seurat
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("Seurat", "Seurat"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {

    assay.s <- DefaultAssay(object = subject)
    assay.q <- DefaultAssay(object = query)
    subject <- granges(x = subject[[assay.s]])
    query <- granges(x = query[[assay.q]])
    callGeneric()
  }
)

## Coverage methods

#' Coverage of a ChromatinAssay object
#'
#' This is the \code{coverage} method for \code{\link{ChromatinAssay}} objects.
#' @param x A \code{\link{ChromatinAssay}} object
#' @param shift How much each range should be shifted before coverage is
#' computed. See \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.
#' @param weight Assigns weight to each range in \code{x}.
#' See \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.
#' @param width Specifies the length of the returned coverage vectors.
#' See \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.
#' @param method See \code{\link[IRanges]{coverage}} in the \pkg{IRanges}
#' package
#'
#' @aliases coverage
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{coverage-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{coverage-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod coverage
#' @describeIn coverage-ChromatinAssay-method method for ChromatinAssay objects
#' @concept coverage
setMethod(
  f = "coverage",
  signature = "ChromatinAssay",
  definition = function(
    x, shift = 0L, width = NULL, weight = 1L,
    method = c("auto", "sort", "hash")
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn coverage-ChromatinAssay-method method for Seurat objects
#' @concept coverage
setMethod(
  f = "coverage",
  signature = "Seurat",
  definition = function(
    x, shift = 0L, width = NULL, weight = 1L,
    method = c("auto", "sort", "hash")
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

## inter-range methods

#' Inter-range transformations for ChromatinAssay objects
#'
#' The \code{range, reduce, gaps, disjoin, isDisjoint, disjointBins} methods
#' are available for \code{\link{ChromatinAssay}} objects.
#'
#' @name inter-range-methods
#' @param x A \code{\link{ChromatinAssay}} object
#' @param ... Additional arguments
#' @param with.revmap See \code{\link[IRanges]{inter-range-methods}} in the
#' \pkg{IRanges} packages
#' @param na.rm Ignored
#'
#' @aliases range range,ChromatinAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{inter-range-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{inter-range-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod range
#' @concept inter_range
setMethod(
  f = "range",
  signature = "ChromatinAssay",
  definition = function(
    x, ..., with.revmap = FALSE, na.rm = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for Seurat objects
#' @concept inter_range
setMethod(
  f = "range",
  signature = "Seurat",
  definition = function(
    x, ..., with.revmap = FALSE, na.rm = FALSE
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @param drop.empty.ranges See \code{?\link{IRanges}{inter-range-methods}}
#' @aliases reduce
#' @describeIn inter-range-methods method for ChromatinAssay objects
#' @exportMethod reduce
#' @concept inter_range
setMethod(
  f = "reduce",
  signature = "ChromatinAssay",
  definition = function(
    x, drop.empty.ranges = FALSE, ...
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for Seurat objects
#' @concept inter_range
setMethod(
  f = "reduce",
  signature = "Seurat",
  definition = function(
    x, drop.empty.ranges = FALSE, ...
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @param start,end See \code{?\link{IRanges}{inter-range-methods}}
#' @aliases gaps
#' @describeIn inter-range-methods method for ChromatinAssay objects
#' @exportMethod gaps
#' @concept inter_range
setMethod(
  f = "gaps",
  signature = "ChromatinAssay",
  definition = function(
    x, start = NA, end = NA
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for Seurat objects
#' @concept inter_range
setMethod(
  f = "gaps",
  signature = "Seurat",
  definition = function(
    x, start = NA, end = NA
  ) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @aliases disjoin
#' @describeIn inter-range-methods method for ChromatinAssay objects
#' @exportMethod disjoin
#' @concept inter_range
setMethod(
  f = "disjoin",
  signature = "ChromatinAssay",
  definition = function(x, ...) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for Seurat objects
#' @concept inter_range
setMethod(
  f = "disjoin",
  signature = "Seurat",
  definition = function(x, ...) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @aliases isDisjoint
#' @describeIn inter-range-methods method for ChromatinAssay objects
#' @exportMethod isDisjoint
#' @concept inter_range
setMethod(
  f = "isDisjoint",
  signature = "ChromatinAssay",
  definition = function(x, ...) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for Seurat objects
#' @concept inter_range
setMethod(
  f = "isDisjoint",
  signature = "Seurat",
  definition = function(x, ...) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

#' @aliases disjointBins
#' @describeIn inter-range-methods method for ChromatinAssay objects
#' @exportMethod disjointBins
#' @concept inter_range
setMethod(
  f = "disjointBins",
  signature = "ChromatinAssay",
  definition = function(x, ...) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for Seurat objects
#' @concept inter_range
setMethod(
  f = "disjointBins",
  signature = "Seurat",
  definition = function(x, ...) {
    assay <- DefaultAssay(object = x)
    x <- granges(x = x[[assay]])
    callGeneric()
  }
)

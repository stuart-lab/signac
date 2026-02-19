#' @include generics.R
#' @importFrom methods callGeneric
#' @importFrom SeuratObject DefaultAssay
#' @importFrom IRanges precede follow nearest distance distanceToNearest
#' findOverlaps countOverlaps coverage
#' reduce disjoin gaps isDisjoint disjointBins
NULL


## Nearest methods

# precede

#' Find the nearest range neighbors for GRangesAssay objects
#'
#' The `precede, follow, nearest, distance, distanceToNearest` methods
#' are available for [GRangesAssay-class] objects.
#'
#' @name nearest-methods
#' @param x A query [GRangesAssay-class] object
#' @param subject The subject [GenomicRanges::GRanges()] or
#' [GRangesAssay-class] object. If missing, `x` is used as the
#' subject.
#' @param select Logic for handling ties. See [GenomicRanges::nearest-methods()]
#' @param ignore.strand Logical argument controlling whether strand information
#' should be ignored.
#' @param ... Additional arguments for methods
#'
#' @aliases precede precede,ANY,GRangesAssay-method
#' @seealso
#'   - [IRanges::nearest()]
#'   - [GenomicRanges::nearest()]
#'   - [GRangesAssay-class]
#' @exportMethod precede
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("ANY", "GRangesAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, ANY
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("GRangesAssay", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, GRangesAssay
#' @concept nearest
setMethod(
  f = "precede",
  signature = c("GRangesAssay", "GRangesAssay"),
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
#' @describeIn nearest-methods method for ANY, GRangesAssay
#' @exportMethod follow
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("ANY", "GRangesAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, ANY
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("GRangesAssay", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, GRangesAssay
#' @concept nearest
setMethod(
  f = "follow",
  signature = c("GRangesAssay", "GRangesAssay"),
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
#' @describeIn nearest-methods method for ANY, GRangesAssay
#' @exportMethod nearest
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("ANY", "GRangesAssay"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, ANY
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("GRangesAssay", "ANY"),
  definition = function(
    x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, GRangesAssay
#' @concept nearest
setMethod(
  f = "nearest",
  signature = c("GRangesAssay", "GRangesAssay"),
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
    x <- granges(x = x[[assay.x]])
    callGeneric()
  }
)

#' @param y For the `distance` method, a
#' [GenomicRanges::GRanges()] object or a [GRangesAssay-class]
#' object
#'
#' @aliases distance
#' @exportMethod distance
#' @describeIn nearest-methods method for ANY, GRangesAssay
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("ANY", "GRangesAssay"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    y <- granges(x = y)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, ANY
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("GRangesAssay", "ANY"),
  definition = function(
    x, y, ignore.strand = FALSE, ...
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, GRangesAssay
#' @concept nearest
setMethod(
  f = "distance",
  signature = c("GRangesAssay", "GRangesAssay"),
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
#' @describeIn nearest-methods method for ANY, GRangesAssay
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("ANY", "GRangesAssay"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, ANY
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("GRangesAssay", "ANY"),
  definition = function(
    x, subject, ignore.strand = FALSE, ...
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn nearest-methods method for GRangesAssay, GRangesAssay
#' @concept nearest
setMethod(
  f = "distanceToNearest",
  signature = c("GRangesAssay", "GRangesAssay"),
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

#' Find overlapping ranges for GRangesAssay objects
#'
#' The `findOverlaps, countOverlaps` methods are available for
#' [GRangesAssay-class] objects. This allows finding overlaps between
#' genomic ranges and the ranges stored in the [GRangesAssay-class]
#'
#' If a [GRangesAssay-class] assay is set as the default assay in a
#' [SeuratObject::Seurat()] object, you can also call `findOverlaps`
#' directly on the Seurat object.
#'
#' @param query,subject A [GRangesAssay-class] object
#' @param maxgap,minoverlap,type,select,ignore.strand See
#' `?[findOverlaps][GenomicRanges::findOverlaps]` in the [GenomicRanges] and
#' [IRanges] packages.
#' @return See [GenomicRanges::findOverlaps()]
#'
#' @name findOverlaps-methods
#' @aliases findOverlaps findOverlaps,Vector,GRangesAssay-method
#' @seealso
#'   - [IRanges::findOverlaps-methods]
#'   - [GenomicRanges::findOverlaps-methods]
#'   - [GRangesAssay-class]
#'
#' @exportMethod findOverlaps
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("Vector", "GRangesAssay"),
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

#' @describeIn findOverlaps-methods method for GRangesAssay, Vector
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("GRangesAssay", "Vector"),
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

#' @describeIn findOverlaps-methods method for GRangesAssay, GRangesAssay
#' @concept overlaps
setMethod(
  f = "findOverlaps",
  signature = c("GRangesAssay", "GRangesAssay"),
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
#' @describeIn findOverlaps-methods method for Vector, GRangesAssay
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("Vector", "GRangesAssay"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    subject <- granges(x = subject)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for GRangesAssay, Vector
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("GRangesAssay", "Vector"),
  definition = function(
    query, subject, maxgap = -1L, minoverlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore.strand = FALSE
  ) {
    query <- granges(x = query)
    callGeneric()
  }
)

#' @describeIn findOverlaps-methods method for GRangesAssay, GRangesAssay
#' @concept overlaps
setMethod(
  f = "countOverlaps",
  signature = c("GRangesAssay", "GRangesAssay"),
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

#' Coverage of a [GRangesAssay-class] object
#'
#' This is the `coverage` method for [GRangesAssay-class] objects.
#' @param x A [GRangesAssay-class] object
#' @param shift How much each range should be shifted before coverage is
#' computed. See [IRanges::coverage()] in the [IRanges] package.
#' @param weight Assigns weight to each range in `x`.
#' See [IRanges::coverage()] in the [IRanges] package.
#' @param width Specifies the length of the returned coverage vectors.
#' See [IRanges::coverage()] in the [IRanges] package.
#' @param method See [IRanges::coverage()] in the [IRanges]
#' package
#'
#' @aliases coverage
#' @seealso
#'   - [IRanges::coverage()]
#'   - [GenomicRanges::coverage()]
#'   - [GRangesAssay-class]
#' @exportMethod coverage
#' @describeIn coverage-GRangesAssay-method method for GRangesAssay objects
#' @concept coverage
setMethod(
  f = "coverage",
  signature = "GRangesAssay",
  definition = function(
    x, shift = 0L, width = NULL, weight = 1L,
    method = c("auto", "sort", "hash")
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn coverage-GRangesAssay-method method for Seurat objects
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

#' Inter-range transformations for [GRangesAssay-class] objects
#'
#' The `range, reduce, gaps, disjoin, isDisjoint, disjointBins` methods
#' are available for [GRangesAssay-class] objects.
#'
#' @name inter-range-methods
#' @param x A [GRangesAssay-class] object
#' @param ... Additional arguments
#' @param with.revmap See [IRanges::inter-range-methods()] in the
#' \pkg{IRanges} packages
#' @param na.rm Ignored
#'
#' @aliases range range,GRangesAssay-method
#' @seealso
#'   - [IRanges::inter-range-methods]
#'   - [GenomicRanges::inter-range-methods]
#'   - [GRangesAssay-class]
#' @exportMethod range
#' @concept inter_range
setMethod(
  f = "range",
  signature = "GRangesAssay",
  definition = function(
    x, ..., with.revmap = FALSE, na.rm = FALSE
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for [SeuratObject::Seurat-class]
#' objects
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

#' @param drop.empty.ranges See
#' `?[inter-range-methods][IRanges::inter-range-methods]`
#' @aliases reduce
#' @describeIn inter-range-methods method for [GRangesAssay-class] objects
#' @exportMethod reduce
#' @concept inter_range
setMethod(
  f = "reduce",
  signature = "GRangesAssay",
  definition = function(
    x, drop.empty.ranges = FALSE, ...
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for [SeuratObject::Seurat-class]
#' objects
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

#' @param start,end See `?[inter-range-methods][IRanges::inter-range-methods]`
#' @aliases gaps
#' @describeIn inter-range-methods method for [GRangesAssay-class] objects
#' @exportMethod gaps
#' @concept inter_range
setMethod(
  f = "gaps",
  signature = "GRangesAssay",
  definition = function(
    x, start = NA, end = NA
  ) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for [SeuratObject::Seurat-class]
#' objects
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
#' @describeIn inter-range-methods method for [GRangesAssay-class] objects
#' @exportMethod disjoin
#' @concept inter_range
setMethod(
  f = "disjoin",
  signature = "GRangesAssay",
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
#' @describeIn inter-range-methods method for [GRangesAssay-class] objects
#' @exportMethod isDisjoint
#' @concept inter_range
setMethod(
  f = "isDisjoint",
  signature = "GRangesAssay",
  definition = function(x, ...) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for [SeuratObject::Seurat-class]
#' objects
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
#' @describeIn inter-range-methods method for [GRangesAssay-class] objects
#' @exportMethod disjointBins
#' @concept inter_range
setMethod(
  f = "disjointBins",
  signature = "GRangesAssay",
  definition = function(x, ...) {
    x <- granges(x = x)
    callGeneric()
  }
)

#' @describeIn inter-range-methods method for [SeuratObject::Seurat-class]
#' objects
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

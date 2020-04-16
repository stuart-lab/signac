#' @include generics.R
#' @importFrom methods callGeneric
#' @importFrom Seurat DefaultAssay
#' @importFrom IRanges precede follow nearest distance distanceToNearest
#' findOverlaps countOverlaps coverage

setOldClass(Classes = "ChromatinAssay")

## Nearest methods

# precede

#' Find the nearest range neighbors for ChromatinAssay objects
#'
#' The \code{precede, follow, nearest, distance, distanceToNearest} methods
#' are available for \code{\link{ChromatinAssay}} objects.
#'
#' @name nearest-methods
#' @aliases precede
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{nearest-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{nearest-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod precede
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

#' @describeIn nearest-methods method for ANY, ChromatinAssay
#' @exportMethod follow
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

# nearest
#' @describeIn nearest-methods method for ANY, ChromatinAssay
#' @exportMethod nearest
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

# distance

#' @exportMethod distance
#' @describeIn nearest-methods method for ANY, ChromatinAssay
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

# distanceToNearest
#' @exportMethod distanceToNearest
#' @describeIn nearest-methods method for ANY, ChromatinAssay
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
#' @name findOverlaps-methods
#' @aliases findOverlaps
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{findOverlaps-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{findOverlaps-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#'
#' @exportMethod findOverlaps
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

# countOverlaps

#' @describeIn findOverlaps-methods method for Vector, ChromatinAssay
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
#' @aliases coverage
#' @seealso
#' \itemize{
#'   \item{\link[IRanges]{coverage-methods} in the \pkg{IRanges} package.}
#'   \item{\link[GenomicRanges]{coverage-methods} in the \pkg{GenomicRanges}
#'   package}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod coverage
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

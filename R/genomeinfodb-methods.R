#' @include generics.R
#' @importFrom methods callGeneric
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomeInfoDb
#' seqinfo seqinfo<-
#' seqnames seqnames<-
#' seqlevels seqlevels<-
#' sortSeqlevels
#' seqlengths seqlengths<-
#' isCircular isCircular<-
#' genome genome<-
NULL

setOldClass(Classes = "ChromatinAssay")

#' Access and modify sequence information for ChromatinAssay objects
#'
#' Methods for accessing and modifying
#' \code{\link[GenomeInfoDb]{Seqinfo}} object information stored in a
#' \code{\link{ChromatinAssay}} object.
#'
#' @name seqinfo-methods
#' @aliases seqinfo
#' @seealso
#' \itemize{
#'   \item{\link[GenomeInfoDb]{seqinfo} in the \pkg{GenomeInfoDb} package.}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod seqinfo
setMethod(
  f = "seqinfo",
  signature = "ChromatinAssay",
  definition = function(x) {
    slot(object = x, name = "seqinfo")
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqinfo<-
setMethod(
  f = "seqinfo<-",
  signature = "ChromatinAssay",
  definition = function(
    x,
    value
  ) {
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = value)
    x
  }
)

#' @aliases seqlevels
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod seqlevels
setMethod(
  f = "seqlevels",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqlevels<-
setMethod(
  f = "seqlevels<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @aliases seqnames
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod seqnames
setMethod(
  f = "seqnames",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqnames<-
setMethod(
  f = "seqnames<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @aliases seqlengths
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod seqlengths
setMethod(
  f = "seqlengths",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqlengths<-
setMethod(
  f = "seqlengths<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @aliases genome
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod genome
setMethod(
  f = "genome",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod genome<-
setMethod(
  f = "genome<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @aliases isCircular
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod isCircular
setMethod(
  f = "isCircular",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod isCircular<-
setMethod(
  f = "isCircular<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
setMethod(
  f = "seqinfo",
  signature = "Seurat",
  definition = function(x) {
    assay <- DefaultAssay(object = x)
    slot(object = x[[assay]], name = "seqinfo")
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
setMethod(
  f = "seqinfo<-",
  signature = "Seurat",
  definition = function(
    x,
    value
  ) {
    assay <- DefaultAssay(object = x)
    seqinfo(x = x[[assay]]) <- value
    x
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
setMethod(
  f = "seqlevels",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
setMethod(
  f = "seqlevels<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
setMethod(
  f = "seqnames",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
setMethod(
  f = "seqnames<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
setMethod(
  f = "seqlengths",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
setMethod(
  f = "seqlengths<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
setMethod(
  f = "genome",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
setMethod(
  f = "genome<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
setMethod(
  f = "isCircular",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
setMethod(
  f = "isCircular<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

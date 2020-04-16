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

setOldClass(Classes = "ChromatinAssay")

#' Access and modify sequence information for ChromatinAssay objects
#'
#' Generic functions for accessing and modifying
#' \code{\link[GenomeInfoDb]{Seqinfo}} information stored in a
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

#' @describeIn seqinfo-methods method for ChromatinAssay objects
#' @exportMethod seqlevels
setMethod(
  f = "seqlevels",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for ChromatinAssay objects
#' @exportMethod seqnames
setMethod(
  f = "seqnames",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for ChromatinAssay objects
#' @exportMethod seqlengths
setMethod(
  f = "seqlengths",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)
#' @describeIn seqinfo-methods method for ChromatinAssay objects
#' @exportMethod genome
setMethod(
  f = "genome",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for ChromatinAssay objects
#' @exportMethod isCircular
setMethod(
  f = "isCircular",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for Seurat objects objects
setMethod(
  f = "seqinfo",
  signature = "Seurat",
  definition = function(x) {
    assay <- DefaultAssay(object = x)
    slot(object = x[[assay]], name = "seqinfo")
  }
)

#' @describeIn seqinfo-methods method for Seurat objects
#' @exportMethod seqlevels
setMethod(
  f = "seqlevels",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for Seurat objects
#' @exportMethod seqnames
setMethod(
  f = "seqnames",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for Seurat objects
#' @exportMethod seqlengths
setMethod(
  f = "seqlengths",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)
#' @describeIn seqinfo-methods method for Seurat objects
#' @exportMethod genome
setMethod(
  f = "genome",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

#' @describeIn seqinfo-methods method for Seurat objects
#' @exportMethod isCircular
setMethod(
  f = "isCircular",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    callGeneric()
  }
)

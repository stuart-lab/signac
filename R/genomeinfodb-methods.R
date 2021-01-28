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

# TODO add keepSeqlevels, dropSeqlevels, renameSeqlevels, etc

#' Access and modify sequence information for ChromatinAssay objects
#'
#' Methods for accessing and modifying
#' \code{\link[GenomeInfoDb]{Seqinfo}} object information stored in a
#' \code{\link{ChromatinAssay}} object.
#'
#' @name seqinfo-methods
#' @param x A \code{\link{ChromatinAssay}} object
#'
#' @aliases seqinfo seqinfo,ChromatinAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[GenomeInfoDb]{seqinfo} in the \pkg{GenomeInfoDb} package.}
#'   \item{\link{ChromatinAssay-class}}
#'  }
#' @exportMethod seqinfo
#' @concept seqinfo
setMethod(
  f = "seqinfo",
  signature = "ChromatinAssay",
  definition = function(x) {
    slot(object = x, name = "seqinfo")
  }
)

#' @param value A \code{\link[GenomeInfoDb]{Seqinfo}} object or name of a UCSC
#' genome to store in the \code{\link{ChromatinAssay}}
#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqinfo<-
#' @concept seqinfo
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
#' @concept seqinfo
setMethod(
  f = "seqlevels",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqlevels<-
#' @concept seqinfo
setMethod(
  f = "seqlevels<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    sinfo <- seqinfo(x = x)
    seqlevels(x = sinfo) <- value
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = sinfo)
    return(x)
  }
)

#' @aliases seqnames
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod seqnames
#' @concept seqinfo
setMethod(
  f = "seqnames",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqnames<-
#' @concept seqinfo
setMethod(
  f = "seqnames<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    sinfo <- seqinfo(x = x)
    seqnames(x = sinfo) <- value
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = sinfo)
    return(x)
  }
)

#' @aliases seqlengths
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod seqlengths
#' @concept seqinfo
setMethod(
  f = "seqlengths",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      if (is.null(x = x)) {
        return(NULL)
      } else {
        callGeneric()
      }
    }
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod seqlengths<-
#' @concept seqinfo
setMethod(
  f = "seqlengths<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    sinfo <- seqinfo(x = x)
    seqlengths(x = sinfo) <- value
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = sinfo)
    return(x)
  }
)

#' @aliases genome
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod genome
#' @concept seqinfo
setMethod(
  f = "genome",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @importFrom Seurat SetAssayData
#' @exportMethod genome<-
#' @concept seqinfo
setMethod(
  f = "genome<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = value)
    return(x)
  }
)

#' @aliases isCircular
#' @describeIn seqinfo-methods get method for ChromatinAssay objects
#' @exportMethod isCircular
#' @concept seqinfo
setMethod(
  f = "isCircular",
  signature = "ChromatinAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for ChromatinAssay objects
#' @exportMethod isCircular<-
#' @concept seqinfo
setMethod(
  f = "isCircular<-",
  signature = "ChromatinAssay",
  definition = function(
    x, value
  ) {
    sinfo <- seqinfo(x = x)
    isCircular(x = sinfo) <- value
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = sinfo)
    return(x)
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqinfo",
  signature = "Seurat",
  definition = function(x) {
    assay <- DefaultAssay(object = x)
    slot(object = x[[assay]], name = "seqinfo")
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
#' @concept seqinfo
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
#' @concept seqinfo
setMethod(
  f = "seqlevels",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqlevels<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    assay <- DefaultAssay(object = x)
    seqlevels(x = x[[assay]]) <- value
    x
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqnames",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqnames<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    assay <- DefaultAssay(object = x)
    seqnames(x = x[[assay]]) <- value
    x
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqlengths",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqlengths<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    assay <- DefaultAssay(object = x)
    seqlengths(x = x[[assay]]) <- value
    x
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "genome",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
#' @importFrom Seurat SetAssayData
#' @concept seqinfo
setMethod(
  f = "genome<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    x <- SetAssayData(object = x, slot = "seqinfo", new.data = value)
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "isCircular",
  signature = "Seurat",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods set method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "isCircular<-",
  signature = "Seurat",
  definition = function(
    x, value
  ) {
    assay <- DefaultAssay(object = x)
    isCircular(x = x[[assay]]) <- value
    x
  }
)

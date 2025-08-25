#' @include generics.R
#' @importFrom methods callGeneric
#' @importFrom SeuratObject DefaultAssay
#' @importFrom GenomeInfoDb
#' seqinfo
#' seqnames
#' seqlevels
#' sortSeqlevels
#' seqlengths
#' isCircular
#' genome
NULL


# TODO add keepSeqlevels, dropSeqlevels, renameSeqlevels, etc

#' Access and modify sequence information for GRangesAssay objects
#'
#' Methods for accessing and modifying
#' \code{\link[GenomeInfoDb]{Seqinfo}} object information stored in a
#' \code{\link{GRangesAssay}} object.
#'
#' @name seqinfo-methods
#' @param x A \code{\link{GRangesAssay}} object
#'
#' @aliases seqinfo seqinfo,GRangesAssay-method
#' @seealso
#' \itemize{
#'   \item{\link[GenomeInfoDb]{seqinfo} in the \pkg{GenomeInfoDb} package.}
#'   \item{\link{GRangesAssay-class}}
#'  }
#' @exportMethod seqinfo
#' @concept seqinfo
setMethod(
  f = "seqinfo",
  signature = "GRangesAssay",
  definition = function(x) {
    seqinfo(x = granges(x = x))
  }
)

#' @aliases seqlevels
#' @describeIn seqinfo-methods get method for GRangesAssay objects
#' @exportMethod seqlevels
#' @concept seqinfo
setMethod(
  f = "seqlevels",
  signature = "GRangesAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @aliases seqnames
#' @describeIn seqinfo-methods get method for GRangesAssay objects
#' @exportMethod seqnames
#' @concept seqinfo
setMethod(
  f = "seqnames",
  signature = "GRangesAssay",
  definition = function(x) {
    x <- granges(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @aliases seqlengths
#' @describeIn seqinfo-methods get method for GRangesAssay objects
#' @exportMethod seqlengths
#' @concept seqinfo
setMethod(
  f = "seqlengths",
  signature = "GRangesAssay",
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

#' @aliases genome
#' @describeIn seqinfo-methods get method for GRangesAssay objects
#' @exportMethod genome
#' @concept seqinfo
setMethod(
  f = "genome",
  signature = "GRangesAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @aliases isCircular
#' @describeIn seqinfo-methods get method for GRangesAssay objects
#' @exportMethod isCircular
#' @concept seqinfo
setMethod(
  f = "isCircular",
  signature = "GRangesAssay",
  definition = function(x) {
    x <- seqinfo(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
  }
)

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqinfo",
  signature = "Seurat",
  definition = function(x) {
    assay <- DefaultAssay(object = x)
    seqinfo(x = x[[assay]])
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

#' @describeIn seqinfo-methods get method for Seurat objects
#' @concept seqinfo
setMethod(
  f = "seqnames",
  signature = "Seurat",
  definition = function(x) {
    x <- granges(x = x)
    if (is.null(x = x)) {
      return(NULL)
    } else {
      callGeneric()
    }
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

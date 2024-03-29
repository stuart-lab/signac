% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{nearest-methods}
\alias{nearest-methods}
\alias{precede}
\alias{precede,ANY,ChromatinAssay-method}
\alias{precede,ChromatinAssay,ANY-method}
\alias{precede,ChromatinAssay,ChromatinAssay-method}
\alias{precede,ANY,Seurat-method}
\alias{precede,Seurat,ANY-method}
\alias{precede,Seurat,Seurat-method}
\alias{follow,ANY,ChromatinAssay-method}
\alias{follow}
\alias{follow,ChromatinAssay,ANY-method}
\alias{follow,ChromatinAssay,ChromatinAssay-method}
\alias{follow,ANY,Seurat-method}
\alias{follow,Seurat,ANY-method}
\alias{follow,Seurat,Seurat-method}
\alias{nearest,ANY,ChromatinAssay-method}
\alias{nearest}
\alias{nearest,ChromatinAssay,ANY-method}
\alias{nearest,ChromatinAssay,ChromatinAssay-method}
\alias{nearest,ANY,Seurat-method}
\alias{nearest,Seurat,ANY-method}
\alias{nearest,Seurat,Seurat-method}
\alias{distance,ANY,ChromatinAssay-method}
\alias{distance}
\alias{distance,ChromatinAssay,ANY-method}
\alias{distance,ChromatinAssay,ChromatinAssay-method}
\alias{distance,ANY,Seurat-method}
\alias{distance,Seurat,ANY-method}
\alias{distance,Seurat,Seurat-method}
\alias{distanceToNearest,ANY,ChromatinAssay-method}
\alias{distanceToNearest}
\alias{distanceToNearest,ChromatinAssay,ANY-method}
\alias{distanceToNearest,ChromatinAssay,ChromatinAssay-method}
\alias{distanceToNearest,ANY,Seurat-method}
\alias{distanceToNearest,Seurat,ANY-method}
\alias{distanceToNearest,Seurat,Seurat-method}
\title{Find the nearest range neighbors for ChromatinAssay objects}
\usage{
\S4method{precede}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ANY,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{Seurat,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{Seurat,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ANY,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{Seurat,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{Seurat,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ANY,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{Seurat,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{Seurat,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{distance}{ANY,ChromatinAssay}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ChromatinAssay,ANY}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ChromatinAssay,ChromatinAssay}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ANY,Seurat}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{Seurat,ANY}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{Seurat,Seurat}(x, y, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ANY,ChromatinAssay}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ChromatinAssay,ANY}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ChromatinAssay,ChromatinAssay}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ANY,Seurat}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{Seurat,ANY}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{Seurat,Seurat}(x, subject, ignore.strand = FALSE, ...)
}
\arguments{
\item{x}{A query \code{\link{ChromatinAssay}} object}

\item{subject}{The subject \code{\link[GenomicRanges]{GRanges}} or
\code{\link{ChromatinAssay}} object. If missing, \code{x} is used as the
subject.}

\item{select}{Logic for handling ties.
See \code{\link[GenomicRanges]{nearest-methods}} in the \pkg{GenomicRanges}
package.}

\item{ignore.strand}{Logical argument controlling whether strand information
should be ignored.}

\item{y}{For the \code{distance} method, a
\code{\link[GenomicRanges]{GRanges}} object or a \code{\link{ChromatinAssay}}
object}

\item{...}{Additional arguments for methods}
}
\description{
The \code{precede, follow, nearest, distance, distanceToNearest} methods
are available for \code{\link{ChromatinAssay}} objects.
}
\section{Functions}{
\itemize{
\item \code{precede(x = ChromatinAssay, subject = ANY)}: method for ChromatinAssay, ANY

\item \code{precede(x = ChromatinAssay, subject = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{precede(x = ANY, subject = Seurat)}: method for ANY, Seurat

\item \code{precede(x = Seurat, subject = ANY)}: method for Seurat, ANY

\item \code{precede(x = Seurat, subject = Seurat)}: method for Seurat, Seurat

\item \code{follow(x = ANY, subject = ChromatinAssay)}: method for ANY, ChromatinAssay

\item \code{follow(x = ChromatinAssay, subject = ANY)}: method for ChromatinAssay, ANY

\item \code{follow(x = ChromatinAssay, subject = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{follow(x = ANY, subject = Seurat)}: method for ANY, Seurat

\item \code{follow(x = Seurat, subject = ANY)}: method for Seurat, ANY

\item \code{follow(x = Seurat, subject = Seurat)}: method for Seurat, Seurat

\item \code{nearest(x = ANY, subject = ChromatinAssay)}: method for ANY, ChromatinAssay

\item \code{nearest(x = ChromatinAssay, subject = ANY)}: method for ChromatinAssay, ANY

\item \code{nearest(x = ChromatinAssay, subject = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{nearest(x = ANY, subject = Seurat)}: method for ANY, Seurat

\item \code{nearest(x = Seurat, subject = ANY)}: method for Seurat, ANY

\item \code{nearest(x = Seurat, subject = Seurat)}: method for Seurat, Seurat

\item \code{distance(x = ANY, y = ChromatinAssay)}: method for ANY, ChromatinAssay

\item \code{distance(x = ChromatinAssay, y = ANY)}: method for ChromatinAssay, ANY

\item \code{distance(x = ChromatinAssay, y = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{distance(x = ANY, y = Seurat)}: method for ANY, Seurat

\item \code{distance(x = Seurat, y = ANY)}: method for Seurat, ANY

\item \code{distance(x = Seurat, y = Seurat)}: method for Seurat, Seurat

\item \code{distanceToNearest(x = ANY, subject = ChromatinAssay)}: method for ANY, ChromatinAssay

\item \code{distanceToNearest(x = ChromatinAssay, subject = ANY)}: method for ChromatinAssay, ANY

\item \code{distanceToNearest(x = ChromatinAssay, subject = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{distanceToNearest(x = ANY, subject = Seurat)}: method for ANY, Seurat

\item \code{distanceToNearest(x = Seurat, subject = ANY)}: method for Seurat, ANY

\item \code{distanceToNearest(x = Seurat, subject = Seurat)}: method for Seurat, Seurat

}}
\seealso{
\itemize{
  \item{\link[IRanges]{nearest-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{nearest-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{nearest}

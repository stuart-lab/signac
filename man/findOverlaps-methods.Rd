% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{findOverlaps-methods}
\alias{findOverlaps-methods}
\alias{findOverlaps}
\alias{findOverlaps,Vector,ChromatinAssay-method}
\alias{findOverlaps,ChromatinAssay,Vector-method}
\alias{findOverlaps,ChromatinAssay,ChromatinAssay-method}
\alias{findOverlaps,Vector,Seurat-method}
\alias{findOverlaps,Seurat,Vector-method}
\alias{findOverlaps,Seurat,Seurat-method}
\alias{countOverlaps,Vector,ChromatinAssay-method}
\alias{countOverlaps}
\alias{countOverlaps,ChromatinAssay,Vector-method}
\alias{countOverlaps,ChromatinAssay,ChromatinAssay-method}
\alias{countOverlaps,Seurat,Vector-method}
\alias{countOverlaps,Vector,Seurat-method}
\alias{countOverlaps,Seurat,Seurat-method}
\title{Find overlapping ranges for ChromatinAssay objects}
\usage{
\S4method{findOverlaps}{Vector,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{ChromatinAssay,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{ChromatinAssay,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{Vector,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{Seurat,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{Seurat,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Vector,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{ChromatinAssay,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{ChromatinAssay,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Seurat,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Vector,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Seurat,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)
}
\arguments{
\item{query, subject}{A \code{\link{ChromatinAssay}} object}

\item{maxgap, minoverlap, type, select, ignore.strand}{See
\code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} and
\pkg{IRanges} packages.}
}
\value{
See \code{\link[GenomicRanges]{findOverlaps}}
}
\description{
The \code{findOverlaps, countOverlaps} methods are available for
\code{\link{ChromatinAssay}} objects. This allows finding overlaps between
genomic ranges and the ranges stored in the ChromatinAssay.
}
\details{
If a ChromatinAssay is set as the default assay in a
\code{\link[SeuratObject]{Seurat}} object, you can also call \code{findOverlaps}
directly on the Seurat object.
}
\section{Functions}{
\itemize{
\item \code{findOverlaps(query = ChromatinAssay, subject = Vector)}: method for ChromatinAssay, Vector

\item \code{findOverlaps(query = ChromatinAssay, subject = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{findOverlaps(query = Vector, subject = Seurat)}: method for Vector, Seurat

\item \code{findOverlaps(query = Seurat, subject = Vector)}: method for Seurat, Vector

\item \code{findOverlaps(query = Seurat, subject = Seurat)}: method for Seurat, Seurat

\item \code{countOverlaps(query = Vector, subject = ChromatinAssay)}: method for Vector, ChromatinAssay

\item \code{countOverlaps(query = ChromatinAssay, subject = Vector)}: method for ChromatinAssay, Vector

\item \code{countOverlaps(query = ChromatinAssay, subject = ChromatinAssay)}: method for ChromatinAssay, ChromatinAssay

\item \code{countOverlaps(query = Seurat, subject = Vector)}: method for Seurat, Vector

\item \code{countOverlaps(query = Vector, subject = Seurat)}: method for Vector, Seurat

\item \code{countOverlaps(query = Seurat, subject = Seurat)}: method for Seurat, Seurat

}}
\seealso{
\itemize{
  \item{\link[IRanges]{findOverlaps-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{findOverlaps-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{overlaps}

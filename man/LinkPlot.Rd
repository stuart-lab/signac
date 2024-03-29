% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{LinkPlot}
\alias{LinkPlot}
\title{Plot linked genomic elements}
\usage{
LinkPlot(
  object,
  region,
  assay = NULL,
  min.cutoff = 0,
  sep = c("-", "-"),
  extend.upstream = 0,
  extend.downstream = 0,
  scale.linewidth = FALSE
)
}
\arguments{
\item{object}{A \code{\link[SeuratObject]{Seurat}} object}

\item{region}{A genomic region to plot}

\item{assay}{Name of assay to use. If NULL, use the default assay.}

\item{min.cutoff}{Minimum absolute score for link to be plotted.}

\item{sep}{Separators to use for strings encoding genomic coordinates. First
element is used to separate the chromosome from the coordinates, second
element is used to separate the start from end coordinate.}

\item{extend.upstream}{Number of bases to extend the region upstream.}

\item{extend.downstream}{Number of bases to extend the region downstream.}

\item{scale.linewidth}{Scale thickness of the line according to link score.}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Display links between pairs of genomic elements within a given region of the
genome.
}
\concept{links}
\concept{visualization}

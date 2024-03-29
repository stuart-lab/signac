% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/footprinting.R
\name{Footprint}
\alias{Footprint}
\alias{Footprint.ChromatinAssay}
\alias{Footprint.Seurat}
\title{Transcription factor footprinting analysis}
\usage{
Footprint(object, ...)

\method{Footprint}{ChromatinAssay}(
  object,
  genome,
  motif.name = NULL,
  key = motif.name,
  regions = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  compute.expected = TRUE,
  in.peaks = FALSE,
  verbose = TRUE,
  ...
)

\method{Footprint}{Seurat}(
  object,
  genome,
  regions = NULL,
  motif.name = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  in.peaks = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat or ChromatinAssay object}

\item{...}{Arguments passed to other methods}

\item{genome}{A \code{BSgenome} object or any other object supported by
\code{getSeq}. Do \code{showMethods("getSeq")} to get the list of all
supported object types.}

\item{motif.name}{Name of a motif stored in the assay to footprint. If not
supplied, must supply a set of regions.}

\item{key}{Key to store positional enrichment information under.}

\item{regions}{A set of genomic ranges containing the motif instances. These
should all be the same width.}

\item{assay}{Name of assay to use}

\item{upstream}{Number of bases to extend upstream}

\item{downstream}{Number of bases to extend downstream}

\item{compute.expected}{Find the expected number of insertions at each
position given the local DNA sequence context and the insertion bias of Tn5}

\item{in.peaks}{Restrict motifs to those that fall in peaks}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Compute the normalized observed/expected Tn5 insertion frequency
for each position surrounding a set of motif instances.
}
\concept{footprinting}

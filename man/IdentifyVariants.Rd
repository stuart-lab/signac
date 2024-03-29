% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/mito.R
\name{IdentifyVariants}
\alias{IdentifyVariants}
\alias{IdentifyVariants.default}
\alias{IdentifyVariants.Assay}
\alias{IdentifyVariants.StdAssay}
\alias{IdentifyVariants.Seurat}
\title{Identify mitochondrial variants}
\usage{
IdentifyVariants(object, ...)

\method{IdentifyVariants}{default}(
  object,
  refallele,
  stabilize_variance = TRUE,
  low_coverage_threshold = 10,
  verbose = TRUE,
  ...
)

\method{IdentifyVariants}{Assay}(object, refallele, ...)

\method{IdentifyVariants}{StdAssay}(object, refallele, ...)

\method{IdentifyVariants}{Seurat}(object, refallele, assay = NULL, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{refallele}{A dataframe containing reference alleles for the
mitochondrial genome.}

\item{stabilize_variance}{Stabilize variance}

\item{low_coverage_threshold}{Low coverage threshold}

\item{verbose}{Display messages}

\item{assay}{Name of assay to use. If NULL, use the default assay.}
}
\value{
Returns a dataframe
}
\description{
Identify mitochondrial variants present in single cells.
}
\examples{
\dontrun{
data.dir <- "path/to/data/directory"
mgatk <- ReadMGATK(dir = data.dir)
variant.df <- IdentifyVariants(
  object = mgatk$counts,
  refallele = mgatk$refallele
)
}
}
\concept{mito}

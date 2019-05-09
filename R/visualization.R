
#' Plot coverage within given region for groups of cells
#'
#' @param object A Seurat object
#' @param coords A set of genomic coordinates to show
#' @param annotation A set of genomic feature annotations. Must be the name of a set of features
#' stored using the \code{SetAnnotations} function
#' @param assay Which assay to plot
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{SetFragments} function
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param cols Vector of colors, where each color corresponds to an identity class. By default, use the ggplot2 colors.
#'
#' @export
#'
CoveragePlot <- function(
  object,
  coords,
  annotation = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  cells = NULL,
  cols = NULL
) {
  # TODO
  return()
}

#' Plot coverage pileup centered on a given genomic feature
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param feature Which genomic feature to center on. Options are: 'TSS', 'TTS'
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{SetFragments} function
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param cols Vector of colors, where each color corresponds to an identity class. By default, use the ggplot2 colors.
#'
#' @export
#'
PileupPlot <- function(
  object,
  assay = NULL,
  feature = 'TSS',
  fragment.path = NULL,
  group.by = NULL,
  cells = NULL,
  cols = NULL
) {
  # TODO
  return()
}

#' Plot fragment length periodicity
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{SetFragments} function
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param cols Vector of colors, where each color corresponds to an identity class. By default, use the ggplot2 colors.
#'
#' @export
#'
PeriodPlot <- function(
  object,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  cells = NULL,
  cols = NULL
) {
  # TODO
  return()
}

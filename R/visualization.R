
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
#' @param range Genomic range to use. Default is fist megabase of chromosome 1.
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#'
#' @importFrom seqminer tabix.read.table
#' @importFrom ggplot2 ggplot geom_histogram theme_bw aes facet_wrap xlim
#' @return Returns a ggplot2 object
#' @export
#'
PeriodPlot <- function(
  object,
  assay = NULL,
  fragment.path = NULL,
  range = 'chr1:1-1000000',
  group.by = NULL,
  cells = NULL
) {
  assay <- assay %||% DefaultAssay(object)
  if (is.null(group.by)) {
    group.by <- Idents(object)
  } else {
    meta.data <- object[[]]
    group.by <- meta.data[[group.by]]
    names(group.by) <- rownames(meta.data)
  }
  if (is.null(fragment.path)) {
    tools <- slot(object = object, name = 'tools')
    if ('fragments' %in% names(tools)) {
      if (assay %in% names(tools$fragments)) {
        fragment.path <- tools$fragments[[assay]]
      } else {
        stop('Fragment file not supplied for the requested assay')
      }
    } else {
      stop('Fragment file not set. Run SetFragments to set the fragment file path.')
    }
  } else if (!(all(file.exists(fragment.path, paste0(fragment.path, '.tbi'))))) {
      stop('Requested file does not exist or is not indexed')
  }
  reads <- tabix.read.table(tabixFile = fragment.path, tabixRange = range)
  colnames(reads) <- c('chrom', 'start', 'stop', 'cell', 'reads')
  reads <- reads[reads$cell %in% names(group.by), ]
  if (!is.null(cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    stop('No cells present in the requested region')
  }
  reads$length <- reads$stop - reads$start
  reads$group <- group.by[reads$cell]
  if (length(unique(reads$group)) == 1) {
    p <- ggplot(data = reads, aes(length)) +
      geom_histogram(bins = 200) +
      xlim(c(0, 800)) +
      theme_bw()
  } else {
    p <- ggplot(data = reads, aes(length, fill = group)) +
      geom_histogram(bins = 200) +
      facet_wrap(~group, scales = 'free_y') +
      xlim(c(0, 800)) +
      theme_bw()
  }
  return(p)
}

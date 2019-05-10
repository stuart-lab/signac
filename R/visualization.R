
#' Plot coverage within given region for groups of cells
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show
#' @param annotation A set of genomic feature annotations. Must be the name of a set of features
#' stored using the \code{SetAnnotations} function
#' @param assay Which assay to plot
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{SetFragments} function
#' @param cells Which cells to plot. Default all cells
#' @param window Smoothing window size
#' @param downsample Fraction of positions to retain in the plot. Default is 0.1 (retain 10 percent, ie every 10th position)
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param cols Vector of colors, where each color corresponds to an identity class. By default, use the ggplot2 colors.
#'
#' @importFrom ggplot2 geom_bar facet_wrap xlab ylab theme_classic aes ylim theme
#' @export
#'
CoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 0.1,
  cells = NULL,
  cols = NULL
) {
  reads <- GetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    group.by = group.by,
    fragment.path = fragment.path,
    verbose = FALSE
  )
  coverages <- CalculateCoverages(
    reads = reads,
    window = window
  )
  if (downsample > 1) {
    warning('Requested downsampling <0%, retaining all positions')
    downsample <- 1
  }
  chromosome <- unlist(strsplit(region, ':'))[[1]]
  stepsize <- 1 / downsample
  total_range <- coverages[nrow(coverages), 'position'] - coverages[1, 'position']
  steps <- ceiling(total_range / stepsize)
  retain_positions <- seq(coverages[1, 'position'], coverages[nrow(coverages), 'position'], by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
  ymax <- round(max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)
  p <- ggplot(downsampled_coverage, aes(position, coverage, fill = group)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~group, strip.position = 'right', ncol = 1) +
    xlab(paste0(chromosome, ' position (bp)')) +
    ylab(paste0('Normalized coverage (range 0 - ', as.character(ymax), ')')) +
    ylim(c(0, ymax)) +
    theme_classic() +
    theme(axis.text.y = element_blank())
  return(p)
  # TODO add optional gene annotation track
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
  # need to go from set of annotations to list of positions (single base)
  # separate into + and - strand
  # get upstream / downstream for each
  # invert the - strand and add to the + strand reads
  # plot histogram with cell groupings (as for PeriodPlot)
  return()
}

#' Plot fragment length periodicity
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{SetFragments} function
#' @param region Genomic range to use. Default is fist megabase of chromosome 1.
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
  region = 'chr1:1-1000000',
  group.by = NULL,
  cells = NULL
) {
  reads <- GetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    group.by = group.by,
    fragment.path = fragment.path,
    verbose = FALSE
  )
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

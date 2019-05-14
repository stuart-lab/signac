
#' SingleCoveragePlot
#'
#' Plot coverage within given region for groups of cells
#'
#' Thanks to Andrew Hill for providing an early version of this function
#'  \url{http://andrewjohnhill.com/blog/2019/04/12/streamlining-scatac-seq-visualization-and-analysis/}
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show
#' @param annotation An Ensembl based annotation package
#' @param assay Which assay to plot
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{SetFragments} function
#' @param cells Which cells to plot. Default all cells
#' @param idents Which identities to include in the plot. Default is all identities.
#' @param window Smoothing window size
#' @param downsample Fraction of positions to retain in the plot. Default is 0.1 (retain 10 percent, ie every 10th position)
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#'
#' @importFrom ggplot2 geom_bar facet_wrap xlab ylab theme_classic aes ylim theme element_blank
#' @importFrom ggbio autoplot
#' @import patchwork
#' @importFrom AnnotationFilter GRangesFilter AnnotationFilterList GeneBiotypeFilter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @export
#'
SingleCoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 0.1,
  cells = NULL,
  idents = NULL
) {
  cells <- cells %||% colnames(x = object)
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
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
  pos <- unlist(strsplit(region, ':'))[[2]]
  start.pos <- as.numeric(unlist(strsplit(pos, '-'))[[1]])
  end.pos <- as.numeric(unlist(strsplit(pos, '-'))[[2]])

  stepsize <- 1 / downsample
  total_range <- end.pos - start.pos
  steps <- ceiling(total_range / stepsize)
  retain_positions <- seq(start.pos, end.pos, by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
  ymax <- signif(max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)

  p <- ggplot(downsampled_coverage, aes(position, coverage, fill = group)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~group, strip.position = 'right', ncol = 1) +
    xlab(paste0(chromosome, ' position (bp)')) +
    ylab(paste0('Normalized coverage (range 0 - ', as.character(ymax), ')')) +
    ylim(c(0, ymax)) +
    theme_classic() +
    theme(axis.text.y = element_blank(), legend.position = 'none')

  if (!is.null(annotation)) {
    gr <- GRanges(
      seqnames = gsub(pattern = 'chr', replacement = '', x = chromosome),
      IRanges(start = start.pos, end = end.pos)
    )
    filters <- AnnotationFilterList(GRangesFilter(value = gr), GeneBiotypeFilter('protein_coding'))
    genes <- suppressMessages(autoplot(EnsDb.Hsapiens.v75, filters, names.expr = 'gene_name'))
    gene.plot <- genes@ggplot +
      xlim(start.pos, end.pos) +
      xlab(paste0(chromosome, ' position (bp)')) +
      theme_classic()
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
      )
    p <- p + gene.plot + plot_layout(ncol = 1, heights = c(4, 1))
  }
  return(p)
}

#' CoveragePlot
#'
#' Plot coverage within given region for groups of cells
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show
#' @param ... Arguments passed to SingleCoveragePlot
#'
#' @importFrom patchwork wrap_plots
#' @export
#'
CoveragePlot <- function(
  object,
  region,
  ...
) {
  if (length(region) > 1) {
    plot.list <- lapply(region, SingleCoveragePlot, object = object, ...)
    return(wrap_plots(plot.list))
  } else {
    return(SingleCoveragePlot(object = object, region = region, ...))
  }
}

#' MotifPlot
#'
#' Plot motifs
#'
#' @param object A Seurat object
#' @param motifs A list of motifs to plot
#' @param assay Name of the assay to use
#' @param ... Additional parameters passed to \code{\link{ggseqlogo}}
#'
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom TFBSTools Matrix
#'
MotifPlot <- function(
  object,
  motifs,
  assay = NULL,
  ...
) {
  data.use <- GetMotifData(object = object, assay = assay, slot = 'pwm')
  if (length(x = data.use) == 0) {
    stop('Position weight matrix list for the requested assay is empty')
  }
  data.use <- data.use[motifs]
  if (class(x = data.use) == "PFMatrixList") {
    pwm <- Matrix(data.use)
    names(x = pwm) <- name(x = data.use)
  } else {
    pwm <- data.use
  }
  p <- ggseqlogo(data = pwm)
  return(p)
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
  annotation,
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

#' CoveragePlot
#'
#' @rdname CoveragePlot
#'
#' @importFrom ggplot2 geom_bar facet_wrap xlab ylab theme_classic aes ylim theme element_blank element_text
#' @importFrom ggbio autoplot
#' @import patchwork
#' @importFrom AnnotationFilter GRangesFilter AnnotationFilterList GeneBiotypeFilter
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
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
  extend.upstream = 0,
  extend.downstream = 0,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-")
) {
  cells <- cells %||% colnames(x = object)
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (class(x = region) != 'GRanges') {
    region <- StringToGRanges(regions = region, sep = sep)
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
    )
  )
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
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  stepsize <- 1 / downsample
  total_range <- end.pos - start.pos
  steps <- ceiling(x = (total_range / stepsize))
  retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
  ymax <- signif(x = max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)

  p <- ggplot(data = downsampled_coverage, mapping = aes(x = position, y = coverage, fill = group)) +
    geom_bar(stat = 'identity') +
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    xlab(label = paste0(chromosome, ' position (bp)')) +
    ylab(label = paste0('Normalized coverage (range 0 - ', as.character(x = ymax), ')')) +
    ylim(c(0, ymax)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      legend.position = 'none',
      strip.text.y = element_text(angle = 0)
    )

  if (!is.null(x = annotation)) {
    gr <- GRanges(
      seqnames = gsub(pattern = 'chr', replacement = '', x = chromosome),
      IRanges(start = start.pos, end = end.pos)
    )
    filters <- AnnotationFilterList(GRangesFilter(value = gr), GeneBiotypeFilter(value = 'protein_coding'))
    if (suppressMessages(expr = nrow(x = select(x = annotation, filters)) > 0)) {
      genes <- suppressMessages(expr = autoplot(object = annotation, filters, names.expr = 'gene_name'))
      gene.plot <- genes@ggplot +
        xlim(start.pos, end.pos) +
        xlab(label = paste0(chromosome, ' position (bp)')) +
        theme_classic()
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()
        )
      p <- p + gene.plot + plot_layout(ncol = 1, heights = c(4, 1))
    }
  }
  return(p)
}

#' CoveragePlot
#'
#' Plot coverage within given regions for groups of cells
#'
#' Thanks to Andrew Hill for providing an early version of this function
#'  \url{http://andrewjohnhill.com/blog/2019/04/12/streamlining-scatac-seq-visualization-and-analysis/}
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show. Can be a GRanges object, a string, or a vector of strings describing the genomic
#' coordinates to plot.
#' @param annotation An Ensembl based annotation package
#' @param assay Name of the  assay to plot
#' @param fragment.path Path to an index fragment file. If NULL, will look for a path stored for the
#' requested assay using the \code{\link{SetFragments}} function
#' @param cells Which cells to plot. Default all cells
#' @param idents Which identities to include in the plot. Default is all identities.
#' @param window Smoothing window size
#' @param downsample Fraction of positions to retain in the plot. Default is 0.1 (retain 10 percent, ie every 10th position)
#' @param extend.upstream Number of bases to extend the region upstream (Default 0)
#' @param extend.downstream Number of bases to extend the region downstream (Default 0)
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param ... Additional arguments passed to \code{\link[patchwork]{wrap_plots}}
#'
#' @importFrom patchwork wrap_plots
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
  extend.upstream = 0,
  extend.downstream = 0,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  ...
) {
  if (length(x = region) > 1) {
    plot.list <- lapply(
      X = region,
      FUN = SingleCoveragePlot,
      object = object,
      annotation = annotation,
      assay = assay,
      fragment.path = fragment.path,
      group.by = group.by,
      window = window,
      downsample = downsample,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      cells = cells,
      idents = idents,
      sep = sep
    )
    return(wrap_plots(plot.list, ...))
  } else {
    return(SingleCoveragePlot(
      object = object,
      region = region,
      annotation = annotation,
      assay = assay,
      fragment.path = fragment.path,
      group.by = group.by,
      window = window,
      downsample = downsample,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      cells = cells,
      idents = idents,
      sep = sep
    ))
  }
}

#' MotifDimPlot
#'
#' Plot motifs in reduced dimesions.
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param group.by A set of identities to group by (present in the Motif object metadata).
#' @param reduction Which dimension reduction to use. Default is tSNE.
#'
#' @importFrom Seurat Embeddings
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme_bw
#'
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
MotifDimPlot <- function(
  object,
  assay = NULL,
  group.by = NULL,
  reduction = 'tSNE',
  ...
) {
  coords.use <- GetMotifData(object = object, assay = assay, slot = 'reductions')
  if (!(reduction %in% names(x = coords.use))) {
    stop("Requested dimension reduction is not present")
  }
  coords.use <- as.data.frame(x = Embeddings(object = coords.use[[reduction]]))
  if (!is.null(x = group.by)) {
    meta.data <- GetMotifData(object = object, slot = 'meta.data')
    if (!(group.by %in% colnames(x = meta.data))) {
      stop("Requested grouping variable not present in Motif metadata")
    }
    coords.use[['ident']] <- meta.data[[group.by]]
  } else {
    coords.use[['ident']] <- 'Motif'
  }
  colnames(x = coords.use) <- c('dim1', 'dim2', 'ident')
  p <- ggplot(data = coords.use, mapping = aes(x = dim1, y = dim2, color = ident)) +
    geom_point() +
    xlab(label = paste0(reduction, '_1')) +
    ylab(label = paste0(reduction, '_2')) +
    theme_bw()
  return(p)
}

#' MotifPlot
#'
#' Plot motifs
#'
#' @param object A Seurat object
#' @param motifs A list of motifs to plot
#' @param assay Name of the assay to use
#' @param ... Additional parameters passed to \code{\link[ggseqlogo]{ggseqlogo}}
#'
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom TFBSTools name Matrix
#'
#' @export
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
    pwm <- Matrix(x = data.use)
    names(x = pwm) <- name(x = data.use)
  } else {
    pwm <- data.use
  }
  p <- ggseqlogo(data = pwm, ...)
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
  assay <- assay %||% DefaultAssay(object = object)
  fragment.path <- fragment.path %||% GetFragments(object = object, assay = assay)
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
#' @param region Genomic range to use. Default is fist two megabases of chromosome 1. Can be a GRanges object, a string, or a vector
#' of strings.
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param ... Additional arguments passed to \code{\link{GetReadsInRegion}}
#'
#' @importFrom ggplot2 ggplot geom_histogram theme_bw aes facet_wrap xlim
#'
#' @return Returns a ggplot2 object
#' @export
#'
PeriodPlot <- function(
  object,
  assay = NULL,
  fragment.path = NULL,
  region = 'chr1-1-2000000',
  group.by = NULL,
  cells = NULL,
  ...
) {
  reads <- GetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    group.by = group.by,
    fragment.path = fragment.path,
    verbose = FALSE,
    ...
  )
  if (length(x = unique(x = reads$group)) == 1) {
    p <- ggplot(data = reads, aes(length)) +
      geom_histogram(bins = 200) +
      xlim(c(0, 800)) +
      theme_bw()
  } else {
    p <- ggplot(data = reads, mapping = aes(x = length, fill = group)) +
      geom_histogram(bins = 200) +
      facet_wrap(~group, scales = 'free_y') +
      xlim(c(0, 800)) +
      theme_bw()
  }
  return(p)
}

#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

globalVariables(names = c('position', 'coverage', 'group'), package = 'Signac')
#' CoveragePlot
#'
#' @rdname CoveragePlot
#'
#' @importFrom ggplot2 geom_area geom_hline facet_wrap xlab ylab theme_classic aes ylim theme element_blank element_text
#' @importFrom ggbio autoplot
#' @importFrom cowplot plot_grid
#' @importFrom AnnotationFilter GRangesFilter AnnotationFilterList GeneBiotypeFilter
#' @importFrom AnnotationDbi select
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom Seurat WhichCells Idents
#' @importFrom Matrix colSums
#' @importFrom methods is
#' @importFrom stats median
#' @importFrom dplyr mutate group_by ungroup
#' @importFrom zoo rollapply
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
  height.tracks = 2,
  extend.upstream = 0,
  extend.downstream = 0,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-")
) {
  cells <- cells %||% colnames(x = object)
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (!is(object = region, class2 = 'GRanges')) {
    region <- StringToGRanges(regions = region, sep = sep)
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  reads.per.group <- AverageCounts(
    object = object,
    group.by = group.by,
    verbose = FALSE
  )
  cells.per.group <- CellsPerGroup(
    object = object,
    group.by = group.by
  )
  cutmat <- CutMatrix(
    object = object,
    region = region,
    cells = cells,
    verbose = FALSE
  )
  group.scale.factors <- reads.per.group * cells.per.group
  scale.factor <- scale.factor %||% median(x = group.scale.factors)
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  coverages <- ApplyMatrixByGroup(
    mat = cutmat,
    fun = colSums,
    groups = obj.groups,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    normalize = TRUE
  )
  if (!is.na(x = window)) {
    coverages <- group_by(.data = coverages, group)
    coverages <- mutate(.data = coverages, coverage = rollapply(
      data = norm.value,
      width = window,
      FUN = mean,
      align = 'center',
      fill = NA
    ))
    coverages <- ungroup(x = coverages)
  } else {
    coverages$coverage <- coverages$norm.value
  }
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  stepsize <- 1 / downsample
  total_range <- end.pos - start.pos
  steps <- ceiling(x = (total_range / stepsize))
  retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
  ymax <- ymax %||% signif(x = max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)
  ymin <- 0
  downsampled_coverage <- downsampled_coverage[!is.na(x = downsampled_coverage$coverage), ]

  p <- ggplot(data = downsampled_coverage, mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = 'identity') +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    xlab(label = paste0(chromosome, ' position (bp)')) +
    ylab(label = paste0('Normalized accessibility \n(range ', as.character(x = ymin), ' - ', as.character(x = ymax), ')')) +
    ylim(c(ymin, ymax)) +
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
      p <- suppressWarnings(plot_grid(
        p, gene.plot,
        ncol = 1,
        axis = 'btlr',
        rel_heights = c(height.tracks, 1),
        align = 'v',
        greedy = FALSE
      ))
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
#' @param height.tracks Height of the accessibility tracks relative to the height of the gene annotation track.
#' Default is 2 (twice as high as annotation track).
#' @param extend.upstream Number of bases to extend the region upstream (Default 0)
#' @param extend.downstream Number of bases to extend the region downstream (Default 0)
#' @param ymax Maximum value for Y axis. If NULL (default) set to the highest value among all the tracks.
#' @param scale.factor Scaling factor for track height. If NULL (default), use the median group scaling factor
#' determined by total number of fragments sequences in each group.
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param sep Separators to use for strings encoding genomic coordinates. First element is used to separate the
#' chromosome from the coordinates, second element is used to separate the start from end coordinate.
#' @param ... Additional arguments passed to \code{\link[cowplot]{plot_grid}}
#'
#' @importFrom cowplot plot_grid
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
  height.tracks = 2,
  extend.upstream = 0,
  extend.downstream = 0,
  scale.factor = NULL,
  ymax = NULL,
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
      ymax = ymax,
      scale.factor = scale.factor,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      cells = cells,
      idents = idents,
      sep = sep
    )
    return(plot_grid(plotlist = plot.list, ...))
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
      height.tracks = height.tracks,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      ymax = ymax,
      scale.factor = scale.factor,
      cells = cells,
      idents = idents,
      sep = sep
    ))
  }
}

globalVariables(names = c('dim1', 'dim2', 'ident'), package = 'Signac')
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
  reduction = 'tSNE'
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
#' @importFrom TFBSTools name
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
  if (is(object = data.use, class2 = "PFMatrixList")) {
    pwm <- TFBSTools::Matrix(x = data.use)
    names(x = pwm) <- name(x = data.use)
  } else {
    pwm <- data.use
  }
  p <- ggseqlogo(data = pwm, ...)
  return(p)
}

globalVariables(names = 'group', package = 'Signac')
#' Plot fragment length periodicity
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param region Genomic range to use. Default is fist two megabases of chromosome 1. Can be a GRanges object, a string, or a vector
#' of strings.
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities
#' @param log.scale Display Y-axis on log scale. Default is FALSE.
#' @param ... Additional arguments passed to \code{\link{GetReadsInRegion}}
#'
#' @importFrom ggplot2 ggplot geom_histogram theme_bw aes facet_wrap xlim scale_y_log10
#'
#' @return Returns a ggplot2 object
#' @export
#'
PeriodPlot <- function(
  object,
  assay = NULL,
  region = 'chr1-1-2000000',
  group.by = NULL,
  cells = NULL,
  log.scale = FALSE,
  ...
) {
  reads <- GetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    group.by = group.by,
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
  if (log.scale) {
    p <- p + scale_y_log10()
  }
  return(p)
}


#' Plot pileup of Tn5 integration sites
#'
#' Plots a pileup of integration sites centered on a set of genomic positions.
#' Each genomic region will be aligned on the region midpoint, and extended upstream
#' and downstream from the midpoint.
#'
#' @param object A Seurat object
#' @param assay Name of the assay to use
#' @param regions A set of GRanges to use
#' @param upstream Number of bases to extend upstream of the region midpoint
#' @param downstream Number of bases to extend downstream of the region midpoint
#' @param group.by A set of identities to group the cells by. Can by anything in the metadata.
#' Default is to use the active identities.
#' @param min.cells Minimum number of cells in the group for the pileup to be displayed for that group.
#' @param ymax Maximum value for the y-axis. If NULL (default), will be set automatically.
#' @param idents Which identities to include in the plot. If NULL (default), include everything with more than
#' \code{min.cells} cells.
#' @param verbose Display messages
#'
#' @importFrom BiocGenerics strand
#' @importFrom Seurat Idents
#' @importFrom Matrix colSums colMeans
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap ylim xlab ylab theme_classic theme element_blank element_text
#' @export
#' @return Returns a \code{\link[ggplot2]{ggplot2}} object
#' @examples
#' \dontrun{
#'
#' }
RegionPileup <- function(
  object,
  regions,
  assay = NULL,
  upstream = 200,
  downstream = 200,
  group.by = NULL,
  min.cells = 100,
  ymax = NULL,
  idents = NULL,
  verbose = TRUE
) {
  full.matrix <- CreateRegionPileupMatrix(
    object = object,
    regions = regions,
    upstream = upstream,
    downstream = downstream,
    assay = assay,
    cells = cells,
    verbose = verbose
  )
  # reads.per.group <- AverageCounts(
  #   object = object,
  #   group.by = group.by,
  #   verbose = FALSE
  # )
  # cells.per.group <- CellsPerGroup(
  #   object = object,
  #   group.by = group.by
  # )
  # group.scale.factors <- reads.per.group * cells.per.group
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  cellcounts <- table(obj.groups)
  obj.groups <- obj.groups[obj.groups %in% names(x = cellcounts[cellcounts > min.cells])]
  if (verbose) {
    message("Computing pileup for each cell group")
  }
  coverages <- ApplyMatrixByGroup(
    mat = full.matrix,
    groups = obj.groups,
    fun = colMeans,
    # group.scale.factors = group.scale.factors,
    # scale.factor = scale.factor,
    normalize = FALSE
  )
  ymin <- 0
  ymax <- ymax %||% signif(x = max(coverages$norm.value, na.rm = TRUE), digits = 2)
  p <- ggplot(data = coverages, mapping = aes(x = position, y = norm.value, color = group)) +
    geom_line(stat = 'identity', size = 0.2) +
    facet_wrap(facets = ~group) +
    xlab(label = paste0('Distance from region midpoint (bp)')) +
    ylab(label = paste0('Mean integration counts\n(0 - ', ymax, ')')) +
    ylim(c(ymin, ymax)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      legend.position = 'none',
      strip.text.y = element_text(angle = 0)
    )
  return(p)
}

#' Plot the enrichment around TSS
#'
#' Plot the normalized TSS enrichment score at each position relative to the TSS.
#' Requires that \code{\link{TSSEnrichment}} has already been run on the assay.
#'
#' @param object A Seurat object
#' @param assay Name of the assay to use. Should have the TSS enrichment information for each cell
#' already computed by running \code{\link{TSSEnrichment}}
#' @param group.by Set of identities to group cells by
#' @param idents Set of identities to include in the plot
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix colMeans
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab theme_minimal
#'
#' @return Returns a \code{\link[ggplot2]{ggplot2}} object
#' @export
TSSPlot <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  # get the normalized TSS enrichment matrix
  misc.slot <- GetAssayData(object = object, assay = assay, slot = 'misc')
  if (!(inherits(x = misc.slot, what = 'list'))) {
    stop("Misc slot does not contain a list")
  }
  if (!("TSS.enrichment.matrix" %in% names(x = misc.slot))) {
    stop("TSS enrichment matrix not present in assay. Run TSSEnrichment.")
  }
  tss.matrix <- misc.slot$TSS.enrichment.matrix

  # average the TSS score per group per base
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  groupmeans <- ApplyMatrixByGroup(
    mat = tss.matrix,
    groups = obj.groups,
    fun = colMeans,
    normalize = FALSE
  )

  p <- ggplot(
    data = groupmeans,
    mapping = aes(x = position, y = norm.value, color = group)
  ) +
    geom_line(stat = 'identity', size = 0.2) +
    facet_wrap(facets = ~group) +
    xlab("Distance from TSS (bp)") +
    ylab(label = 'Mean TSS enrichment score') +
    theme_minimal()
  return(p)
}


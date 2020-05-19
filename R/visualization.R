#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

globalVariables(names = c("Component", "counts"), package = "Signac")
#' Sequencing depth correlation
#'
#' Compute the correlation between total counts and each reduced
#' dimension component.
#'
#' @param object A \code{\link[Seurat]{Seurat}} object
#' @param reduction Name of a dimension reduction stored in the
#' input object
#' @param assay Name of assay to use for sequencing depth. If NULL, use the
#' default assay.
#' @param n Number of components to use. If \code{NULL}, use all components.
#' @param ... Additional arguments passed to \code{\link[stats]{cor}}
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @importFrom Seurat Embeddings DefaultAssay
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous
#' ylab ylim theme_light ggtitle aes
#' @importFrom stats cor
#' @concept visualization
#' @examples
#' DepthCor(object = atac_small)
DepthCor <- function(object, assay = NULL, reduction = 'lsi', n = 10, ...) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  dr <- object[[reduction]]
  embed <- Embeddings(object = dr)
  counts <- object[[paste0("nCount_", assay)]]
  embed <- embed[rownames(x = counts), ]
  n <- SetIfNull(x = n, y = ncol(x = embed))
  embed <- embed[, seq_len(length.out = n)]
  depth.cor <- as.data.frame(cor(x = embed, y = counts, ...))
  depth.cor$counts <- depth.cor[, 1]
  depth.cor$Component <- seq_len(length.out = nrow(x = depth.cor))
  p <- ggplot(depth.cor, aes(Component, counts)) +
    geom_point() +
    scale_x_continuous(n.breaks = n, limits = c(1, n)) +
    ylab("Correlation") +
    ylim(c(-1, 1)) +
    theme_light() +
    ggtitle("Correlation between depth and reduced dimension components",
            subtitle = paste0("Assay: ", assay, "\t", "Reduction: ", reduction))
  return(p)
}

globalVariables(
  names = c("feature", "group", "mn", "norm.value"),
  package = "Signac"
)
#' Plot footprinting results
#'
#' @param object A Seurat object
#' @param features A vector of features to plot
#' @param assay Name of assay to use
#' @param group.by A grouping variable
#' @param idents Set of identities to include in the plot
#' @param show.expected Plot the expected Tn5 integration frequency below the
#' main footprint plot
#' @param normalization Method to normalize for Tn5 DNA sequence bias. Options
#' are "subtract", "divide", or NULL to perform no bias correction.
#' @param label Label groups
#' @param repel Repel labels from each other
#' @param label.top Number of groups to label based on highest accessibility
#' in motif flanking region.
#' @export
#' @concept visualization
#' @concept footprinting
#' @importFrom Seurat DefaultAssay
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap xlab ylab theme_classic
#' theme element_blank geom_label guides guide_legend
#' @importFrom dplyr group_by summarize top_n
#' @import patchwork
PlotFootprint <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  label = TRUE,
  repel = TRUE,
  show.expected = TRUE,
  normalization = "subtract",
  label.top = 3
) {
  # TODO add option to show variance among cells
  # TODO update TSSPlot to use the GetFootprintData function
  plot.data <- GetFootprintData(
    object = object,
    features = features,
    assay = assay,
    group.by = group.by,
    idents = idents
  )
  motif.sizes <- GetMotifSize(
    object = object,
    features = features,
    assay = assay
  )
  obs <- plot.data[plot.data$class == "Observed", ]
  expect <- plot.data[plot.data$class == "Expected", ]

  # flanks are motif edge to 50 bp each side
  # add flank information (T/F)
  base <- ceiling(motif.sizes / 2)
  obs$flanks <- sapply(
    X = seq_len(length.out = nrow(x = obs)),
    FUN = function(x) {
      pos <- abs(obs[x, "position"])
      size <- base[[obs[x, "feature"]]]
      return((pos > size) & (pos < (size + 50)))
    })

  if (!is.null(normalization)) {
    # need to group by position and motif
    correction.vec <- expect$norm.value
    names(correction.vec) <- paste(expect$position, expect$feature)
    if (normalization == "subtract") {
      obs$norm.value <- obs$norm.value - correction.vec[
        paste(obs$position, obs$feature)
        ]
    } else if (normalization == "divide") {
      obs$norm.value <- obs$norm.value / correction.vec[
        paste(obs$position, obs$feature)
        ]
    } else {
      stop("Unknown normalization method requested")
    }
  }

  # find flanking accessibility for each group and each feature
  flanks <- obs[obs$flanks, ]
  flanks <- group_by(.data = flanks, feature, group)
  flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))

  # find top n groups for each feature
  topmean <- top_n(x = flankmeans, n = label.top, wt = mn)

  # find the top for each feature to determine axis limits
  ymax <- top_n(x = flankmeans, n = 1, wt = mn)
  ymin <- top_n(x = flankmeans, n = 1, wt =-mn)

  # make df for labels
  label.df <- data.frame()
  sub <- obs[obs$position == 50, ]
  for (i in seq_along(along.with = features)) {
    groups.use <- topmean[topmean$feature == features[[i]], ]$group
    df.sub <- sub[
      (sub$feature == features[[i]]) &
        (sub$group %in% groups.use), ]
    label.df <- rbind(label.df, df.sub)
  }
  obs$label <- NA
  label.df$label <- label.df$group
  obs <- rbind(obs, label.df)

  plotlist <- list()
  for (i in seq_along(along.with = features)) {
    # plot each feature separately rather than using facet
    # easier to manage the "expected" track
    df <- obs[obs$feature == features[[i]], ]
    min.use <- ifelse(test = normalization == "subtract", yes = -0.5, no = 0.5)
    axis.min <- min(min.use, ymin[ymin$feature == features[[i]], ]$mn)
    axis.max <- ymax[ymax$feature == features[[i]], ]$mn + 0.5

    p <- ggplot(
      data = df,
      mapping = aes(
        x = position,
        y = norm.value,
        color = group,
        label = label)
    )
    p <- p +
      geom_line(size = 0.2) +
      xlab("Distance from motif") +
      ylab(label = "Tn5 insertion\nenrichment") +
      theme_classic() +
      ggtitle(label = features[[i]]) +
      ylim(c(axis.min, axis.max)) +
      guides(color = guide_legend(override.aes = list(size = 1)))
    if (label) {
      if (repel) {
        p <- p + geom_label_repel(box.padding = 0.5, show.legend = FALSE)
      } else {
        p <- p + geom_label(show.legend = FALSE)
      }
    }
    if (show.expected) {
      df <- expect[expect$feature == features[[i]], ]
      p1 <- ggplot(
        data = df,
        mapping = aes(x = position, y = norm.value)
      ) +
        geom_line(size = 0.2) +
        xlab("Distance from motif") +
        ylab(label = "Expected\nTn5 enrichment") +
        theme_classic()

      # remove x-axis labels from top plot
      p <- p + theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x.bottom = element_blank()
      )
      p <- p + p1 + plot_layout(ncol = 1, heights = c(3, 1))
      plotlist[[i]] <- p
    }
  }
  plots <- wrap_plots(plotlist)
  return(plots)
}

globalVariables(
  names = c("position", "coverage", "group", "gene_name", "direction"),
  package = "Signac"
)
#' @rdname CoveragePlot
#' @importFrom ggplot2 geom_area geom_hline facet_wrap xlab ylab theme_classic
#' aes ylim theme element_blank element_text geom_segment scale_color_identity
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importMethodsFrom GenomicRanges start end
#' @importFrom Seurat WhichCells Idents DefaultAssay
#' @importFrom Matrix colSums
#' @importFrom methods is
#' @importFrom stats median
#' @importFrom dplyr mutate group_by ungroup
#' @importFrom zoo rollapply
#' @importFrom grid unit
#' @importFrom gggenes geom_gene_arrow geom_gene_label
#' @import patchwork
SingleCoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  peaks = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 0.1,
  height.tracks = 10,
  extend.upstream = 0,
  extend.downstream = 0,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-")
) {
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = annotation, what = "logical")) {
    annotation <- SetIfNull(
      x = annotation, y = Annotation(object = object[[assay]])
    )
  } else {
    annotation <- NULL
  }
  if (!inherits(x = peaks, what = "logical")) {
    peaks <- SetIfNull(x = peaks, y = granges(x = object[[assay]]))
  } else {
    peaks <- NULL
  }
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  if (!is(object = region, class2 = "GRanges")) {
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
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- start(x = region):end(x = region)
  group.scale.factors <- reads.per.group * cells.per.group
  scale.factor <- SetIfNull(
    x = scale.factor, y = median(x = group.scale.factors)
  )
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
      align = "center",
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
  retain_positions <- seq(from = start.pos, to = end.pos, by = stepsize)
  downsampled_coverage <- coverages[coverages$position %in% retain_positions, ]
  ymax <- SetIfNull(x = ymax, y = signif(
    x = max(downsampled_coverage$coverage, na.rm = TRUE), digits = 2)
  )
  ymin <- 0
  downsampled_coverage <- downsampled_coverage[!is.na(
    x = downsampled_coverage$coverage
  ), ]

  gr <- GRanges(
    seqnames = chromosome,
    IRanges(start = start.pos, end = end.pos)
  )
  p <- ggplot(
    data = downsampled_coverage,
    mapping = aes(x = position, y = coverage, fill = group)
    ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "right", ncol = 1) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("Normalized accessibility \n(range ",
                        as.character(x = ymin), " - ",
                        as.character(x = ymax), ")")) +
    ylim(c(ymin, ymax)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      legend.position = "none",
      strip.text.y = element_text(angle = 0)
    )
  if (!is.null(x = peaks)) {
    # subset to covered range
    peak.intersect <- subsetByOverlaps(x = peaks, ranges = gr)
    peak.df <- as.data.frame(x = peak.intersect)
    if (nrow(x = peak.df) > 0) {
      peak.plot <- ggplot(data = peak.df, mapping = aes(color = "darkgrey")) +
        geom_segment(aes(x = start, y = 0, xend = end, yend = 0, size = 2),
                     data = peak.df)
    } else {
      # no peaks present in region, make empty panel
      peak.plot <- ggplot(data = peak.df)
    }
    peak.plot <- peak.plot + theme_classic() +
      ylab(label = "Peaks") +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "none") +
      xlab(label = paste0(chromosome, " position (bp)")) +
      xlim(c(start.pos, end.pos)) +
      scale_color_identity()
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  } else {
    peak.plot <- NULL
  }
  if (!is.null(x = annotation)) {
    if (!inherits(x = annotation, what = 'GRanges')) {
      stop("Annotation must be a GRanges object or EnsDb object.")
    }
    annotation.subset <- subsetByOverlaps(x = annotation, ranges = gr)
    annotation.df <- as.data.frame(x = annotation.subset)
    # adjust coordinates so within the plot
    annotation.df$start[annotation.df$start < start.pos] <- start.pos
    annotation.df$end[annotation.df$end > end.pos] <- end.pos
    annotation.df$direction <- ifelse(
      test = annotation.df$strand == "-", yes = -1, no = 1
    )
    if (nrow(x = annotation.df) > 0) {
      gene.plot <- ggplot(
        data = annotation.df,
        mapping = aes(
          xmin = start,
          xmax = end,
          y = strand,
          fill = strand,
          label = gene_name,
          forward = direction)
        ) +
        geom_gene_arrow(
          arrow_body_height = unit(x = 4, units = "mm"),
          arrowhead_height = unit(x = 4, units = "mm"),
          arrowhead_width = unit(x = 5, units = "mm")) +
        geom_gene_label(
          grow = TRUE,
          reflow = TRUE,
          height = unit(x = 4, units = "mm")
        ) +
        xlim(start.pos, end.pos) +
        xlab(label = paste0(chromosome, " position (bp)")) +
        ylab("Genes") +
        theme_classic() +
        theme(legend.position = "none",
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank())

        # remove axis from coverage plot
        p <- p + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank()
        )
        if (!is.null(x = peak.plot)) {
          peak.plot <- peak.plot + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.ticks.x.bottom = element_blank()
          )
          p <- p +
            peak.plot +
            gene.plot +
            plot_layout(ncol = 1, heights = c(height.tracks, 1, 1))
        } else {
          p <- p +
            gene.plot +
            plot_layout(ncol = 1, heights = c(height.tracks, 1))
        }
    } else {
      if (!is.null(peak.plot)) {
        p <- p +
          peak.plot +
          plot_layout(ncol = 1, heights = c(height.tracks, 1))
      }
    }
  } else {
    if (!is.null(peak.plot)) {
      p <- p +
        peak.plot +
        plot_layout(ncol = 1, heights = c(height.tracks, 1))
    }
  }
  return(p)
}

#' Plot Tn5 insertion sites over a region
#'
#' Plot fragment coverage (frequence of Tn5 insertion) within given regions
#' for groups of cells.
#'
#' Thanks to Andrew Hill for providing an early version of this function
#' \url{http://andrewjohnhill.com/blog/2019/04/12/streamlining-scatac-seq-visualization-and-analysis/}
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show. Can be a GRanges object,
#' a string, or a vector of strings describing the genomic
#' coordinates to plot.
#' @param annotation A GRanges object containing genomic annotations. If NULL,
#' use the annotations stored in the assay. If NA or FALSE, don't plot
#' annotations.
#' @param peaks A GRanges object containing peak coordinates. If NULL, use the
#' genomic ranges associated with the assay. If NA or FALSE, don't plot genomic
#' ranges.
#' @param assay Name of the  assay to plot
#' @param fragment.path Path to an index fragment file. If NULL, will look for a
#' path stored in the fragments slot of the ChromatinAssay object
#' @param cells Which cells to plot. Default all cells
#' @param idents Which identities to include in the plot. Default is all
#' identities.
#' @param window Smoothing window size
#' @param downsample Fraction of positions to retain in the plot.
#' @param height.tracks Height of the accessibility tracks relative to the
#' height of the gene annotation track.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' @param ymax Maximum value for Y axis. If NULL (default) set to the highest
#' value among all the tracks.
#' @param scale.factor Scaling factor for track height. If NULL (default),
#' use the median group scaling factor determined by total number of fragments
#' sequences in each group.
#' @param group.by Name of one or more metadata columns to group (color) the
#' cells by. Default is the current cell identities
#' @param sep Separators to use for strings encoding genomic coordinates. First
#' element is used to separate the chromosome from the coordinates, second
#' element is used to separate the start from end coordinate.
#' @param ... Additional arguments passed to \code{\link[patchwork]{wrap_plots}}
#'
#' @importFrom patchwork wrap_plots
#' @export
#' @concept visualization
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' CoveragePlot(object = atac_small, region = c("chr1-713500-714500"))
#' }
CoveragePlot <- function(
  object,
  region,
  annotation = NULL,
  peaks = NULL,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  window = 100,
  downsample = 0.1,
  height.tracks = 10,
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
      X = seq_along(region),
      FUN = function(x) {
        SingleCoveragePlot(
          object = object,
          region = region[x],
          annotation = annotation,
          peaks = peaks,
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
      }
    )
    return(wrap_plots(plot.list, ...))
  } else {
    return(SingleCoveragePlot(
      object = object,
      region = region,
      annotation = annotation,
      peaks = peaks,
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

#' MotifPlot
#'
#' Plot motifs
#'
#' @param object A Seurat object
#' @param motifs A list of motifs to plot
#' @param assay Name of the assay to use
#' @param use.names Use motif names stored in the motif object
#' @param ... Additional parameters passed to \code{\link[ggseqlogo]{ggseqlogo}}
#'
#' @importFrom ggseqlogo ggseqlogo
#' @export
#' @concept visualization
#' @concept motifs
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @examples
#' \donttest{
#' motif.obj <- Seurat::GetAssayData(atac_small, slot = "motifs")
#' MotifPlot(atac_small, motifs = head(colnames(motif.obj)))
#' }
MotifPlot <- function(
  object,
  motifs,
  assay = NULL,
  use.names = TRUE,
  ...
) {
  data.use <- GetMotifData(object = object, assay = assay, slot = "pwm")
  if (length(x = data.use) == 0) {
    stop("Position weight matrix list for the requested assay is empty")
  }
  data.use <- data.use[motifs]
  if (use.names) {
    names(x = data.use) <- GetMotifData(
      object = object, assay = assay, slot = "motif.names"
    )[motifs]
  }
  p <- ggseqlogo(data = data.use, ...)
  return(p)
}

globalVariables(names = "group", package = "Signac")
#' Plot fragment length histogram
#'
#' @param object A Seurat object
#' @param assay Which assay to use. Default is the active assay.
#' @param region Genomic range to use. Default is fist two megabases of
#' chromosome 1. Can be a GRanges object, a string, or a vector
#' of strings.
#' @param cells Which cells to plot. Default all cells
#' @param group.by Name of one or more metadata columns to group (color) the
#' cells by. Default is the current cell identities
#' @param log.scale Display Y-axis on log scale. Default is FALSE.
#' @param ... Arguments passed to other functions
#'
#' @importFrom ggplot2 ggplot geom_histogram theme_bw aes facet_wrap xlim
#' scale_y_log10
#'
#' @export
#' @concept visualization
#' @concept qc
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' FragmentHistogram(object = atac_small, region = "chr1-10245-780007")
#' }
FragmentHistogram <- function(
  object,
  assay = NULL,
  region = "chr1-1-2000000",
  group.by = NULL,
  cells = NULL,
  log.scale = FALSE,
  ...
) {
  reads <- MultiGetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    verbose = FALSE,
    ...
  )
  # add group information
  if (is.null(x = group.by)) {
    groups <- Idents(object = object)
  } else {
    md <- object[[]]
    groups <- object[[group.by]]
    groups <- groups[, 1]
    names(x = groups) <- rownames(x = md)
  }
  reads$group <- groups[reads$cell]
  if (length(x = unique(x = reads$group)) == 1) {
    p <- ggplot(data = reads, aes(length)) +
      geom_histogram(bins = 200) +
      xlim(c(0, 800)) +
      theme_bw()
  } else {
    p <- ggplot(data = reads, mapping = aes(x = length, fill = group)) +
      geom_histogram(bins = 200) +
      facet_wrap(~group, scales = "free_y") +
      xlim(c(0, 800)) +
      theme_bw()
  }
  if (log.scale) {
    p <- p + scale_y_log10()
  }
  return(p)
}

globalVariables(names = "norm.value", package = "Signac")
#' Plot the enrichment of Tn5 integration sites centered on a set of genomic
#' regions.
#'
#' @param object A Seurat object
#' @param enrichment.key Name of a position enrichment matrix stored in the
#' \code{positionEnrichment} slot
#' of a \code{\link{ChromatinAssay}}.
#' @param assay Name of the assay to use. Must be a \code{\link{ChromatinAssay}}
#' and have the enrichment information for each cell stored in the
#' \code{positionEnrichment} slot.
#' @param group.by Set of identities to group cells by
#' @param idents Set of identities to include in the plot
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix colMeans
#' @importFrom ggplot2 ggplot aes geom_line
#'
#' @return Returns a \code{\link[ggplot2]{ggplot2}} object
#' @export
#' @concept visualization
EnrichmentPlot <- function(
  object,
  enrichment.key,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  # get the normalized TSS enrichment matrix
  positionEnrichment <- GetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment"
  )
  if (!(enrichment.key %in% names(x = positionEnrichment))) {
    stop("Position enrichment matrix not present in assay")
  }
  enrichment.matrix <- positionEnrichment[[enrichment.key]]

  # average the signal per group per base
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  groupmeans <- ApplyMatrixByGroup(
    mat = enrichment.matrix,
    groups = obj.groups,
    fun = colMeans,
    normalize = FALSE
  )

  p <- ggplot(
    data = groupmeans,
    mapping = aes(x = position, y = norm.value, color = group)
  ) +
    geom_line(stat = "identity", size = 0.2) +
    facet_wrap(facets = ~group)
  return(p)
}

#' Plot the enrichment around TSS
#'
#' Plot the normalized TSS enrichment score at each position relative to the
#' TSS. Requires that \code{\link{TSSEnrichment}} has already been run on the
#' assay.
#'
#' Wrapper for the \code{\link{EnrichmentPlot}} function
#'
#' @param object A Seurat object
#' @param assay Name of the assay to use. Should have the TSS enrichment
#' information for each cell
#' already computed by running \code{\link{TSSEnrichment}}
#' @param group.by Set of identities to group cells by
#' @param idents Set of identities to include in the plot
#'
#' @importFrom ggplot2 xlab ylab theme_minimal
#'
#' @return Returns a \code{\link[ggplot2]{ggplot2}} object
#' @export
#' @concept visualization
#' @concept qc
TSSPlot <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  p <- EnrichmentPlot(
    object = object,
    assay = assay,
    group.by = group.by,
    idents = idents,
    enrichment.key = "TSS"
  )
  p <- p +
    xlab("Distance from TSS (bp)") +
    ylab(label = "Mean TSS enrichment score") +
    theme_minimal()
  return(p)
}

#' Combine genome region plots
#'
#' This can be used to combine coverage plots, peak region plots, gene
#' annotation plots, and linked element plots. The different tracks are stacked
#' on top of each other and the x-axis combined.
#'
#' @param plotlist A list of plots to combine. Must be from the same genomic
#' region.
#' @return Returns a patchworked ggplot2 object
#' @export
#' @concept visualization
CombineTracks <- function(plotlist) {
  return()
}

#' Plot peaks in a genomic region
#'
#' Display the genomic ranges in a \code{\link{ChromatinAssay}} object that fall
#' in a given genomic region
#'
#' @param object A \code{\link[Seurat]{Seurat}} object
#' @param region A genomic region to plot
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @concept visualization
#' @importFrom GenomicRanges start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 ggplot aes geom_segment theme_classic element_blank
#' theme xlab ylab scale_color_identity
#' @examples
#' PeakPlot(atac_small, region = "chr1-710000-715000")
PeakPlot <- function(object, region) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  # get ranges from object
  peaks <- granges(x = object)
  # subset to covered range
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)

  if (nrow(x = peak.df) > 0) {
    peak.plot <- ggplot(data = peak.df, mapping = aes(color = "darkgrey")) +
      geom_segment(aes(x = start, y = 0, xend = end, yend = 0, size = 2),
                   data = peak.df)
  } else {
    # no peaks present in region, make empty panel
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() +
    ylab(label = "Peaks") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start.pos, end.pos)) +
    scale_color_identity()
  return(peak.plot)
}

#' Plot linked genomic elements
#'
#' Display links between pairs of genomic elements within a given region of the
#' genome.
#'
#' @param object A \code{\link[Seurat]{Seurat}} object
#' @param region A genomic region to plot
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 ggplot geom_hline geom_curve aes theme_classic ylim xlim
#' ylab theme element_blank
#' @concept visualization
LinkPlot <- function(object, region) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  chromosome <- seqnames(x = region)

  # extract link information
  links <- Links(object = object)

  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)

  # convert to dataframe
  link.df <- as.data.frame(x = links.keep)
  link.df$group <- as.factor(link.df$group)

  # plot
  p <- ggplot(data = link.df) +
    geom_hline(yintercept = 0, color = 'grey') +
    geom_curve(
      mapping = aes(x = start, y = 0, xend = end, yend = 0, alpha = score),
      curvature = 1/2
    ) +
    theme_classic() +
    ylim(c(-1, 0)) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    ylab("Links") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
  return(p)
}

#' Plot gene annotations
#'
#' Display gene annotations in a given region of the genome.
#'
#' @param object A \code{\link[Seurat]{Seurat}} object
#' @param region A genomic region to plot
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @concept visualization
AnnotationPlot <- function(object, region) {
  return()
}

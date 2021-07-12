#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

globalVariables(names = c("bin", "score"), package = "Signac")
#' Plot data from BigWig
#'
#' Create a BigWig track. Note that this function does not work on windows.
#'
#' @param region GRanges object specifying region to plot
#' @param bigwig Path to a bigwig file
#' @param smooth Number of bases to smooth data over (rolling mean). If NULL,
#' do not apply smoothing.
#' @param type Plot type. Can be one of "line", "heatmap", or "coverage"
#' @param y_label Y-axis label
#' @param max.downsample Minimum number of positions kept when downsampling.
#' Downsampling rate is adaptive to the window size, but this parameter will set
#' the minimum possible number of positions to include so that plots do not
#' become too sparse when the window size is small.
#' @param downsample.rate Fraction of positions to retain when downsampling.
#' Retaining more positions can give a higher-resolution plot but can make the
#' number of points large, resulting in larger file sizes when saving the plot
#' and a longer period of time needed to draw the plot.
#'
#' @importFrom ggplot2 ggplot aes_string geom_line geom_tile xlab ylab geom_area
#' scale_fill_viridis_c
#' @importFrom RcppRoll roll_mean
#' @importFrom GenomicRanges start end seqnames width
#' @importFrom dplyr slice_sample group_by mutate ungroup
#' @concept visualization
#'
#' @export
BigwigTrack <- function(
  region,
  bigwig,
  smooth = 200,
  type = "coverage",
  y_label = "Score",
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  possible_types <- c("line", "heatmap", "coverage")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  if (.Platform$OS.type == "windows") {
    message("BigwigTrack not supported on Windows")
    return(NULL)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  region_data <- rtracklayer::import(
    con = bigwig,
    which = region,
    as = "NumericList"
  )[[1]]
  if (!is.null(x = smooth)) {
    region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
  }
  region_data <- data.frame(
    position = start(x = region):end(x = region),
    score = region_data,
    stringsAsFactors = FALSE
  )
  window.size = width(x = region)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = region_data, n = sampling)
  p <- ggplot(
    data = coverages,
    mapping = aes_string(x = "position", y = "score")
  )
  if (type == "line") {
    p <- p + geom_line()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    region_data$bin <- floor(x = region_data$position / smooth)
    region_data <- group_by(region_data, bin)
    region_data <- mutate(region_data, score = mean(x = score))
    region_data <- ungroup(region_data)
    region_data <- unique(x = region_data[, c("bin", "score")])
    p <- ggplot(
      data = region_data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c()
  } else if (type == "coverage") {
    p <- p + geom_area()
  }
  chromosome <- as.character(x = seqnames(x = region))
  p <- p + theme_browser() +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label)
  return(p)
}

globalVariables(names = c("Component", "counts"), package = "Signac")
#' Plot sequencing depth correlation
#'
#' Compute the correlation between total counts and each reduced
#' dimension component.
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
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
#' \donttest{
#' DepthCor(object = atac_small)
#' }
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
#' Plot motif footprinting results
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
#' @param label TRUE/FALSE value to control whether groups are labeled.
#' @param repel Repel labels from each other
#' @param label.top Number of groups to label based on highest accessibility
#' in motif flanking region.
#' @param label.idents Vector of identities to label. If supplied,
#' \code{label.top} will be ignored.
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
  label.top = 3,
  label.idents = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  # TODO add option to show variance among cells
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
  sub <- obs[obs$position == 75, ]
  for (i in seq_along(along.with = features)) {
    if (is.null(x = label.idents)) {
      # determine which idents to label based on flanking accessibility
      groups.use <- topmean[topmean$feature == features[[i]], ]$group
    } else {
      # supplied list of idents to label
      groups.use <- label.idents
    }
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
#' @importFrom ggplot2 ylab scale_fill_manual unit element_text theme
#' @importMethodsFrom GenomicRanges start end
#' @importFrom Seurat WhichCells Idents DefaultAssay Idents<-
SingleCoveragePlot <- function(
  object,
  region,
  features = NULL,
  assay = NULL,
  show.bulk = FALSE,
  expression.assay = NULL,
  expression.slot = "data",
  annotation = TRUE,
  peaks = TRUE,
  peaks.group.by = NULL,
  ranges = NULL,
  ranges.group.by = NULL,
  ranges.title = "Ranges",
  region.highlight = NULL,
  links = TRUE,
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  group.by = NULL,
  window = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  heights = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("Requested assay is not a ChromatinAssay.")
  }
  if (!is.null(x = group.by)) {
    Idents(object = object) <- group.by
  }
  if (!is.null(x = idents)) {
    ident.cells <- WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
  }
  region <- FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
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
  group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
  scale.factor <- SetIfNull(
    x = scale.factor, y = median(x = group.scale.factors)
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  p <- CoverageTrack(
    cutmat = cutmat,
    region = region,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    window = window,
    ymax = ymax,
    obj.groups = obj.groups,
    region.highlight = region.highlight,
    downsample.rate = downsample.rate,
    max.downsample = max.downsample
  )
  # create bigwig tracks
  if (!is.null(x = bigwig)) {
    if (!inherits(x = bigwig, what = "list")) {
      warning("BigWig should be a list of file paths")
      bigwig <- list("bigWig" = bigwig)
    }
    if (length(x = bigwig.type == 1)) {
      bigwig.type <- rep(x = bigwig.type, length(x = bigwig))
    } else if (length(x = bigwig.type) != length(x = bigwig)) {
      stop("Must supply a bigWig track type for each bigWig file")
    }
    # iterate over list of files
    bigwig.tracks <- list()
    for (i in seq_along(along.with = bigwig)) {
      bigwig.tracks[[i]] <- BigwigTrack(
        region = region,
        bigwig = bigwig[[i]],
        y_label = names(x = bigwig)[[i]],
        type = bigwig.type[[i]]
      )
    }
  } else {
    bigwig.tracks <- NULL
  }
  if (!is.null(x = features)) {
    ex.plot <- ExpressionPlot(
      object = object,
      features = features,
      assay = expression.assay,
      idents = idents,
      group.by = group.by,
      slot = expression.slot
    )
    widths <- c(10, length(x = features))
  } else {
    ex.plot <- NULL
    widths <- NULL
  }
  if (annotation) {
    gene.plot <- AnnotationPlot(object = object[[assay]], region = region)
  } else {
    gene.plot <- NULL
  }
  if (links) {
    link.plot <- LinkPlot(object = object[[assay]], region = region)
  } else {
    link.plot <- NULL
  }
  if (peaks) {
    peak.plot <- PeakPlot(
      object = object,
      assay = assay,
      region = region,
      group.by = peaks.group.by
    )
  } else {
    peak.plot <- NULL
  }
  if (!is.null(x = ranges)) {
    range.plot <- PeakPlot(
      object = object,
      assay = assay,
      region = region,
      peaks = ranges,
      group.by = ranges.group.by,
      color = "brown3") +
      ylab(ranges.title)
  } else {
    range.plot <- NULL
  }
  if (tile) {
    # reuse cut matrix
    tile.df <- ComputeTile(
      cutmatrix = cutmat,
      groups = obj.groups,
      window = tile.size,
      n = tile.cells,
      order = "total"
    )
    tile.plot <- CreateTilePlot(
      df = tile.df,
      n = tile.cells
    )
  } else {
    tile.plot <- NULL
  }
  if (show.bulk) {
    object$bulk <- "All cells"
    reads.per.group <- AverageCounts(
      object = object,
      group.by = "bulk",
      verbose = FALSE
    )
    cells.per.group <- CellsPerGroup(
      object = object,
      group.by = "bulk"
    )
    bulk.scale.factor <- suppressWarnings(reads.per.group * cells.per.group)
    bulk.groups <- rep(x = "All cells", length(x = obj.groups))
    names(x = bulk.groups) <- names(x = obj.groups)

    bulk.plot <- CoverageTrack(
      cutmat = cutmat,
      region = region,
      group.scale.factors = bulk.scale.factor,
      scale.factor = scale.factor,
      window = window,
      ymax = ymax,
      obj.groups = bulk.groups,
      downsample.rate = downsample.rate,
      max.downsample = max.downsample
    ) +
      scale_fill_manual(values = "grey") +
      ylab("")
  } else {
    bulk.plot <- NULL
  }
  nident <- length(x = unique(x = obj.groups))
  bulk.height <- (1 / nident) * 10
  bw.height <- bulk.height * length(x = bigwig.tracks)
  heights <- SetIfNull(
    x = heights, y = c(10, bulk.height, bw.height, 10, 2, 1, 1, 3)
  )
  # pre-combine bigwig tracks
  if (!is.null(x = bigwig.tracks)) {
    bigwig.tracks <- CombineTracks(
      plotlist = bigwig.tracks,
      heights = rep(x = 1, length(x = bigwig.tracks))
    )
  }
  p <- CombineTracks(
    plotlist = list(p, bulk.plot, bigwig.tracks, tile.plot, gene.plot,
                    peak.plot, range.plot, link.plot),
    expression.plot = ex.plot,
    heights = heights,
    widths = widths
  ) & theme(
    legend.key.size = unit(x = 1/2, units = "lines"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
  return(p)
}

# Coverage Track
#
#' @importFrom ggplot2 geom_area geom_hline facet_wrap xlab ylab theme_classic
#' aes ylim theme element_blank element_text geom_segment scale_color_identity
#' scale_fill_manual geom_rect aes_string
#' @importFrom IRanges IRanges width
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Matrix colSums
#' @importFrom stats median
#' @importFrom dplyr mutate group_by ungroup slice_sample
#' @importFrom RcppRoll roll_sum
#' @importFrom methods is
#' @importFrom GenomicRanges GRanges
#' @importFrom scales hue_pal
#' @importFrom S4Vectors mcols
#' @importMethodsFrom GenomicRanges start end
CoverageTrack <- function(
  cutmat,
  region,
  group.scale.factors,
  scale.factor,
  obj.groups,
  ymax,
  downsample.rate,
  region.highlight = NULL,
  window = 100,
  max.downsample = 3000
) {
  window.size <- width(x = region)
  levels.use <- levels(x = obj.groups)
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
    coverages <- mutate(.data = coverages, coverage = roll_sum(
      x = norm.value, n = window, fill = NA, align = "center"
    ))
    coverages <- ungroup(x = coverages)
  } else {
    coverages$coverage <- coverages$norm.value
  }
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  coverages <- coverages[!is.na(x = coverages$coverage), ]
  coverages <- group_by(.data = coverages, group)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = coverages, n = sampling)

  # restore factor levels
  if (!is.null(x = levels.use)) {
    colors_all <- hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
    coverages$group <- factor(x = coverages$group, levels = levels.use)
  }
  ymax <- SetIfNull(x = ymax, y = signif(
    x = max(coverages$coverage, na.rm = TRUE), digits = 2)
  )
  ymin <- 0

  gr <- GRanges(
    seqnames = chromosome,
    IRanges(start = start.pos, end = end.pos)
  )
  p <- ggplot(
    data = coverages,
    mapping = aes(x = position, y = coverage, fill = group)
  ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0("Normalized accessibility \n(range ",
                        as.character(x = ymin), " - ",
                        as.character(x = ymax), ")")) +
    ylim(c(ymin, ymax)) +
    theme_browser(legend = FALSE) +
    theme(panel.spacing.y = unit(x = 0, units = "line"))
  if (!is.null(x = levels.use)) {
    p <- p + scale_fill_manual(values = colors_all)
  }
  if (!is.null(x = region.highlight)) {
    if (!inherits(x = region.highlight, what = "GRanges")) {
      warning("region.highlight must be a GRanges object")
    } else {
      md <- mcols(x = region.highlight)
      if ("color" %in% colnames(x = md)) {
        color.use <- md$color
      } else {
        color.use <- rep(x = "grey", length(x = region.highlight))
      }
      df <- data.frame(
        "start" = start(x = region.highlight),
        "end" = end(x = region.highlight),
        "color" = color.use
      )
      df$start <- ifelse(
        test = df$start < start.pos,
        yes = start.pos,
        no = df$start
      )
      df$end <- ifelse(
        test = df$end > end.pos,
        yes = end.pos,
        no = df$end
      )
      p <- p +
        geom_rect(
          data = df,
          inherit.aes = FALSE,
          aes_string(
            xmin = "start",
            xmax = "end",
            ymin = 0,
            ymax = ymax),
          fill = rep(x = df$color, length(x = unique(x = coverages$group))),
          color = "transparent",
          alpha = 0.2
        )
    }
  }
  return(p)
}

#' Plot Tn5 insertion frequency over a region
#'
#' Plot frequency of Tn5 insertion events for different groups of cells within
#' given regions of the genome.
#'
#' Thanks to Andrew Hill for providing an early version of this function.
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show. Can be a GRanges object,
#' a string encoding a genomic position, a gene name, or a vector of strings
#' describing the genomic coordinates or gene names to plot. If a gene name is
#' supplied, annotations must be present in the assay.
#' @param features A vector of features present in another assay to plot
#' alongside accessibility tracks (for example, gene names).
#' @param assay Name of the assay to plot
#' @param show.bulk Include coverage track for all cells combined (pseudo-bulk).
#' Note that this will plot the combined accessibility for all cells included in
#' the plot (rather than all cells in the object).
#' @param expression.assay Name of the assay containing expression data to plot
#' alongside accessibility tracks. Only needed if supplying \code{features}
#' argument.
#' @param expression.slot Name of slot to pull expression data from. Only needed
#' if supplying the \code{features} argument.
#' @param annotation Display gene annotations
#' @param peaks Display peaks
#' @param peaks.group.by Grouping variable to color peaks by. Must be a variable
#' present in the feature metadata. If NULL, do not color peaks by any variable.
#' @param ranges Additional genomic ranges to plot
#' @param ranges.group.by Grouping variable to color ranges by. Must be a
#' variable present in the metadata stored in the \code{ranges} genomic ranges.
#' If NULL, do not color by any variable.
#' @param ranges.title Y-axis title for ranges track. Only relevant if
#' \code{ranges} parameter is set.
#' @param region.highlight Region to highlight on the plot. Should be a GRanges
#' object containing the coordinates to highlight. By default, regions will be
#' highlighted in grey. To change the color of the highlighting, include a
#' metadata column in the GRanges object named "color" containing the color to
#' use for each region.
#' @param links Display links
#' @param tile Display per-cell fragment information in sliding windows.
#' @param tile.size Size of the sliding window for per-cell fragment tile plot
#' @param tile.cells Number of cells to display fragment information for in tile
#' plot.
#' @param bigwig List of bigWig file paths to plot data from. Files can be
#' remotely hosted. The name of each element in the list will determine the
#' y-axis label given to the track.
#' @param bigwig.type Type of track to use for bigWig files ("line", "heatmap",
#' or "coverage"). Should either be a single value, or a list of values giving
#' the type for each individual track in the provided list of bigwig files.
#' @param cells Which cells to plot. Default all cells
#' @param idents Which identities to include in the plot. Default is all
#' identities.
#' @param window Smoothing window size
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
#' @param heights Relative heights for each track (accessibility, gene
#' annotations, peaks, links).
#' @param max.downsample Minimum number of positions kept when downsampling.
#' Downsampling rate is adaptive to the window size, but this parameter will set
#' the minimum possible number of positions to include so that plots do not
#' become too sparse when the window size is small.
#' @param downsample.rate Fraction of positions to retain when downsampling.
#' Retaining more positions can give a higher-resolution plot but can make the
#' number of points large, resulting in larger file sizes when saving the plot
#' and a longer period of time needed to draw the plot.
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
#'
#' # Basic coverage plot
#' CoveragePlot(object = atac_small, region = c("chr1-713500-714500"))
#'
#' # Show additional ranges
#' ranges.show <- StringToGRanges("chr1-713750-714000")
#' CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), ranges = ranges.show)
#'
#' # Highlight region
#' CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), region.highlight = ranges.show)
#'
#' # Change highlight color
#' ranges.show$color <- "orange"
#' CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), region.highlight = ranges.show)
#'
#' # Show expression data
#' CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), features = "ELK1")
#' }
CoveragePlot <- function(
  object,
  region,
  features = NULL,
  assay = NULL,
  show.bulk = FALSE,
  expression.assay = "RNA",
  expression.slot = "data",
  annotation = TRUE,
  peaks = TRUE,
  peaks.group.by = NULL,
  ranges = NULL,
  ranges.group.by = NULL,
  ranges.title = "Ranges",
  region.highlight = NULL,
  links = TRUE,
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  heights = NULL,
  group.by = NULL,
  window = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  scale.factor = NULL,
  ymax = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  max.downsample = 3000,
  downsample.rate = 0.1,
  ...
) {
  if (length(x = region) > 1) {
    plot.list <- lapply(
      X = seq_along(region),
      FUN = function(x) {
        SingleCoveragePlot(
          object = object,
          region = region[x],
          features = features,
          expression.assay = expression.assay,
          expression.slot = expression.slot,
          show.bulk = show.bulk,
          annotation = annotation,
          peaks = peaks,
          peaks.group.by = peaks.group.by,
          ranges = ranges,
          ranges.group.by = ranges.group.by,
          ranges.title = ranges.title,
          region.highlight = region.highlight,
          assay = assay,
          links = links,
          tile = tile,
          tile.size = tile.size,
          tile.cells = tile.cells,
          bigwig = bigwig,
          bigwig.type = bigwig.type,
          group.by = group.by,
          window = window,
          ymax = ymax,
          scale.factor = scale.factor,
          extend.upstream = extend.upstream,
          extend.downstream = extend.downstream,
          cells = cells,
          idents = idents,
          sep = sep,
          heights = heights,
          max.downsample = max.downsample,
          downsample.rate = downsample.rate
        )
      }
    )
    return(wrap_plots(plot.list, ...))
  } else {
    return(SingleCoveragePlot(
      object = object,
      region = region,
      annotation = annotation,
      features = features,
      expression.assay = expression.assay,
      expression.slot = expression.slot,
      show.bulk = show.bulk,
      peaks = peaks,
      peaks.group.by = peaks.group.by,
      ranges = ranges,
      ranges.group.by = ranges.group.by,
      ranges.title = ranges.title,
      region.highlight = region.highlight,
      assay = assay,
      links = links,
      tile = tile,
      tile.size = tile.size,
      tile.cells = tile.cells,
      bigwig = bigwig,
      bigwig.type = bigwig.type,
      group.by = group.by,
      window = window,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      ymax = ymax,
      scale.factor = scale.factor,
      cells = cells,
      idents = idents,
      sep = sep,
      heights = heights,
      max.downsample = max.downsample,
      downsample.rate = downsample.rate
    ))
  }
}

#' Plot DNA sequence motif
#'
#' Plot position weight matrix or position frequency matrix for different DNA
#' sequence motifs.
#'
#' @param object A Seurat object
#' @param motifs A list of motifs to plot
#' @param assay Name of the assay to use
#' @param use.names Use motif names stored in the motif object
#' @param ... Additional parameters passed to \code{\link[ggseqlogo]{ggseqlogo}}
#'
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom Seurat DefaultAssay
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
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
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
#' Plot the frequency that fragments of different lengths are present for
#' different groups of cells.
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
#' @importFrom ggplot2 ggplot geom_histogram theme_classic aes facet_wrap xlim
#' scale_y_log10 theme element_blank
#' @importFrom Seurat DefaultAssay
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
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
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
      geom_histogram(bins = 200)
  } else {
    p <- ggplot(data = reads, mapping = aes(x = length, fill = group)) +
      geom_histogram(bins = 200) +
      facet_wrap(~group, scales = "free_y")
  }
  p <- p + xlim(c(0, 800)) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank()
    ) +
    xlab("Fragment length (bp)") +
    ylab("Count")
  if (log.scale) {
    p <- p + scale_y_log10()
  }
  return(p)
}

globalVariables(names = "norm.value", package = "Signac")
#' Plot signal enrichment around TSSs
#'
#' Plot the normalized TSS enrichment score at each position relative to the
#' TSS. Requires that \code{\link{TSSEnrichment}} has already been run on the
#' assay.
#'
#' @param object A Seurat object
#' @param assay Name of the assay to use. Should have the TSS enrichment
#' information for each cell
#' already computed by running \code{\link{TSSEnrichment}}
#' @param group.by Set of identities to group cells by
#' @param idents Set of identities to include in the plot
#'
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom Matrix colMeans
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab theme_classic ggtitle
#' theme element_blank
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
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  # get the normalized TSS enrichment matrix
  positionEnrichment <- GetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment"
  )
  if (!("TSS" %in% names(x = positionEnrichment))) {
    stop("Position enrichment matrix not present in assay")
  }
  enrichment.matrix <- positionEnrichment[["TSS"]]

  # remove motif and expected
  if (nrow(x = enrichment.matrix) == (ncol(x = object) + 2)) {
    enrichment.matrix <- enrichment.matrix[1:(nrow(x = enrichment.matrix) - 2), ]
  }

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
    facet_wrap(facets = ~group) +
    xlab("Distance from TSS (bp)") +
    ylab(label = "Mean TSS enrichment score") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank()
    ) +
    ggtitle("TSS enrichment")
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
#' @param expression.plot Plot containing gene expression information. If
#' supplied, this will be placed to the left of the coverage tracks and aligned
#' with each track
#' @param heights Relative heights for each plot. If NULL, the first plot will
#' be 8x the height of the other tracks.
#' @param widths Relative widths for each plot. Only required if adding a gene
#' expression panel. If NULL, main plots will be 8x the width of the gene
#' expression panel
#' @return Returns a patchworked ggplot2 object
#' @export
#' @importFrom ggplot2 theme element_blank
#' @importFrom patchwork wrap_plots plot_layout guide_area
#' @concept visualization
#' @examples
#' \donttest{
#' p1 <- PeakPlot(atac_small, region = "chr1-29554-39554")
#' p2 <- AnnotationPlot(atac_small, region = "chr1-29554-39554")
#' CombineTracks(plotlist = list(p1, p2), heights = c(1, 1))
#' }
CombineTracks <- function(
  plotlist,
  expression.plot = NULL,
  heights = NULL,
  widths = NULL
) {
  # remove any that are NULL
  nullplots <- sapply(X = plotlist, FUN = is.null)
  plotlist <- plotlist[!nullplots]
  heights <- heights[!nullplots]

  if (length(x = plotlist) == 1) {
    return(plotlist[[1]])
  }

  # remove x-axis from all but last plot
  for (i in 1:(length(x = plotlist) - 1)) {
    plotlist[[i]] <- plotlist[[i]] + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
  }

  # combine plots
  if (is.null(x = heights)) {
    # set height of first element to 10x more than other elements
    n.plots <- length(x = plotlist)
    heights <- c(8, rep(1, n.plots - 1))
  } else {
    if (length(x = heights) != length(x = plotlist)) {
      stop("Relative height must be supplied for each plot")
    }
  }
  if (!is.null(x = expression.plot)) {
    # align expression plot with the first element in plot list
    p <- (plotlist[[1]] + expression.plot) +
      plot_layout(widths = widths)

    n <- length(x = plotlist)
    heights.2 <- heights[2:n]
    p2 <- wrap_plots(plotlist[2:n], ncol = 1, heights = heights.2)

    p <- p + p2 + guide_area() + plot_layout(
      ncol = 2, heights = c(heights[[1]], sum(heights.2)),
      guides = "collect")
  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
  }
  return(p)
}

#' Plot peaks in a genomic region
#'
#' Display the genomic ranges in a \code{\link{ChromatinAssay}} object that fall
#' in a given genomic region
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param region A genomic region to plot
#' @param peaks A GRanges object containing peak coordinates. If NULL, use
#' coordinates stored in the Seurat object.
#' @param group.by Name of variable in feature metadata (if using ranges in the
#' Seurat object) or genomic ranges metadata (if using supplied ranges) to color
#' ranges by. If NULL, do not color by any metadata variable.
#' @param color Color to use. If \code{group.by} is not NULL, this can be a
#' custom color scale (see examples).
#'
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @concept visualization
#' @importFrom Seurat DefaultAssay
#' @importFrom S4Vectors mcols<-
#' @importFrom GenomicRanges start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 ggplot aes geom_segment theme_classic element_blank
#' theme xlab ylab scale_color_manual
#' @examples
#' \donttest{
#' # plot peaks in assay
#' PeakPlot(atac_small, region = "chr1-710000-715000")
#'
#' # manually set color
#' PeakPlot(atac_small, region = "chr1-710000-715000", color = "red")
#'
#' # color by a variable in the feature metadata
#' PeakPlot(atac_small, region = "chr1-710000-715000", group.by = "count")
#' }
PeakPlot <- function(
  object,
  region,
  assay = NULL,
  peaks = NULL,
  group.by = NULL,
  color = "dimgrey"
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }

  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  if (is.null(x = peaks)) {
    peaks <- granges(x = object[[assay]])
    md <- object[[assay]][[]]
    mcols(x = peaks) <- md
  }
  # subset to covered range
  peak.intersect <- subsetByOverlaps(x = peaks, ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)

  if (nrow(x = peak.df) > 0) {
    if (!is.null(x = group.by)) {
      if (!(group.by %in% colnames(x = peak.df))) {
        warning("Requested grouping variable not found")
        group.by <- NULL
      }
    }
    peak.df$start[peak.df$start < start.pos] <- start.pos
    peak.df$end[peak.df$end > end.pos] <- end.pos
    peak.plot <- ggplot(
      data = peak.df,
      aes_string(color = SetIfNull(x = group.by, y = "color"))
    ) +
      geom_segment(aes(x = start, y = 0, xend = end, yend = 0),
                   size = 2,
                   data = peak.df)
  } else {
    # no peaks present in region, make empty panel
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() +
    ylab(label = "Peaks") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start.pos, end.pos))
  if (is.null(x = group.by)) {
    # remove legend, change color
    peak.plot <- peak.plot +
      scale_color_manual(values = color) +
      theme(legend.position = "none")
  }
  return(peak.plot)
}

#' Plot linked genomic elements
#'
#' Display links between pairs of genomic elements within a given region of the
#' genome.
#'
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param region A genomic region to plot
#' @param min.cutoff Minimum absolute score for link to be plotted.
#'
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 ggplot geom_hline aes theme_classic xlim
#' ylab theme element_blank scale_color_gradient2 aes_string
#' @importFrom ggforce geom_bezier
#' @concept visualization
#' @concept links
LinkPlot <- function(object, region, min.cutoff = 0) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  chromosome <- seqnames(x = region)

  # extract link information
  links <- Links(object = object)

  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }

  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)

  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  link.df <- link.df[abs(x = link.df$score) > min.cutoff, ]

  # remove links outside region
  link.df <- link.df[link.df$start >= start(x = region) & link.df$end <= end(x = region), ]

  # plot
  if (nrow(x = link.df) > 0) {
    # convert to format for geom_bezier
    link.df$group <- seq_len(length.out = nrow(x = link.df))
    df <- data.frame(
      x = c(link.df$start,
            (link.df$start + link.df$end) / 2,
            link.df$end),
      y = c(rep(x = 0, nrow(x = link.df)),
            rep(x = -1, nrow(x = link.df)),
            rep(x = 0, nrow(x = link.df))),
      group = rep(x = link.df$group, 3),
      score = rep(link.df$score, 3)
    )
    p <- ggplot(data = df) +
      geom_bezier(
        mapping = aes_string(x = "x", y = "y", group = "group", color = "score")
      ) +
      geom_hline(yintercept = 0, color = 'grey') +
      scale_color_gradient2(low = "red", mid = "grey", high = "blue")
  } else {
    p <- ggplot(data = link.df)
  }
  p <- p +
    theme_classic() +
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
#' @param object A \code{\link[SeuratObject]{Seurat}} object
#' @param region A genomic region to plot
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 theme_classic ylim xlim ylab xlab
#' geom_segment geom_text aes_string scale_color_manual
#' @importFrom grid arrow
#' @importFrom S4Vectors split
#' @importFrom fastmatch fmatch
#' @concept visualization
#' @examples
#' \donttest{
#' AnnotationPlot(object = atac_small, region = c("chr1-29554-39554"))
#' }
AnnotationPlot <- function(object, region) {
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)

  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_name)
  annotation.subset <- annotation[
    fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
  ]

  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 5
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.05, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.05, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + (n_stack * 0.2)
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes_string(x = "position", y = "dodge", label = "gene_name"),
      size = 3
    )
    y_limit <- c(0.9, n_stack + (n_stack * 0.5))
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (kb)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}

globalVariables(names = "gene", package = "Signac")
#' Plot gene expression
#'
#' Display gene expression values for different groups of cells and different
#' genes. Genes will be arranged on the x-axis and different groups stacked on
#' the y-axis, with expression value distribution for each group shown as a
#' violin plot. This is designed to work alongside a genomic coverage track,
#' and the plot will be able to be aligned with coverage tracks for the same
#' groups of cells.
#'
#' @param object A Seurat object
#' @param features A list of features to plot
#' @param assay Name of the assay storing expression information
#' @param group.by A grouping variable to group cells by. If NULL, use the
#' current cell identities
#' @param idents A list of identities to include in the plot. If NULL, include
#' all identities
#' @param slot Which slot to pull expression data from
#'
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom ggplot2 ggplot geom_violin facet_wrap aes theme_classic theme
#' element_blank scale_y_discrete scale_x_continuous scale_fill_manual
#' @importFrom scales hue_pal
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges start end
#' @importFrom patchwork wrap_plots
#' @importFrom fastmatch fmatch
#'
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' ExpressionPlot(atac_small, features = "TSPAN6", assay = "RNA")
#' }
ExpressionPlot <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  slot = "data"
) {
  # get data
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  data.plot <- GetAssayData(
    object = object,
    assay = assay,
    slot = slot
  )[features, ]
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = NULL
  )
  # if levels set, define colors based on all groups
  levels.use <- levels(x = obj.groups)
  if (!is.null(x = levels.use)) {
    colors_all <- hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
  }
  if (!is.null(x = idents)) {
    cells.keep <- names(x = obj.groups)[
      fmatch(x = obj.groups, table = idents, nomatch = 0L) > 0
    ]
    if (length(x = features) > 1) {
      data.plot <- data.plot[, cells.keep]
    } else {
      data.plot <- data.plot[cells.keep]
    }
    obj.groups <- obj.groups[cells.keep]
  }
  # construct data frame
  if (length(x = features) == 1) {
    df <- data.frame(
      gene = features,
      expression = data.plot,
      group = obj.groups
    )
  } else {
    df <- data.frame()
    for (i in features) {
      df.1 <- data.frame(
        gene = i,
        expression = data.plot[i, ],
        group = obj.groups
      )
      df <- rbind(df, df.1)
    }
  }
  p.list <- list()
  for (i in seq_along(along.with = features)) {
    df.use <- df[df$gene == features[[i]], ]
    p <- ggplot(data = df.use, aes(x = expression, y = gene, fill = group)) +
      geom_violin(size = 1/4) +
      facet_wrap(~group, ncol = 1, strip.position = "right") +
      theme_classic() +
      scale_y_discrete(position = "top") +
      scale_x_continuous(position = "bottom", limits = c(0, NA)) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none"
      )
    if (!is.null(x = levels.use)) {
      p <- p + scale_fill_manual(values = colors_all)
    }
    p.list[[i]] <- p
  }
  p <- wrap_plots(p.list, ncol = length(x = p.list))
  return(p)
}

#' Genome browser
#'
#' Interactive version of the \code{\link{CoveragePlot}} function. Allows
#' altering the genome position interactively. The current view at any time can
#' be saved to a list of \code{\link[ggplot2]{ggplot}} objects using the "Save
#' plot" button, and this list of plots will be returned after ending the
#' browser by pressing the "Done" button.
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates
#' @param assay Name of assay to use
#' @param sep Separators for genomic coordinates if region supplied as a string
#' rather than GRanges object
#' @param ... Parameters passed to \code{\link{CoveragePlot}}
#'
#' @return Returns a list of ggplot objects
#'
#' @export
#' @concept visualization
CoverageBrowser <- function(
  object,
  region,
  assay = NULL,
  sep = c("-", "-"),
  ...
) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Please install shiny. https://shiny.rstudio.com/")
  }
  if (!requireNamespace("miniUI", quietly = TRUE)) {
    stop("Please install miniUI. https://github.com/rstudio/miniUI")
  }

  if (!is(object = region, class2 = "GRanges")) {
    if (all(sapply(X = sep, FUN = grepl, x = region))) {
      region <- StringToGRanges(regions = region, sep = sep)
    } else {
      region <- LookupGeneCoords(object = object, assay = assay, gene = region)
    }
  }

  # work out group_by options
  grouping_options <- object[[]]
  group_types <- sapply(X = grouping_options, FUN = class)
  grouping_options <- colnames(x = grouping_options)[
    group_types %in% c("factor", "character")
  ]
  grouping_options <- c("Idents", grouping_options)

  startpos <- start(x = region)
  endpos <- end(x = region)
  chrom <- seqnames(x = region)

  ui <- miniUI::miniPage(

    miniUI::miniTabstripPanel(

      miniUI::miniTabPanel(

        title = "View",
        icon = shiny::icon("area-chart"),

        miniUI::miniButtonBlock(

          shiny::actionButton(
            inputId = "save",
            label = "Save plot",
            style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
          ),

          shiny::actionButton(
            inputId = "up_large",
            label = "",
            icon = shiny::icon(name = "angle-double-left")
          ),

          shiny::actionButton(
            inputId = "up_small",
            label = "",
            icon = shiny::icon(name = "angle-left")
          ),

          shiny::actionButton(
            inputId = "minus",
            label = "-"
          ),

          shiny::actionButton(
            inputId = "plus",
            label = "+"
          ),

          shiny::actionButton(
            inputId = "down_small",
            label = "",
            icon = shiny::icon(name = "angle-right")
          ),

          shiny::actionButton(
            inputId = "down_large",
            label = "",
            icon = shiny::icon(name = "angle-double-right")
          ),

          shiny::actionButton(
            inputId = "done",
            label = "Done",
            style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
          ),

          border = "bottom"
        ),

        miniUI::miniContentPanel(
          shiny::plotOutput(outputId = "access", height = "100%")
        )
      ),

      miniUI::miniTabPanel(

        title = "Parameters",
        icon = shiny::icon(name = "sliders"),

        miniUI::miniContentPanel(
          shiny::fillCol(
            flex = c(NA, 1), height = "10%",
            shiny::actionButton(
              inputId = "go",
              label = "Update plot",
              style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
            )
          ),

          shiny::fillRow(
            shiny::fillCol(
              flex = c(1, 4),

              shiny::titlePanel(title = "Navigation"),

              shiny::wellPanel(
                shiny::textInput(
                  inputId = "chrom",
                  label = "Chromosome",
                  value = chrom
                ),

                shiny::numericInput(
                  inputId = "startpos",
                  label = "Start",
                  value = startpos,
                  step = 5000
                ),

                shiny::numericInput(
                  inputId = "endpos",
                  label = "End",
                  value = endpos,
                  step = 5000
                ),

                shiny::textInput(
                  inputId = "gene_lookup",
                  label = "Gene",
                  value = NULL
                )
              )
            ),

          shiny::fillCol(
            flex = c(1, 4),

            shiny::titlePanel(title = "Tracks"),

            shiny::wellPanel(
              shiny::checkboxGroupInput(
                inputId = "tracks",
                label = "Display",
                choices = c("Annotations" = "genes",
                            "Peaks" = "peaks",
                            "Links" = "links",
                            "Tile" = "tile"),
                selected = c("genes", "peaks", "links")
              ),

              shiny::selectInput(
                inputId = "group_by",
                label = "Group by",
                choices = grouping_options,
                selected = NULL
              )

            )
          ),

          shiny::fillCol(
            flex = c(1, 4),

            shiny::titlePanel(title = "Gene expression"),

            shiny::wellPanel(
              shiny::textInput(
                inputId = "gene",
                label = "Gene",
                value = NULL
              ),

              shiny::textInput(
                inputId = "expression_assay",
                label = "Assay",
                value = "RNA"
              ),

              shiny::selectInput(
                inputId = "expression_slot",
                label = "Slot",
                choices = c("data", "scale.data", "counts")
                )
            )
            )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {

    # listen for change in any of the inputs
    changed <- shiny::reactive({
      paste(
        input$up_small,
        input$up_large,
        input$minus,
        input$plus,
        input$down_small,
        input$down_large,
        input$go
      )
    })

    # determine which tracks to show
    tracks <- shiny::reactiveValues(
      annotation = TRUE,
      peaks = TRUE,
      links = TRUE
    )

    shiny::observeEvent(
      eventExpr = input$tracks,
      handlerExpr = {
        tracks$annotation <- "genes" %in% input$tracks
        tracks$peaks <- "peaks" %in% input$tracks
        tracks$links <- "links" %in% input$tracks
        tracks$tile <- "tile" %in% input$tracks
      }
    )

    # list of reactive values for storing and modifying current coordinates
    coords <- shiny::reactiveValues(
      chromosome = chrom,
      startpos = startpos,
      endpos = endpos,
      width = endpos - startpos
    )

    # manually set coordinates
    shiny::observeEvent(
      eventExpr = input$chrom,
      handlerExpr = {
        coords$chromosome <- input$chrom
      }
    )
    shiny::observeEvent(
      eventExpr = input$startpos,
      handlerExpr = {
        coords$startpos <- input$startpos
        coords$width <- coords$endpos - coords$startpos
      }
    )
    shiny::observeEvent(
      eventExpr = input$endpos,
      handlerExpr = {
        coords$endpos <- input$endpos
        coords$width <- coords$endpos - coords$startpos
      }
    )

    # set coordinates based on gene name
    shiny::observeEvent(
      eventExpr = input$go,
      handlerExpr = {
        if (input$gene_lookup != "") {
          region <- LookupGeneCoords(
            object = object,
            assay = assay,
            gene = input$gene_lookup
          )
          if (!is.null(x = region)) {
            coords$chromosome <- seqnames(x = region)
            coords$startpos <- start(x = region)
            coords$endpos <- end(x = region)
            coords$width <- coords$endpos - coords$startpos
          }
        }
      }
    )

    # set gene expression panel
    gene_expression <- shiny::reactiveValues(
      genes = NULL,
      assay = "RNA",
      slot = "data"
    )
    shiny::observeEvent(
      eventExpr = input$go,
      handlerExpr = {
        if (is.null(x = input$gene)) {
          gene_expression$genes <- NULL
        } else if (nchar(x = input$gene) == 0) {
          gene_expression$genes <- NULL
        } else {
          gene_expression$genes <- unlist(x = strsplit(x = input$gene, split = ", "))
        }
      }
    )
    shiny::observeEvent(
      eventExpr = input$expression_assay,
      handlerExpr = {
        gene_expression$assay <- input$expression_assay
      }
    )
    shiny::observeEvent(
      eventExpr = input$expression_slot,
      handlerExpr = {
        gene_expression$slot <- input$expression_slot
      }
    )

    # scroll upstream
    shiny::observeEvent(
      eventExpr = input$up_large,
      handlerExpr = {
        coords$startpos <- coords$startpos - coords$width
        coords$endpos <- coords$endpos - coords$width
        coords$width <- coords$endpos - coords$startpos
      }
    )
    shiny::observeEvent(
      eventExpr = input$up_small,
      handlerExpr = {
        coords$startpos <- coords$startpos - (coords$width / 2)
        coords$endpos <- coords$endpos - (coords$width / 2)
        coords$width <- coords$endpos - coords$startpos
      }
    )

    # scroll downstream
    shiny::observeEvent(
      eventExpr = input$down_large,
      handlerExpr = {
        coords$startpos <- coords$startpos + coords$width
        coords$endpos <- coords$endpos + coords$width
        coords$width <- coords$endpos - coords$startpos
      }
    )
    shiny::observeEvent(
      eventExpr = input$down_small,
      handlerExpr = {
        coords$startpos <- coords$startpos + (coords$width / 2)
        coords$endpos <- coords$endpos + (coords$width / 2)
        coords$width <- coords$endpos - coords$startpos
      }
    )

    # zoom
    shiny::observeEvent(
      eventExpr = input$minus,
      handlerExpr = {
        coords$startpos <- coords$startpos - (coords$width / 4)
        coords$endpos <- coords$endpos + (coords$width / 4)
        coords$width <- coords$endpos - coords$startpos
      }
    )
    shiny::observeEvent(
      eventExpr = input$plus,
      handlerExpr = {
        coords$startpos <- coords$startpos + (coords$width / 4)
        coords$endpos <- coords$endpos - (coords$width / 4)
        coords$width <- coords$endpos - coords$startpos
      }
    )

    # update current region
    current_region <- shiny::eventReactive(
        eventExpr = changed(),
      valueExpr = GRanges(
        seqnames = coords$chromosome,
        ranges = IRanges(start = coords$startpos, end = coords$endpos)
      ),
      ignoreNULL = FALSE
    )

    # grouping
    group_by <- shiny::reactiveVal(value = NULL)
    shiny::observeEvent(
      eventExpr = input$go,
      handlerExpr = {
        if (input$group_by == "Idents") {
          group_by(NULL)
        } else {
          group_by(input$group_by)
        }
      }
    )

    current_plot <- shiny::reactiveVal(value = NULL)
    plot_list <- shiny::reactiveVal(value = list())

    shiny::observeEvent(
      eventExpr = input$save,
      handlerExpr = {
        n <- length(x = plot_list()) + 1
        p.list <- plot_list()
        p.list[[n]] <- current_plot()
        plot_list(p.list)
      }
    )

    shiny::observeEvent(
      eventExpr = changed(),
      handlerExpr = {
        p <- CoveragePlot(
          object = object,
          region = current_region(),
          group.by = group_by(),
          annotation = tracks$annotation,
          peaks = tracks$peaks,
          links = tracks$links,
          tile = tracks$tile,
          features = gene_expression$genes,
          expression.slot = gene_expression$slot,
          expression.assay = gene_expression$assay,
          max.downsample = 2000,
          downsample.rate = 0.02,
          ...
        )
        current_plot(p)
      }
    )

    output$access <- shiny::renderPlot(expr = {
      current_plot()
    })

    shiny::observeEvent(
      eventExpr = input$done,
      handlerExpr = {
        shiny::stopApp(returnValue = plot_list())
      })
  }

  shiny::runGadget(
    app = ui,
    server = server
  )
}

#' Plot strand concordance vs. VMR
#'
#' Plot the Pearson correlation between allele frequencies on each strand
#' versus the log10 mean-variance ratio for the allele.
#'
#' @param variants A dataframe containing variant information. This should be
#' computed using \code{\link{IdentifyVariants}}
#' @param min.cells Minimum number of high-confidence cells detected with the
#' variant for the variant to be displayed.
#' @param concordance.threshold Strand concordance threshold
#' @param vmr.threshold Mean-variance ratio threshold
#'
#' @concept mito
#' @concept visualization
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point labs scale_y_log10
#' geom_vline geom_hline theme_classic scale_color_manual theme
#' @importFrom scales comma
VariantPlot <- function(
  variants,
  min.cells = 2,
  concordance.threshold = 0.65,
  vmr.threshold = 0.01
) {
  high.conf <- variants[variants$n_cells_conf_detected >= min.cells, ]
  high.conf$pos <- high.conf$vmr > vmr.threshold &
    high.conf$strand_correlation > concordance.threshold
  p <- ggplot(
    data = high.conf,
    mapping = aes_string(x = "strand_correlation", y = "vmr", color = "pos")
    ) +
    geom_point() +
    labs(x = "Strand concordance", y = "Variance-mean ratio") +
    geom_vline(
      xintercept = concordance.threshold, color = "black", linetype = 2
      ) +
    geom_hline(
      yintercept = vmr.threshold, color = "black", linetype = 2
      ) +
    scale_color_manual(values = c("black", "firebrick")) +
    scale_y_log10(labels = comma) +
    theme_classic() +
    theme(legend.position = "none")
  return(p)
}


#' Plot integration sites per cell
#'
#' Plots the presence/absence of Tn5 integration sites for each cell
#' within a genomic region.
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show. Can be a GRanges object,
#' a string encoding a genomic position, a gene name, or a vector of strings
#' describing the genomic coordinates or gene names to plot. If a gene name is
#' supplied, annotations must be present in the assay.
#' @param assay Name of assay to use
#' @param group.by Name of grouping variable to group cells by. If NULL, use the
#' current cell identities
#' @param idents List of cell identities to include in the plot. If NULL, use
#' all identities.
#' @param sep Separators to use for strings encoding genomic coordinates. First
#' element is used to separate the chromosome from the coordinates, second
#' element is used to separate the start from end coordinate.
#' @param tile.size Size of the sliding window for per-cell fragment tile plot
#' @param tile.cells Number of cells to display fragment information for in tile
#' plot.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' @param order.by Option for determining how cells are chosen from each group.
#' Options are "total" or "random". "total" will select the top cells based on
#' total number of fragments in the region, "random" will select randomly.
#' @param cells Which cells to plot. Default all cells
#'
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
#' @importFrom Seurat DefaultAssay
#' @importFrom ggplot2 xlab
#' @importFrom GenomeInfoDb seqnames
#'
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' TilePlot(object = atac_small, region = c("chr1-713500-714500"))
#' }
TilePlot <- function(
  object,
  region,
  sep = c("-", "-"),
  tile.size = 100,
  tile.cells = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  assay = NULL,
  cells = NULL,
  group.by = NULL,
  order.by = "total",
  idents = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("Requested assay is not a ChromatinAssay.")
  }
  region <- FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  cutmat <- CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- start(x = region):end(x = region)
  tile.df <- ComputeTile(
    cutmatrix = cutmat,
    groups = obj.groups,
    window = tile.size,
    n = tile.cells,
    order = order.by
  )
  tile.plot <- CreateTilePlot(df = tile.df, n = tile.cells)
  tile.plot <- tile.plot +
    xlab(label = paste0(seqnames(x = region), " position (bp)"))
  return(tile.plot)
}

# Create tile plot for a region, given a pre-computed cutmatrix
# this is designed for use inside the main CoveragePlot function,
# avoids re-computing the cut matrix
#
# @param cutmatrix A sparse matrix of Tn5 integration sites per cell
# @param groups A grouping vector describing the group that each cell belong to
# @param idents Identities to include. If NULL, use all groups
# @param order Option for determining how cells are chosen from each group.
# Options are "total" or "random". "total" will select the top cells based on
# total number of fragments in the region, "random" will select randomly.
# @param n Number of cells to choose from each group
# @param window Size of sliding window. Integration events are summed within
# each window.
#
# @return Returns a ggplot2 object
#
#' @importFrom RcppRoll roll_sum
#' @importFrom Matrix rowSums
#' @importFrom tidyr pivot_longer
#' @importFrom utils head
ComputeTile <- function(
  cutmatrix,
  groups,
  window = 200,
  n = 100,
  idents = NULL,
  order = "total"
) {
  # for each group, choose n cells based on total counts
  totals <- rowSums(x = cutmatrix)
  unique.groups <- unique(x = groups)
  cells.use <- vector(mode = "character")
  cell.idx <- vector(mode = "numeric")
  for (i in seq_along(along.with = unique.groups)) {
    tot.use <- totals[names(x = groups[groups == unique.groups[[i]]])]
    cell.keep <- names(x = head(x = sort(x = tot.use, decreasing = TRUE), n))
    cell.idx <- c(cell.idx, seq_along(along.with = cell.keep))
    cells.use <- c(cells.use, cell.keep)
  }
  names(x = cell.idx) <- cells.use
  cutmatrix <- cutmatrix[cells.use, ]

  # create sliding window sum of integration sites using RcppRoll
  # note that this coerces to a dense matrix
  smoothed <- apply(
    X = cutmatrix,
    MARGIN = 1,
    FUN = roll_sum,
    n = window,
    by = window
  )

  # create dataframe
  smoothed <- as.data.frame(x = smoothed)

  # add extra column as bin ID
  smoothed$bin <- seq_len(length.out = nrow(x = smoothed))
  smoothed <- pivot_longer(
    data = smoothed,
    cols = cells.use
  )

  smoothed$group <- groups[smoothed$name]
  smoothed$idx <- cell.idx[smoothed$name]
  smoothed$bin <- smoothed$bin + as.numeric(x = colnames(x = cutmatrix)[[1]])
  return(smoothed)
}

#' @importFrom ggplot2 ggplot aes_string geom_raster ylab scale_fill_gradient
#' scale_y_reverse guides guide_legend geom_hline
CreateTilePlot <- function(df, n, legend = TRUE) {
  # create plot
  p <- ggplot(
    data = df,
    aes_string(x = "bin", y = "idx", fill = "value")) +
    facet_wrap(
      facets = ~group,
      scales = "free_y",
      ncol = 1,
      strip.position = "left"
    ) +
    geom_raster() +
    theme_browser(legend = legend) +
    geom_hline(yintercept = c(0, n), size = 0.1) +
    ylab(paste0("Fragments (", n, " cells)")) +
    scale_fill_gradient(low = "white", high = "darkred") +
    scale_y_reverse() +
    guides(fill = guide_legend(
      title = "Fragment\ncount",
      keywidth = 1/2, keyheight = 1
      )
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  return(p)
}

#' Genome browser theme
#'
#' Theme applied to panels in the \code{\link{CoveragePlot}} function.
#'
#' @param ... Additional arguments
#' @param legend Display plot legend
#'
#' @importFrom ggplot2 theme theme_classic element_blank element_text
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' PeakPlot(atac_small, region = "chr1-710000-715000") + theme_browser()
#' }
theme_browser <- function(..., legend = TRUE) {
  browser.theme <- theme_classic() +
    theme(
      axis.text.y = element_blank(),
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0)
    )
  if (!legend) {
    browser.theme <- browser.theme +
      theme(
        legend.position = "none"
      )
  }
  return(browser.theme)
}

####################
### Not exported ###
####################

split_body <- function(df, width = 1000) {
  wd <- df$end - df$start
  nbreak <- wd / width
  if (nbreak > 1) {
    steps <- 0:(nbreak)
    starts <- (width * steps) + df$start
    starts[starts > df$end] <- NULL
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    gene_name = df$gene_name[[1]],
    gene_biotype = df$gene_biotype[[1]],
    type = "arrow"
  )
  return(breaks)
}

reformat_annotations <- function(
  annotation,
  start.pos,
  end.pos
) {
  annotation <- annotation[annotation$type == "exon"]
  exons <- as.data.frame(x = annotation)
  annotation <- split(
    x = annotation,
    f = annotation$gene_name
  )
  annotation <- lapply(X = annotation, FUN = as.data.frame)

  # add gene total start / end
  gene_bodies <- list()
  for (i in seq_along(annotation)) {
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      gene_name = annotation[[i]]$gene_name[[1]],
      gene_biotype = annotation[[i]]$gene_biotype[[1]],
      type = "body"
    )
    # trim any that extend beyond region
    df$start <- ifelse(
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
    )
    df$end <- ifelse(
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
    )
    breaks <- split_body(df = df)
    df <- rbind(df, breaks)
    gene_bodies[[i]] <- df
  }
  gene_bodies <- do.call(what = rbind, args = gene_bodies)

  # record if genes overlap
  overlap_idx <- record_overlapping(annotation = gene_bodies, min.gapwidth = 1000)
  gene_bodies$dodge <- overlap_idx[gene_bodies$gene_name]
  exons$dodge <- overlap_idx[exons$gene_name]

  label_df <- gene_bodies[gene_bodies$type == "body", ]
  label_df$width <- label_df$end - label_df$start
  label_df$position <- label_df$start + (label_df$width / 2)

  onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
  onminus <- gene_bodies[gene_bodies$strand == "-", ]

  return(
    list(
      "labels" = label_df,
      "exons" = exons,
      "plus" = onplus,
      "minus" = onminus
    )
  )
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
record_overlapping <- function(annotation, min.gapwidth = 1000) {
  # convert back to granges
  annotation.stash <- annotation
  annotation$strand <- "*"
  gr <- makeGRangesFromDataFrame(
    df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
  )
  # work out which ranges overlap
  collapsed <- reduce(
    x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
  )$revmap
  idx <- seq_along(gr)
  for (i in seq_along(collapsed)) {
    mrg <- collapsed[[i]]
    for (j in seq_along(mrg)) {
      idx[[mrg[[j]]]] <- j
    }
  }
  names(x = idx) <- gr$gene_name
  return(idx)
}

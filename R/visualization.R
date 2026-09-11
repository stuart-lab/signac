#' @include generics.R
#' @importFrom utils globalVariables
#'
NULL

globalVariables(names = c("bin", "score", "bw"), package = "Signac")
#' Plot multiple coverage plots
#' @param object A [SeuratObject::Seurat] object
#' @param regions List of genes or a
#' [GenomicRanges::GRanges()] object to plot
#' @param region_names List of plot titles for each region.
#' If NULL, no title is added.
#' @param extend.upstream Number of bases to extend the region
#' upstream. Default is 0. Can be a vector the same length as
#' `regions`, or a single integer. If supplied with a single
#' integer, extend upstream the same number of bases for
#' every region in `regions`.
#' @param extend.downstream Number of bases to extend the
#' region downstream. Follows the same logic as
#' `extend.upstream`, but applying to downstream bases.
#' @param region.highlight List of regions to highlight for
#' each plotted region, should be the same length as `regions`.
#' Each item should be a [GenomicRanges::GRanges] object
#' containing the coordinates to highlight, or `NULL` if no
#' highlight is needed for that particular region. By default,
#' regions will be highlighted in grey. To change the color of
#' the highlighting, include a metadata column in the GRanges
#' object named "color" containing the color to use for each
#' region.
#' @param assay Name of the assay to plot. If a list of assays is provided,
#' data from each assay will be shown overlaid on each track. The first assay in
#' the list will define the assay used for gene annotations, links, and peaks
#' (if shown). The order of assays given defines the plotting order.
#' @param split.assays When plotting data from multiple assays, display each
#' assay as a separate track. If FALSE, data from different assays are overlaid
#' on a single track with transparency applied.
#' @param assay.scale Scaling to apply to data from different
#' assays. Can be:
#'  - common: plot all assays on a common scale (default)
#'  - separate: plot each assay on a separate scale ranging
#' from zero to the maximum value for that assay within the
#' plotted region
#' @param annotation Display gene annotations. Set to `TRUE`
#' or `FALSE` to control whether genes models are displayed,
#' or choose "transcript" to display all transcript isoforms,
#' or "gene" to display gene models only (same as setting
#' TRUE).
#' @param peaks Display peaks
#' @param group.by Name of one or more metadata columns to
#' group (color) the cells by. Default is the current cell
#' identities
#' @param idents Which identities to include in the plot.
#' Default is all identities.
#' @param split.by A metadata variable to split the tracks by.
#' For example, grouping by "celltype" and splitting by
#' "batch" will create separate tracks for each combination of
#' celltype and batch.
#' @param cells Which cells to plot. Default all cells
#' @param show.bulk Include coverage track for all cells
#' combined (pseudo-bulk). Note that this will plot the
#' combined accessibility for all cells included in the plot
#' (rather than all cells in the object).
#' @param ranges_list List of additional genomic ranges to
#' plot for each region, should be the same length as
#' `regions`. Each item should be a [GenomicRanges::GRanges]
#' object or `NULL` if no additional ranges is needed for a
#' region. Can also be a single GRanges object to be plotted
#' for every region.
#' @param ranges.title Y-axis title for ranges track. Only
#' relevant if `ranges_list` parameter is set.
#' @param max.downsample Maximum amount of downsampling to apply. Corresponds to
#' the minimum number of positions to be kept when downsampling. The
#' downsampling rate is adaptive to the window size, but this parameter will set
#' the minimum possible number of positions to include so that plots do not
#' become too sparse when the window size is small.
#' @param downsample.rate Fraction of positions to retain when
#' downsampling. Retaining more positions can give a
#' higher-resolution plot but can make the number of points
#' large, resulting in larger file sizes when saving the plot
#' and a longer period of time needed to draw the plot.
#' @param scale.factor Scaling factor for track height. If
#' NULL (default), use the median group scaling factor
#' determined by total number of fragments sequences in each
#' group.
#' @param ymax Maximum value for Y axis. Can be one of:
#'  - `NULL`: set to the highest value among all the tracks
#' (default)
#'  - qXX: clip the maximum value to the XX quantile (for
#' example, q95 will set the maximum value to 95% of the
#' maximum value in the data). This can help remove the effect
#' of extreme values that may otherwise distort the scale.
#'  - numeric: manually define a Y-axis limit
#' @param window Smoothing window size
#' @param bigwig A named list or named character vector of bigWig file paths to
#' plot data from. Files can be remotely hosted. The name of each element
#' determines the y-axis label given to the track.
#' @param bigwig.type Type of track to use for bigWig files
#' ("line" or "coverage"). Should either be a single value,
#' or a list of values giving the type for each individual
#' track in the provided list of bigwig files.
#' @param bigwig.scale Same as `assay.scale` parameter,
#' except for bigWig files when plotted with
#' `bigwig.type="coverage"`
#' @param title_size Size of plot title, only applicable if
#' `region_names` is not NULL. Default is 6
#' @param subtitle_size Size of plot subtitle (coverage plot
#' range). Default is 5
#' @param coords_size Size of plot coordinate text. Default
#' is 6
#' @param links Character vector containing the keys of link 
#' information present in the assay to display. Default is 
#' "linkpeaks" which is the default key for peak-gene links 
#' stored using the [LinkPeaks()] function. If `NULL`, links
#' will not be displayed. Color gradient will take into account 
#' all links in `regions` and a common legend will be shown.
#' @return Returns a ggplot object
#' @export
#' @concept visualization
MultiCoveragePlot <- function(
  object,
  regions,
  region_names = NULL,
  extend.upstream = NULL,
  extend.downstream = NULL,
  region.highlight = NULL,
  assay = "peaks",
  assay.scale = "common",
  split.assays = FALSE,
  annotation = TRUE,
  peaks = TRUE,
  group.by = NULL,
  idents = NULL,
  split.by = NULL,
  cells = NULL,
  show.bulk = FALSE,
  ranges_list = NULL,
  ranges.title = "Ranges",
  max.downsample = 3000,
  downsample.rate = 0.1,
  scale.factor = NULL,
  ymax = NULL,
  window = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  bigwig.scale = "common",
  title_size = 6,
  subtitle_size = 5,
  coords_size = 6,
  links = "linkpeaks"
) {
  # check valid.assay.scale
  valid.assay.scale <- c("common", "separate")
  if (!(assay.scale %in% valid.assay.scale)) {
    stop(
      "Unknown assay.scale requested. Please choose from: ",
      paste(valid.assay.scale, collapse = ", ")
    )
  }

  # check bigwig.type. This can be one value per bigwig file, so every element
  # is checked rather than relying on a length-1 condition
  enabled.bigwig.type <- c("line", "coverage")
  if (!all(unlist(x = bigwig.type) %in% enabled.bigwig.type)) {
    stop(
      paste0(
        "Unknown bigwig.type requested. Please choose from: ",
        paste(enabled.bigwig.type, collapse = ", "),
        ". bigwig.type=`heatmap` is not enabled ",
        "for MultiCoveragePlot"
      )
    )
  }

  # check assay
  cells <- cells %||% Cells(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = assay, what = "list")) {
    assay <- list(assay)
  }
  lapply(X = assay, FUN = function(x) {
    if (!inherits(x = object[[x]], what = "ChromatinAssay5")) {
      stop("Requested assay is not a ChromatinAssay5.")
    }
  })

  # check if requested link key is present
  if (!is.null(links)) {
    obj.links <- Links(object = object[[assay[[1]]]])
    if (is.null(obj.links) || is.null(obj.links[[links]])) {
      links <- NULL
    }
  }

  # get region highlights from regions
  if (!is.null(x = region.highlight)) {
    stopifnot(
      "region.highlight must be a list of the same length as regions" =
        length(x = region.highlight) == length(regions)
    )
    region.to.highlight <- region.highlight
  } else {
    region.to.highlight <- vector("list", length(regions))
  }

  # get plotting regions from regions
  regions.to.plot <- c()
  for (i in seq_along(regions)) {
    regions.to.plot[[i]] <- FindRegion(
      object = object,
      region = regions[[i]],
      assay = assay[[1]]
    )
  }

  # get additional ranges from ranges_list
  if (!is.null(ranges_list)) {
    if (is(ranges_list, "GRanges")) {
      ranges.to.plot <- vector("list", length(regions))
      for (i in seq_along(ranges.to.plot)) {
        ranges.to.plot[[i]] <- ranges_list
      }
    } else if (is(ranges_list, "list")) {
      stopifnot(
        "ranges_list must be a list of the same length as regions" =
          length(ranges_list) == length(regions)
      )
      for (i in seq_along(ranges_list)) {
        if (is.null(ranges_list[[i]])) {
          ranges_list[[i]] <- GRanges()
        }
      }
      ranges.to.plot <- ranges_list
    } else {
      stop(
        "Invalid ranges_list, choose either a single ",
        "GRanges object to plot for every plot or a ",
        "list of GRanges objects for each plot"
      )
    }
  } else {
    ranges.to.plot <- vector("list", length(regions))
  }

  # check region_names
  if (!is.null(region_names) &&
        length(region_names) != length(regions)) {
    stop("Supplied region list and region names differ in length.")
  }

  # upstream/downstream extension
  if (is.null(extend.upstream)) {
    extend.upstream <- rep(0, length(regions))
  } else if (length(extend.upstream) == 1) {
    extend.upstream <- rep(extend.upstream, length(regions))
  }
  if (length(extend.upstream) != length(regions)) {
    stop(
      "Supplied region list and extend.upstream ",
      "vector differ in length."
    )
  }

  if (is.null(extend.downstream)) {
    extend.downstream <- rep(0, length(regions))
  } else if (length(extend.downstream) == 1) {
    extend.downstream <- rep(extend.downstream, length(regions))
  }
  if (length(extend.downstream) != length(regions)) {
    stop(
      "Supplied region list and extend.downstream ",
      "vector differ in length."
    )
  }

  # call SingleCoveragePlot
  single.plots <- c()
  region.names_list <- as.list(region_names)
  for (i in seq_along(regions.to.plot)) {
    single.plots[[i]] <- SingleCoveragePlot(
      object = object,
      region = regions.to.plot[[i]],
      extend.upstream = extend.upstream[[i]],
      extend.downstream = extend.downstream[[i]],
      region.highlight = region.to.highlight[[i]],
      assay = assay,
      split.assays = split.assays,
      assay.scale = assay.scale,
      annotation = annotation,
      peaks = peaks,
      group.by = group.by,
      idents = idents,
      split.by = split.by,
      cells = cells,
      show.bulk = show.bulk,
      ranges = ranges.to.plot[[i]],
      ranges.title = ranges.title,
      max.downsample = max.downsample,
      downsample.rate = downsample.rate,
      scale.factor = scale.factor,
      ymax = ymax,
      window = window,
      bigwig = bigwig,
      bigwig.type = bigwig.type,
      bigwig.scale = bigwig.scale,
      links = links
    )
  }

  # check number of plots
  plot_params <- list(
    show.bulk = show.bulk,
    bigwig = !is.null(bigwig),
    annotation = !identical(annotation, FALSE),
    peaks = peaks,
    ranges.to.plot = !is.null(ranges_list),
    links = !is.null(links)
  )
  n_plots <- 1 + sum(unlist(lapply(plot_params, isTRUE)))

  # check bigwig plot position
  if (!is.null(bigwig)) {
    n_true <- sum(unlist(lapply(plot_params, isTRUE)))
    if (isTRUE(plot_params$show.bulk)) {
      if (n_true == 2) {
        bw_idx <- 1
      } else {
        bw_idx <- 3
      }
    } else if (n_true == 1) {
      bw_idx <- 1
    } else {
      bw_idx <- 2
    }
  }

  # rearrange plots
  arranged.plots <- c()
  for (i in seq_along(single.plots)) {
    if (n_plots > 1) {
      ## coverage.track
      covplot <- single.plots[[i]]$patches$plots[[1]]

      y_label <- covplot@labels$y
      range <- sub(".*range ([^)]*).*", "\\1", y_label)
      if (i == 1) {
        range <- paste0(range, " (range)")
      }

      # move region name to title, range to subtitle
      covplot <- covplot +
        labs(
          title = if (!is.null(region_names)) {
            region.names_list[[i]]
          },
          subtitle = range
        ) +
        theme(
          plot.title = element_text(
            size = title_size, hjust = 0.5
          ),
          plot.subtitle = element_text(
            size = subtitle_size, hjust = 1
          )
        )

      # adjust y axis label
      covplot@labels$y <- "Normalized accessibility"

      if (!is.null(links)) {
        covplot <- covplot + guides(color = "none", fill = "none")
      }
        
      single.plots[[i]]$patches$plots[[1]] <- covplot

      # remove any plot legends
      for (j in 1:(n_plots - 1)) {
        currentplot <- single.plots[[i]]$patches$plots[[j]]
        currentplot <- currentplot + theme(
          legend.position = "none"
        )
        single.plots[[i]]$patches$plots[[j]] <- currentplot
      }
      single.plots[[i]] <- single.plots[[i]] + theme(
        legend.position = "none"
      )

    } else if (n_plots == 1) {
      y_label <- single.plots[[i]]@labels$y
      range <- sub(".*range ([^)]*).*", "\\1", y_label)
      if (i == 1) {
        range <- paste0(range, " (range)")
      }

      # adjust covplot title, subtitle
      single.plots[[i]] <- single.plots[[i]] +
        labs(
          title = if (!is.null(region_names)) {
            region.names_list[[i]]
          },
          subtitle = range
        ) +
        theme(
          plot.title = element_text(
            size = title_size, hjust = 0.5
          ),
          plot.subtitle = element_text(
            size = subtitle_size, hjust = 1
          )
        )

      # adjust y axis label
      single.plots[[i]]@labels$y <- "Normalized accessibility"
    }

    arranged.plots[[i]] <- single.plots[[i]]
  }

  # remove y axis text from 2nd plot on
  for (i in seq_len(max(0, length(arranged.plots) - 1)) + 1) {
    if (n_plots > 1) {
      for (j in 1:(n_plots - 1)) {
        currentplot <- arranged.plots[[i]]$patches$plots[[j]]
        currentplot <- currentplot + theme(
          axis.title.y = element_blank(),
          strip.text.y.left = element_blank(),
          axis.text.y = element_blank(),
          strip.background = element_blank(),
          axis.ticks.y = element_blank(),
          line = element_blank()
        )

        arranged.plots[[i]]$patches$plots[[j]] <- currentplot
      }
      arranged.plots[[i]] <- arranged.plots[[i]] + theme(
        axis.title.y = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
    } else if (n_plots == 1) {
      arranged.plots[[i]] <- arranged.plots[[i]] + theme(
        axis.title.y = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
    }
  }

  ## move bigwig plot range to sub title
  if (!is.null(bigwig)) {
    for (i in seq_along(arranged.plots)) {
      if (bw_idx != 1) {
        bwplot <- arranged.plots[[i]]$patches$plots[[bw_idx]]
      } else {
        bwplot <- arranged.plots[[i]]
      }

      # get min max range for bigwig plot
      max_range <- max(bwplot@data$score)
      min_range <- min(bwplot@data$score)
      bw_range <- paste0(min_range, " - ", round(max_range))
      if (i == 1) {
        bw_range <- paste0(bw_range, " (score range)")
      }

      # add as sub title
      bwplot <- bwplot +
        labs(subtitle = bw_range) +
        theme(
          plot.subtitle = element_text(
            size = subtitle_size, hjust = 1
          ),
          axis.text.y = element_blank()
        )

      if (bw_idx != 1) {
        arranged.plots[[i]]$patches$plots[[bw_idx]] <- bwplot
      } else {
        arranged.plots[[i]] <- bwplot
      }
    }
  }
  
  ## make common legend for links
  if (!is.null(links)) {
    # get global score range
    max_score <- 0
    min_score <- 0
    for (i in seq_along(arranged.plots)) {
      scores <- arranged.plots[[i]]@data$score
      max_score <- max(max(scores), max_score)
      min_score <- min(min(scores), min_score)
    }

    # adjust color gradient
    for (i in seq_along(arranged.plots)) {
      arranged.plots[[i]] <- arranged.plots[[i]] +
        scale_color_gradient2(
          low = "red", mid = "grey", high = "blue",
          limits = c(min_score, max_score),
          n.breaks = 3)

      arranged.plots[[i]] <- arranged.plots[[i]] +
        theme(legend.position = "right")
    }
  }

  # adjust plot text, dimensions & add grid lines
  for (i in seq_along(arranged.plots)) {
    # remove chr position numbers
    arranged.plots[[i]] <- arranged.plots[[i]] + theme(
      axis.text.x = element_blank()
    )

    # change x axis text to chr start-end
    if (n_plots > 1) {
      startpos <- min(
        arranged.plots[[i]]$patches$plots[[1]]$data$position
      )
      endpos <- max(
        arranged.plots[[i]]$patches$plots[[1]]$data$position
      )
    } else if (n_plots == 1) {
      startpos <- min(arranged.plots[[i]]$data$position)
      endpos <- max(arranged.plots[[i]]$data$position)
    }

    x_label <- arranged.plots[[i]]@labels$x
    chrpos <- sub(" position \\(bp\\)", "", x_label)
    arranged.plots[[i]]@labels$x <- paste0(
      chrpos, "\n", startpos, "-\n", endpos
    )

    # change font size
    arranged.plots[[i]] <- arranged.plots[[i]] + theme(
      axis.title.x.bottom = element_text(
        size = coords_size, hjust = 0.5
      )
    )
  }

  # combine plots
  multi.plot <- wrap_plots(
    arranged.plots, ncol = length(arranged.plots)
  ) &
    theme(
      axis.line.y = element_line(
        linewidth = 1 / 2, color = "darkgrey"
      ),
      plot.margin = margin(2, 0, 2, 0, "pt")
    )

  if (!is.null(links)) {
    multi.plot <- (multi.plot | guide_area()) +
      plot_layout(guides = "collect", widths = c(15, 1)) & 
      theme(
        legend.position = "inside",
        legend.position.inside = c(1, 0),
        legend.justification = c(1, 0)
      )
  }

  return(multi.plot)
}


#' Plot data from BigWig files
#'
#' Create coverage tracks, heatmaps, or line plots from bigwig files.
#'
#' Note that this function does not work on windows.
#'
#' @param region GRanges object specifying region to plot
#' @param bigwig A named list or named character vector of bigwig file paths.
#' The name of each element will be displayed alongside the track in the final
#' plot.
#' @param smooth Number of bases to smooth data over (rolling mean). If NULL,
#' do not apply smoothing.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' @param type Plot type. Can be one of "line", "heatmap", or "coverage"
#' @param y_label Y-axis label
#' @param bigwig.scale Scaling to apply to data from different bigwig files.
#' Can be:
#' \itemize{
#' \item{common: plot each bigwig on a common scale (default)}
#' \item{separate: plot each bigwig on a separate scale ranging from zero to the
#' maximum value for that bigwig file within the plotted region}
#' }
#' @param ymax Maximum value for Y axis. Can be one of:
#'  - `NULL`: set to the highest value among all the tracks (default)
#'  - `qXX`: clip the maximum value to the XX quantile (for example, q95 will
#'  set the maximum value to 95% of the maximum value in the data). This can
#'  help remove the effect of extreme values that may otherwise distort the
#'  scale.
#'  - numeric: manually define a Y-axis limit
#' @param max.downsample Maximum amount of downsampling to apply. Corresponds to
#' the minimum number of positions to be kept when downsampling. The
#' downsampling rate is adaptive to the window size, but this parameter will set
#' the minimum possible number of positions to include so that plots do not
#' become too sparse when the window size is small.
#' @param downsample.rate Fraction of positions to retain when downsampling.
#' Retaining more positions can give a higher-resolution plot but can make the
#' number of points large, resulting in larger file sizes when saving the plot
#' and a longer period of time needed to draw the plot.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_tile xlab ylab geom_area scale_fill_viridis_c scale_color_grey scale_fill_grey facet_wrap
#' @importFrom RcppRoll roll_mean
#' @importFrom GenomicRanges start end seqnames width
#' @importFrom dplyr slice_sample group_by mutate ungroup
#' @concept visualization
#' @return Returns a ggplot object
#'
#' @export
BigwigTrack <- function(
  region,
  bigwig,
  smooth = 200,
  extend.upstream = 0,
  extend.downstream = 0,
  type = "coverage",
  y_label = "bigWig",
  bigwig.scale = "common",
  ymax = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1
) {
  if (!inherits(x = bigwig, what = "list")) {
    # accept a (named) character vector of file paths
    bigwig <- as.list(x = bigwig)
  }
  if (is.null(x = names(x = bigwig))) {
    names(x = bigwig) <- paste0("bigWig_", seq_along(along.with = bigwig))
  }
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
    message(
      "Please install rtracklayer. ",
      "http://www.bioconductor.org/packages/rtracklayer/"
    )
    return(NULL)
  }
  region <- FindRegion(
    object = NULL,
    region = region,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  all.data <- data.frame()
  for (i in seq_along(bigwig)) {
    region_data <- rtracklayer::import(
      con = bigwig[[i]],
      which = region,
      as = "NumericList"
    )[[1]]
    if (!is.null(x = smooth)) {
      region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
    }
    region_data <- data.frame(
      position = start(x = region):end(x = region),
      score = region_data,
      stringsAsFactors = FALSE,
      bw = names(x = bigwig)[[i]]
    )
    if (bigwig.scale == "separate") {
      # scale to fraction of max for each separately
      file.max <- max(region_data$score, na.rm = TRUE)
      region_data$score <- region_data$score / file.max
    }
    all.data <- rbind(all.data, region_data)
  }
  all.data$bw <- factor(x = all.data$bw, levels = names(x = bigwig))
  window.size <- width(x = region)
  sampling <- min(max.downsample, window.size * downsample.rate)
  coverages <- slice_sample(.data = all.data, n = as.integer(x = sampling))

  covmax <- signif(x = max(coverages$score, na.rm = TRUE), digits = 2)
  if (is.null(x = ymax)) {
    ymax <- covmax
  } else if (is.character(x = ymax)) {
    if (!startsWith(x = ymax, prefix = "q")) {
      stop("Unknown ymax requested. Must be NULL, a numeric value, or
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
    }
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
    ) / 100
    ymax <- covmax * percentile.use
  }

  # perform clipping
  coverages$score[coverages$score > ymax] <- ymax

  if (type == "line") {
    p <- ggplot(
      data = coverages,
      mapping = aes(x = .data[["position"]], y = .data[["score"]], color = .data[["bw"]])
    ) +
      geom_line() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1) +
      scale_color_grey()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    all.data$bin <- floor(x = all.data$position / smooth)
    all.data <- group_by(all.data, bin, bw)
    all.data <- mutate(all.data, score = mean(x = score))
    all.data <- ungroup(all.data)
    all.data <- unique(x = all.data[, c("bin", "score", "bw")])
    p <- ggplot(
      data = all.data,
      mapping = aes(x = .data[["bin"]], y = 1, fill = .data[["score"]])
    ) +
      geom_tile() +
      scale_fill_viridis_c() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1)
  } else if (type == "coverage") {
    p <- ggplot(
      data = coverages,
      mapping = aes(x = .data[["position"]], y = .data[["score"]], fill = .data[["bw"]])
    ) +
      geom_area() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1) +
      scale_fill_grey()
  }
  chromosome <- as.character(x = seqnames(x = region))
  p <- p + theme_browser(axis.text.y = TRUE) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label)
  return(p)
}

#' Plot sequencing depth correlation
#'
#' Compute the correlation between total counts and each reduced
#' dimension component.
#'
#' @param object A [SeuratObject::Seurat()] object
#' @param reduction Name of a dimension reduction stored in the
#' input object
#' @param assay Name of assay to use for sequencing depth. If NULL, use the
#' default assay.
#' @param n Number of components to use. If `NULL`, use all components.
#' @param ... Additional arguments passed to [stats::cor()]
#' @return Returns a [ggplot2::ggplot()] object
#' @export
#' @importFrom SeuratObject Embeddings DefaultAssay
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous ylab ylim theme_light ggtitle
#' @importFrom stats cor
#' @concept visualization
#' @examples
#' \donttest{
#' DepthCor(object = atac_small)
#' }
DepthCor <- function(object, assay = NULL, reduction = "lsi", n = 10, ...) {
  assay <- assay %||% DefaultAssay(object = object)
  dr <- object[[reduction]]
  embed <- Embeddings(object = dr)
  counts <- object[[paste0("nCount_", assay)]]
  embed <- embed[rownames(x = counts), ]
  n <- n %||% ncol(x = embed)
  embed <- embed[, seq_len(length.out = n)]
  depth.cor <- as.data.frame(cor(x = embed, y = counts, ...))
  depth.cor$counts <- depth.cor[, 1]
  depth.cor$Component <- seq_len(length.out = nrow(x = depth.cor))
  p <- ggplot(depth.cor, aes(x = .data[["Component"]], y = .data[["counts"]])) +
    geom_point() +
    scale_x_continuous(n.breaks = n, limits = c(1, n)) +
    ylab("Correlation") +
    ylim(c(-1, 1)) +
    theme_light() +
    ggtitle("Correlation between depth and reduced dimension components",
      subtitle = paste0("Assay: ", assay, "\t", "Reduction: ", reduction)
    )
  return(p)
}

# Get density of points in 2 dimensions.
#
# Modified from original code by Kamil Slowikowski
# https://slowkow.com/notes/ggplot2-color-by-density/
#
# @param x A numeric vector.
# @param y A numeric vector.
# @param n_sub Number of points to sample
# @return The density within each square.
get_density <- function(x, y, n_sub = 50000, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install MASS: install.packages('MASS')")
  }
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Please install fields: install.packages('fields')")
  }

  n <- length(x)
  if (n > n_sub) {
    idx <- sample.int(n, n_sub)
    xs <- x[idx]
    ys <- y[idx]
  } else {
    xs <- x
    ys <- y
  }

  # KDE on subsample
  dens <- MASS::kde2d(xs, ys, n = 1000, ...)

  # Interpolate back onto full data
  out <- fields::interp.surface(
    obj = list(x = dens$x, y = dens$y, z = dens$z),
    loc = cbind(x, y)
  )
  return(out)
}

#' Create GWAS locus zoom track
#'
#' This function will plot the p-values associated with variants in a given
#' region of the genome (genome position: x-axis; -log10(p): y-axis).
#'
#' If an LD file is provided using the `ld.file` parameter, or a column
#' named `r2` is present in the input GWAS data, the variants will be
#' colored according to their LD value.
#'
#' If a fine mapping file is provided using the `credset.file` parameter,
#' or a column named `in_credset` is present in the input GWAS data,
#' variants in the credible set will be denoted by shape. If LD information is
#' not provided, credible variants will also be denoted by color.
#'
#' @param region Genomic region ([GenomicRanges::GRanges] or a string that can
#' be converted to `GRanges` like "chr10:112900000-113100000")
#' @param gwas Path to GWAS summary statistics file, or a dataframe
#' containing the gwas data in the GWAS-SSF format
#' @param ld.file Path to LD file. Optional.
#' @param ld.lead.snp Lead SNP for LD (required if ld.file provided)
#' @param credset.file Path to credible set file. Optional.
#' @param credset.threshold PIP threshold for credible sets (default: 0.01)
#' @param p.threshold Genome-wide significance threshold (default: 5e-8)
#' @param ymax Maximum y-axis value (default: auto)
#' @param point.size Point size (default: 1)
#' @param point.color Point color when no LD (default: "steelblue")
#' @param show.axis Show x-axis (default: TRUE)
#' @return ggplot2 object
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_hline theme_classic labs theme element_blank element_line element_text scale_shape_manual scale_size_manual scale_color_manual scale_y_continuous
#' @importFrom Seqinfo seqnames
#' @importFrom GenomicRanges start end
#' @export
#' @concept visualization
GWASTrack <- function(
  gwas,
  region,
  ld.file = NULL,
  ld.lead.snp = NULL,
  credset.file = NULL,
  credset.threshold = 0.01,
  p.threshold = 5e-8,
  ymax = NULL,
  point.size = 1,
  point.color = "steelblue",
  show.axis = TRUE
) {
  # Load GWAS data
  if (is.character(x = gwas)) {
    gwas <- LoadGWAS(gwas.file = gwas)
  }

  if (!inherits(x = region, what = "GRanges")) {
    region <- GRanges(region)
  }

  # subset to region. Compare chromosome names with the "chr" prefix removed
  # from both sides, so that a GWAS file naming chromosomes "1" matches an
  # object naming them "chr1"
  chromosome <- as.character(x = seqnames(x = region))
  gwas[["chrom.key"]] <- NormalizeChromosome(x = gwas[["chromosome"]])
  chromosomes.present <- unique(x = gwas[["chromosome"]])
  on.chromosome <- gwas[["chrom.key"]] == NormalizeChromosome(x = chromosome)
  gwas <- gwas[
    on.chromosome &
      gwas[["base_pair_location"]] >= start(x = region) &
      gwas[["base_pair_location"]] <= end(x = region),
  ]
  gwas <- gwas[!is.na(x = gwas[["p_value"]]), ]

  if (nrow(x = gwas) == 0) {
    if (!any(on.chromosome)) {
      stop(
        "No GWAS data found on ", chromosome,
        ". Chromosomes present in the GWAS data: ",
        paste(head(x = chromosomes.present), collapse = ", ")
      )
    }
    stop("No GWAS data found in region")
  }
  gwas[["log10p"]] <- -log10(x = gwas[["p_value"]])

  # Validate LD parameters
  if (!is.null(x = ld.file) && is.null(x = ld.lead.snp)) {
    stop("ld.lead.snp required when ld.file provided")
  }

  # LocusZoom colors
  ld_colors <- c(
    "r2_0-0.2" = "#0000CD",
    "r2_0.2-0.4" = "#00CED1",
    "r2_0.4-0.6" = "#32CD32",
    "r2_0.6-0.8" = "#FFA500",
    "r2_0.8-1.0" = "#FF0000"
  )

  # Merge LD data
  if (!is.null(x = ld.file)) {
    ld_data <- LoadLDData(ld.file = ld.file)
    ld_data[["chrom.key"]] <- NormalizeChromosome(
      x = ld_data[["chromosome"]]
    )
    ld_data[["chromosome"]] <- NULL
    gwas <- merge(
      x = gwas,
      y = ld_data,
      by = c("chrom.key", "base_pair_location"),
      all.x = TRUE
    )
    gwas[["r2"]] <- as.numeric(x = gwas[["r2"]])
  }

  if ("r2" %in% colnames(x = gwas)) {
    gwas[["ld_category"]] <- cut(
      gwas[["r2"]],
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c(
        "r2_0-0.2", "r2_0.2-0.4",
        "r2_0.4-0.6", "r2_0.6-0.8", "r2_0.8-1.0"
      ),
      include.lowest = TRUE
    )
  }

  # Merge credible set data
  if (!is.null(x = credset.file)) {
    credset_data <- LoadCredibleSets(
      credset.file = credset.file, credset.threshold = credset.threshold
    )
    credset_data[["chrom.key"]] <- NormalizeChromosome(
      x = credset_data[["chromosome"]]
    )
    credset_data[["chromosome"]] <- NULL
    gwas <- merge(
      x = gwas,
      y = credset_data,
      by = c("chrom.key", "base_pair_location"),
      all.x = TRUE
    )
    gwas[["in_credset"]] <- !is.na(x = gwas[["pip"]])
  }

  # the join key is internal, and should not appear in the returned plot data
  gwas[["chrom.key"]] <- NULL

  # Y-axis limit
  if (is.null(x = ymax)) {
    ymax <- max(gwas[["log10p"]], na.rm = TRUE) * 1.1
  }

  # Build plot
  if ("in_credset" %in% colnames(x = gwas)) {
    if ("ld_category" %in% colnames(x = gwas)) {
      # LD + credible sets
      p <- ggplot(data = gwas, mapping = aes(
        x = .data[["base_pair_location"]],
        y = .data[["log10p"]],
        color = .data[["ld_category"]]
      )) +
        geom_point(
          aes(shape = .data[["in_credset"]], size = .data[["in_credset"]]),
          alpha = 0.6
        ) +
        scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 18),
          guide = "none"
        ) +
        scale_size_manual(
          values = c("FALSE" = point.size, "TRUE" = point.size * 2),
          guide = "none"
        ) +
        scale_color_manual(
          values = ld_colors,
          name = expression(LD ~ (r^2)),
          na.value = "grey50"
        )
    } else {
      # Credible sets only
      p <- ggplot(data = gwas, mapping = aes(
        x = .data[["base_pair_location"]],
        y = .data[["log10p"]]
      )) +
        geom_point(
          mapping = aes(
            shape = .data[["in_credset"]],
            size = .data[["in_credset"]],
            color = .data[["in_credset"]]
          ),
          alpha = 0.6
        ) +
        scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 18),
          guide = "none"
        ) +
        scale_size_manual(
          values = c("FALSE" = point.size, "TRUE" = point.size * 2),
          guide = "none"
        ) +
        scale_color_manual(
          values = c("FALSE" = point.color, "TRUE" = "#FF0000"),
          name = "Credible set",
          labels = c("FALSE" = "No", "TRUE" = "Yes")
        )
    }
  } else if ("ld_category" %in% colnames(x = gwas)) {
    # LD only
    p <- ggplot(
      data = gwas, mapping = aes(
        x = .data[["base_pair_location"]],
        y = .data[["log10p"]],
        color = .data[["ld_category"]]
      )
    ) +
      geom_point(size = point.size, alpha = 0.6) +
      scale_color_manual(
        values = ld_colors,
        name = expression(LD ~ (r^2)),
        na.value = "grey50"
      )
  } else {
    # Basic plot
    p <- ggplot(
      data = gwas, mapping = aes(
        x = .data[["base_pair_location"]],
        y = .data[["log10p"]]
      )
    ) +
      geom_point(color = point.color, size = point.size, alpha = 0.6)
  }

  # Common elements
  p <- p +
    geom_hline(
      yintercept = -log10(x = p.threshold),
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    theme_classic() +
    labs(
      x = paste0(chromosome, " position (bp)"),
      y = expression(-log[10](italic(P)))
    ) +
    theme(
      axis.title.x = if (show.axis) element_text() else element_blank(),
      axis.text.x = if (show.axis) element_text() else element_blank(),
      axis.line.x = if (show.axis) element_line() else element_blank(),
      axis.ticks.x = if (show.axis) element_line() else element_blank()
    ) +
    xlim(c(start(x = region), end(x = region)))

  return(p)
}

globalVariables(".data")
#' Scatterplot colored by point density
#'
#' Create a scatterplot using variables in the object metadata
#' and color cells by the density of points in the x-y space.
#'
#' @param object A Seurat object
#' @param x Name of metadata variable to plot on x axis
#' @param y Name of metadata variable to plot on y axis
#' @param log_x log10 transform x values
#' @param log_y log10 transform y values
#' @param quantiles Vector of quantiles to display
#' for x and y data distribution. Must be integer values
#' between 0 and 100.
#' TRUE can be passed as a shorthand way to set
#' `c(5, 10, 90, 95)`. If FALSE or NULL, no quantile
#' information is displayed
#' @param raster Convert points to raster format. If NULL, points will
#' automatically be rasterized if plotting more than 100,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots.
#' @return Returns a ggplot object
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c theme_bw scale_x_log10 scale_y_log10 geom_vline geom_hline labs
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @export
#' @concept visualization
DensityScatter <- function(
  object,
  x,
  y,
  log_x = FALSE,
  log_y = FALSE,
  quantiles = NULL,
  raster = NULL,
  raster.dpi = c(512, 512)
) {
  md <- object[[]]
  if (!(x %in% colnames(x = md))) {
    stop(x, " not found")
  }
  if (!(y %in% colnames(x = md))) {
    stop(y, " not found")
  }
  log10p <- function(x) {
    return(log10(x = x + 1))
  }
  null_fn <- function(x) {
    return(x)
  }
  logfnx <- ifelse(test = log_x, yes = log10p, no = null_fn)
  logfny <- ifelse(test = log_y, yes = log10p, no = null_fn)

  md$Density <- get_density(
    x = logfnx(md[[x]]),
    y = logfny(md[[y]]),
    h = c(1, 1)
  )
  md <- md[order(md$Density), ]

  # resolve the automatic case first, so raster is always TRUE or FALSE below
  raster <- raster %||% (nrow(x = md) > 100000)

  if (isTRUE(x = raster) &&
      !requireNamespace(package = "scattermore", quietly = TRUE)) {
    warning("scattermore is not installed, plot cannot be rasterized")
    raster <- FALSE
  }

  # quantiles
  use_quantile <- FALSE
  if (!is.null(x = quantiles)) {
    use_quantile <- TRUE
    # set default if TRUE passed
    if (is.logical(x = quantiles)) {
      if (quantiles) {
        quantiles <- c(5, 10, 90, 95)
      } else {
        use_quantile <- FALSE
      }
    }
    # make sure integers between 0 and 100
    if (!(all(quantiles >= 0) && all(quantiles <= 100))) {
      warning("Quantile values must be between 0 and 100",
        immediate. = TRUE
      )
      use_quantile <- FALSE
    }
    if (any(sapply(X = quantiles, FUN = function(x) x %% 1 != 0))) {
      warning("Quantile values must be integers",
        immediate. = TRUE
      )
      use_quantile <- FALSE
    }
  }
  if (use_quantile) {
    # convert to string
    quantiles <- paste0(quantiles, "%")

    # quantiles for x and y
    x_quant <- quantile(x = md[[x]], probs = seq(0, 1, 0.01))
    y_quant <- quantile(x = md[[y]], probs = seq(0, 1, 0.01))
    xlines <- x_quant[quantiles]
    ylines <- y_quant[quantiles]

    # round
    xlines <- round(x = xlines, digits = 2)
    ylines <- round(x = ylines, digits = 2)
  }
  p <- ggplot(
    data = md,
    mapping = aes(x = .data[[x]], y = .data[[y]], color = .data[["Density"]])
  )
  if (isTRUE(x = raster)) {
    p <- p + scattermore::geom_scattermore(pixels = raster.dpi, pointsize = 3.2)
  } else {
    p <- p + geom_point(size = 1)
  }
  p <- p + scale_color_viridis_c(option = "B") + theme_bw()
  if (log_x) {
    p <- p + scale_x_log10()
  }
  if (log_y) {
    p <- p + scale_y_log10()
  }
  if (use_quantile) {
    p <- p +
      geom_vline(xintercept = unname(obj = xlines), color = "red") +
      geom_hline(yintercept = unname(obj = ylines), color = "red") +
      labs(
        title = "Quantiles",
        subtitle = paste0(
          x, ": ", paste0(names(x = xlines), ":",
            xlines,
            collapse = " "
          ), "\n",
          y, ": ", paste0(names(x = ylines), ":",
            ylines,
            collapse = " "
          )
        )
      )
  }
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
#' @param split.by A metadata variable to split the plot by. For example,
#' grouping by "celltype" and splitting by "batch" will create separate plots
#' for each celltype and batch.
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
#' `label.top` will be ignored.
#' @export
#' @concept visualization
#' @concept footprinting
#' @importFrom SeuratObject DefaultAssay
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap xlab ylab theme element_blank geom_label guides guide_legend theme_classic
#' @importFrom dplyr group_by summarize top_n
#' @import patchwork
PlotFootprint <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  split.by = NULL,
  idents = NULL,
  label = TRUE,
  repel = TRUE,
  show.expected = TRUE,
  normalization = "subtract",
  label.top = 3,
  label.idents = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  splitby_str <- "__signac_tmp__"
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("The requested assay is not a ChromatinAssay5.")
  }
  if (!is.null(x = split.by)) {
    if (is.null(x = group.by)) {
      grouping.var <- Idents(object = object)
    } else {
      grouping.var <- object[[group.by]][, 1]
    }
    # combine split.by and group.by information
    combined.var <- paste0(object[[split.by]][, 1], splitby_str, grouping.var)
    object$grouping_tmp <- combined.var
    group.by <- "grouping_tmp"
    if (!is.null(x = idents)) {
      # adjust idents parameter with new split.by information
      idents.keep <- combined.var[grouping.var %in% idents]
      idents <- unique(x = idents.keep)
    }
  }
  # TODO add option to show variance among cells
  plot.data <- GetFootprintData(
    object = object,
    features = features,
    assay = assay,
    group.by = group.by,
    idents = idents
  )
  if (length(x = plot.data) == 1) {
    stop("Footprinting data not found")
  }
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
  flank.pos <- abs(x = obs$position)
  flank.size <- unname(obj = base[obs$feature])
  obs$flanks <- (flank.pos > flank.size) & (flank.pos < (flank.size + 50))

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

  # split back into group.by and split.by
  if (!is.null(x = split.by)) {
    splitvar <- strsplit(x = obs[["group"]], split = splitby_str)
    obs[["group"]] <- sapply(X = splitvar, FUN = `[[`, 2)
    obs[["split"]] <- sapply(X = splitvar, FUN = `[[`, 1)
  }

  # find flanking accessibility for each group and each feature
  flanks <- obs[obs$flanks, ]
  flanks <- group_by(.data = flanks, feature, group)
  flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))

  # find top n groups for each feature
  topmean <- top_n(x = flankmeans, n = label.top, wt = mn)

  # find the top for each feature to determine axis limits
  ymax <- top_n(x = flankmeans, n = 1, wt = mn)
  ymin <- top_n(x = flankmeans, n = 1, wt = -mn)

  # make df for labels. Labels sit at position 75 by default, but the plotted
  # region can be narrower than that, so use the closest position available
  all.positions <- unique(x = obs$position)
  label.position <- all.positions[which.min(x = abs(x = all.positions - 75))]
  label.df <- data.frame()
  sub <- obs[obs$position == label.position, ]
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
        (sub$group %in% groups.use),
    ]
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
        x = .data[["position"]],
        y = .data[["norm.value"]],
        color = .data[["group"]],
        label = .data[["label"]]
      )
    )
    p <- p +
      geom_line(linewidth = 0.2) +
      xlab("Distance from motif") +
      ylab(label = "Tn5 insertion\nenrichment") +
      theme_classic() +
      ggtitle(label = features[[i]]) +
      ylim(c(axis.min, axis.max)) +
      guides(color = guide_legend(override.aes = list(linewidth = 1)))
    if (!is.null(x = split.by)) {
      p <- p + facet_wrap(facets = ~split)
    }
    if (label) {
      if (repel) {
        if (!requireNamespace(package = "ggrepel", quietly = TRUE)) {
          warning(
            "Please install ggrepel to enable repel=TRUE: ",
            "install.packages('ggrepel')"
          )
          p <- p + geom_label(show.legend = FALSE)
        } else {
          p <- p + ggrepel::geom_label_repel(
            box.padding = 0.5, show.legend = FALSE
          )
        }
      } else {
        p <- p + geom_label(show.legend = FALSE)
      }
    }
    if (show.expected) {
      if (!is.null(x = split.by)) {
        warning("Cannot plot expected enrichment with split.by")
      } else {
        df <- expect[expect$feature == features[[i]], ]
        p1 <- ggplot(
          data = df,
          mapping = aes(x = .data[["position"]], y = .data[["norm.value"]])
        ) +
          geom_line(linewidth = 0.2) +
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
      }
    }
    plotlist[[i]] <- p
  }
  plots <- wrap_plots(plotlist)
  return(plots)
}

globalVariables(
  names = c("group"),
  package = "Signac"
)
#' Region heatmap
#'
#' Plot fragment counts within a set of regions.
#'
#' @param object The output of [RegionMatrix()]: a list containing two elements:
#' - `matrix`: a named list of region by position matrices, one for each group
#' of cells, with the name of each element corresponding to the group identity.
#' - `parameters`: a list of function parameters "upstream", "downstream" and
#' "cells".
#' Optionally, a list of such lists can be supplied for multi-assay plotting. In
#' this case, the name of each element should correspond to the assay name.
#' @param window Smoothing window to apply
#' @param normalize Normalize by number of cells in each group
#' @param order Order regions by the total number of fragments in the region
#' across all included identities
#' @param upstream Number of bases to include upstream of region. If NULL, use
#' all bases that were included in the `RegionMatrix` function call. Note
#' that this value cannot be larger than the value for `upstream` given in
#' the original `RegionMatrix` function call. If NULL, use parameters that
#' were given in the `RegionMatrix` function call
#' @param downstream Number of bases to include downstream of region. See
#' documentation for `upstream`
#' @param max.cutoff Maximum cutoff value. Data above this value will be clipped
#' to the maximum value. A quantile maximum can be specified in the form of
#' "q##" where "##" is the quantile (eg, "q90" for 90th quantile). If NULL, no
#' cutoff will be set
#' @param cols Vector of colors to use as the maximum value of the color scale.
#' One color must be supplied for each assay. If NULL, the default ggplot2
#' colors are used.
#' @param min.counts Minimum total counts to display region in plot
#' @param idents Cell identities to include. Note that cells cannot be
#' regrouped, this will require re-running `RegionMatrix` to generate a
#' new set of matrices
#' @param group.order Order of groups to be shown in the plot. This should be a
#' character vector. If NULL, the group order will not be changed.
#' @param nrow Number of rows to use when creating plot. If NULL, chosen
#' automatically by ggplot2
#'
#' @seealso RegionMatrix
#'
#' @return Returns a ggplot2 object
#'
#' @importFrom RcppRoll roll_sum
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes facet_wrap geom_raster guides theme .data vars element_blank element_text scale_fill_gradient ylab guide_legend xlab
#' @importFrom scales hue_pal
#' @importFrom patchwork wrap_plots
#'
#' @export
#' @concept visualization
#' @concept heatmap
RegionHeatmap <- function(
  object,
  idents = NULL,
  group.order = NULL,
  normalize = TRUE,
  upstream = 3000,
  downstream = 3000,
  max.cutoff = "q95",
  cols = NULL,
  min.counts = 1,
  window = (upstream + downstream) / 30,
  order = TRUE,
  nrow = NULL
) {
  # for multiassay support
  if ("parameters" %in% names(x = object)) {
    object <- list("default" = object)
  }
  assay <- names(x = object)

  if (is.null(x = cols)) {
    colors_all <- hue_pal()(length(x = assay))
    names(x = colors_all) <- assay
  } else {
    if (length(x = cols) != length(x = assay)) {
      stop("Wrong number of colors supplied. Must give one color per assay")
    }
    colors_all <- cols
    if (is.null(x = names(x = colors_all))) {
      names(x = colors_all) <- assay
    } else if (!all(names(x = colors_all) %in% assay)) {
      names(x = colors_all) <- assay
    }
  }

  all.assay <- data.frame()
  for (j in seq_along(along.with = object)) {
    upstream.max <- object[[j]]$parameters$upstream
    downstream.max <- object[[j]]$parameters$downstream
    matlist <- object[[j]]$matrix
    cells.per.group <- object[[j]]$parameters$cells

    upstream <- upstream %||% upstream.max
    downstream <- downstream %||% downstream.max

    # define clipping
    cols.keep <- (upstream.max - upstream + 1):(upstream.max + downstream + 1)

    if (!is.null(x = idents)) {
      valid.idents <- intersect(x = idents, y = names(x = matlist))
      matlist <- matlist[valid.idents]
    }
    if (length(x = matlist) == 0) {
      stop("None of the requested idents found")
    }

    if (j == 1) {
      rsums <- lapply(X = matlist, FUN = rowSums)
      rsums <- Reduce(f = `+`, x = rsums)
      rows.retain <- rsums >= min.counts
    }

    # remove low count
    matlist <- lapply(X = matlist, FUN = function(x) {
      x[rows.retain, ]
    })

    if (order && j == 1) {
      rsums <- lapply(X = matlist, FUN = rowSums)
      rsums <- Reduce(f = `+`, x = rsums)
      order.use <- base::order(rsums)
    }

    for (i in seq_along(along.with = matlist)) {
      grp.name <- names(x = matlist)[[i]]
      m <- matlist[[i]]

      # clip up/downstream
      m <- m[, cols.keep]
      colnames(m) <- seq_len(length.out = ncol(x = m))

      if (order) {
        m <- m[order.use, ]
      }

      if (normalize) {
        m <- m / cells.per.group[[grp.name]]
        guide.label <- "Fragment counts\nper cell"
      } else {
        guide.label <- "Fragment\ncount"
      }

      smoothed <- apply(
        X = m,
        MARGIN = 1,
        FUN = roll_sum,
        n = window,
        by = window
      )
      # create dataframe
      smoothed <- as.data.frame(x = smoothed)
      colnames(smoothed) <- seq_len(length.out = ncol(x = smoothed))

      # values are clipped once across all groups after the loop, so that the
      # shared colour scale corresponds to a single threshold

      # add extra column as bin ID
      regions <- colnames(x = smoothed)
      smoothed$bin <- seq_len(length.out = nrow(x = smoothed))
      smoothed <- pivot_longer(
        data = smoothed,
        cols = all_of(regions)
      )
      smoothed$group <- grp.name
      if (i == 1) {
        df <- smoothed
      } else {
        df <- rbind(df, smoothed)
      }
    }

    # fix bin label
    df$bin <- (df$bin - (upstream / window)) * window
    df$name <- as.numeric(x = df$name)
    df$assay <- assay[[j]]

    all.assay <- rbind(all.assay, df)
  }

  if (!is.null(x = max.cutoff) && !is.na(x = max.cutoff)) {
    cutoff <- SetQuantile(cutoff = max.cutoff, data = all.assay$value)
    all.assay$value[all.assay$value > cutoff] <- cutoff
  }

  maxval <- max(all.assay$value)

  if (!is.null(x = group.order)) {
    if (length(x = group.order) != length(x = unique(x = all.assay$group))) {
      warning(
        "Incorrect number of groups provided in group.order parameter.",
        " Groups will not be reordered"
      )
    } else {
      all.assay$group <- factor(x = all.assay$group, levels = group.order)
    }
  }

  # create separate heatmap for each assay so that color scales are different
  plist <- list()
  for (i in seq_along(along.with = assay)) {
    data.use <- all.assay[all.assay$assay == assay[[i]], ]
    pp <- ggplot(
      data = data.use,
      mapping = aes(
        x = .data[["bin"]],
        y = .data[["name"]],
        fill = .data[["value"]]
      )
    ) +
      facet_wrap(
        facets = vars(.data[["group"]]),
        scales = "free_y",
        nrow = nrow
      ) +
      geom_raster() +
      theme_browser(legend = TRUE) +
      ylab("Region") +
      ggtitle(assay[[i]]) +
      scale_fill_gradient(
        low = "white",
        high = colors_all[[assay[[i]]]],
        limits = c(0, maxval)
      ) +
      guides(
        fill = guide_legend(
          title = ifelse(
            test = length(x = assay) > 1,
            yes = assay[[i]],
            no = guide.label
          ),
          keywidth = 1 / 2,
          keyheight = 1
        )
      ) +
      theme(
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank()
      )
    plist[[i]] <- pp
  }
  p <- wrap_plots(plist, guides = "collect")
  return(p)
}

#' Region plot
#'
#' Plot fragment counts within a set of regions.
#'
#' @param object The output of [RegionMatrix()]: a list containing two elements:
#' - `matrix`: a named list of region by position matrices, one for each group
#' of cells, with the name of each element corresponding to the group identity.
#' - `parameters`: a list of function parameters "upstream", "downstream" and
#' "cells".
#' Optionally, a list of such lists can be supplied for multi-assay plotting. In
#' this case, the name of each element should correspond to the assay name.
#' @param window Smoothing window to apply
#' @param normalize Normalize by number of cells in each group
#' @param upstream Number of bases to include upstream of region. If NULL, use
#' all bases that were included in the `RegionMatrix` function call. Note
#' that this value cannot be larger than the value for `upstream` given in
#' the original `RegionMatrix` function call. If NULL, use parameters that
#' were given in the `RegionMatrix` function call
#' @param downstream Number of bases to include downstream of region. See
#' documentation for `upstream`
#' @param idents Cell identities to include. Note that cells cannot be
#' regrouped, this will require re-running [RegionMatrix()] to generate a
#' new set of matrices
#' @param group.order Order of groups to be shown in the plot. This should be a
#' character vector. If NULL, the group order will not be changed.
#' @param nrow Number of rows to use when creating plot. If NULL, chosen
#' automatically by ggplot2
#'
#' @seealso RegionMatrix
#'
#' @return Returns a ggplot2 object
#'
#' @importFrom RcppRoll roll_sum
#' @importFrom tidyselect all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes facet_wrap guides theme theme_classic element_blank element_text ylab xlab geom_line .data vars
#'
#' @export
#' @concept visualization
#' @concept heatmap
RegionPlot <- function(
  object,
  idents = NULL,
  group.order = NULL,
  normalize = TRUE,
  upstream = NULL,
  downstream = NULL,
  window = (upstream + downstream) / 500,
  nrow = NULL
) {
  # for multiassay support
  if ("parameters" %in% names(x = object)) {
    object <- list("default" = object)
  }
  assay <- names(x = object)

  all.assay <- data.frame()
  for (j in seq_along(along.with = object)) {
    upstream.max <- object[[j]]$parameters$upstream
    downstream.max <- object[[j]]$parameters$downstream
    matlist <- object[[j]]$matrix
    cells.per.group <- object[[j]]$parameters$cells

    upstream <- upstream %||% upstream.max
    downstream <- downstream %||% downstream.max

    # define clipping
    cols.keep <- (upstream.max - upstream + 1):(upstream.max + downstream + 1)

    if (!is.null(x = idents)) {
      valid.idents <- intersect(x = idents, y = names(x = matlist))
      matlist <- matlist[valid.idents]
    }
    if (length(x = matlist) == 0) {
      stop("None of the requested idents found")
    }

    for (i in seq_along(along.with = matlist)) {
      grp.name <- names(x = matlist)[[i]]
      m <- matlist[[i]]

      # clip up/downstream
      m <- m[, cols.keep]
      colnames(m) <- seq_len(length.out = ncol(x = m))

      if (normalize) {
        m <- m / cells.per.group[[grp.name]]
        guide.label <- "Fragment counts\nper cell"
      } else {
        guide.label <- "Fragment\ncount"
      }

      totals <- colSums(x = m)
      smoothed <- roll_sum(x = totals, n = window, by = window)
      grp.name <- names(x = matlist)[[i]]

      if (i == 1) {
        df <- data.frame(
          "data" = smoothed,
          "group" = grp.name,
          "bin" = seq_along(along.with = smoothed)
        )
      } else {
        df <- rbind(
          df,
          data.frame(
            "data" = smoothed,
            "group" = grp.name,
            "bin" = seq_along(along.with = smoothed)
          )
        )
      }
    }

    # fix bin label
    df$bin <- (df$bin - (upstream / window)) * window

    df$assay <- assay[[j]]
    all.assay <- rbind(all.assay, df)
  }

  if (!is.null(x = group.order)) {
    if (length(x = group.order) != length(x = unique(x = all.assay$group))) {
      warning(
        "Incorrect number of groups provided in group.order parameter.",
        " Groups will not be reordered"
      )
    } else {
      all.assay$group <- factor(x = all.assay$group, levels = group.order)
    }
  }

  p <- ggplot(
    data = all.assay,
    aes(x = .data[["bin"]], y = .data[["data"]], color = .data[["assay"]])
  ) +
    facet_wrap(facets = vars(.data[["group"]])) +
    geom_line() +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    ylab(label = guide.label) +
    xlab("Distance from center (bp)")

  return(p)
}

#' @importFrom ggplot2 ylab scale_fill_manual unit element_text theme
#' @importMethodsFrom GenomicRanges start end
#' @importFrom SeuratObject WhichCells Idents DefaultAssay Idents<-
SingleCoveragePlot <- function(
  object,
  region,
  features = NULL,
  assay = NULL,
  split.assays = FALSE,
  assay.scale = "common",
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
  links = "linkpeaks",
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  bigwig.scale = "common",
  group.by = NULL,
  split.by = NULL,
  window = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  ymax = NULL,
  scale.factor = NULL,
  cells = NULL,
  idents = NULL,
  heights = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1,
  gwas = NULL,
  gwas.ld.file = NULL,
  gwas.ld.lead.snp = NULL,
  gwas.credset.file = NULL,
  gwas.credset.threshold = 0.01,
  variants = NULL,
  avi = NULL,
  avi.features = c("Expression", "Accessibility", "Total"),
  avi.type = "stacked",
  avi.bins = 200,
  avi.api.key = NULL,
  avi.max.workers = 4,
  avi.verbose = FALSE
) {
  valid.assay.scale <- c("common", "separate")
  if (!(assay.scale %in% valid.assay.scale)) {
    stop(
      "Unknown assay.scale requested. Please choose from: ",
      paste(valid.assay.scale, collapse = ", ")
    )
  }
  cells <- cells %||% Cells(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = assay, what = "list")) {
    assay <- list(assay)
  }
  lapply(X = assay, FUN = function(x) {
    if (!inherits(x = object[[x]], what = "ChromatinAssay5")) {
      stop("Requested assay is not a ChromatinAssay5.")
    }
  })
  is.granges <- inherits(x = object[[assay[[1]]]], what = "GRangesAssay")
  if (length(colnames(object)) > length(colnames(object[[assay[[1]]]]))) {
    object <- UpdateChromatinObject(
      object = object,
      chromatin.assay = assay,
      expression.assay = expression.assay,
      features = features
    )
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
    assay = assay[[1]],
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  if (!is.null(x = split.by)) {
    # combine split.by and group.by information
    grouping.var <- Idents(object = object)
    combined.var <- paste0(object[[split.by]][, 1], "_", grouping.var)
    object$grouping_tmp <- combined.var
    Idents(object = object) <- "grouping_tmp"
    group.by <- "grouping_tmp"
    if (!is.null(x = idents)) {
      # adjust idents parameter with new split.by information
      idents.keep <- combined.var[grouping.var %in% idents]
      idents <- unique(x = idents.keep)
    }
  }
  cells.per.group <- CellsPerGroup(
    object = object,
    group.by = group.by
  )

  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )

  # subset to used cells
  obj.groups <- obj.groups[cells]

  cm.list <- list()
  sf.list <- list()
  gsf.list <- list()
  for (i in seq_along(along.with = assay)) {
    reads.per.group <- AverageCounts(
      object = object,
      assay = assay[[i]],
      group.by = group.by,
      verbose = FALSE
    )
    cutmat <- CutMatrix(
      object = object,
      region = region,
      assay = assay[[i]],
      cells = cells,
      verbose = FALSE
    )
    colnames(cutmat) <- start(x = region):end(x = region)
    group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
    scale.factor.use <- scale.factor %||% median(x = group.scale.factors)
    cm.list[[i]] <- cutmat
    sf.list[[i]] <- scale.factor.use
    gsf.list[[i]] <- group.scale.factors
  }
  names(x = cm.list) <- unlist(x = assay)
  p <- CoverageTrack(
    cutmat = cm.list,
    region = region,
    group.scale.factors = gsf.list,
    scale.factor = sf.list,
    window = window,
    ymax = ymax,
    split.assays = split.assays,
    assay.scale = assay.scale,
    obj.groups = obj.groups,
    region.highlight = region.highlight,
    downsample.rate = downsample.rate,
    max.downsample = max.downsample
  )
  # create bigwig tracks
  if (!is.null(x = bigwig)) {
    if (!inherits(x = bigwig, what = "list")) {
      # accept a (named) character vector of file paths
      bigwig <- as.list(x = bigwig)
    }
    if (is.null(x = names(x = bigwig))) {
      names(x = bigwig) <- paste0("bigWig_", seq_along(along.with = bigwig))
    }
    if (length(x = bigwig.type) == 1) {
      bigwig.type <- rep(x = bigwig.type, length(x = bigwig))
    } else if (length(x = bigwig.type) != length(x = bigwig)) {
      stop("Must supply a bigWig track type for each bigWig file")
    }
    unique.types <- unique(x = bigwig.type)
    bw.all <- list()
    for (i in seq_along(unique.types)) {
      bw.use <- which(x = bigwig.type == unique.types[[i]])
      bw.all[[i]] <- BigwigTrack(
        region = region,
        bigwig = bigwig[bw.use],
        type = unique.types[[i]],
        bigwig.scale = bigwig.scale,
        ymax = ymax,
        max.downsample = max.downsample
      )
    }
    bigwig.tracks <- CombineTracks(
      plotlist = bw.all,
      heights = table(unlist(x = bigwig.type))
    )
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
  if (is.logical(x = annotation)) {
    if (annotation) {
      gene.plot <- AnnotationPlot(
        object = object[[assay[[1]]]],
        region = region,
        mode = "gene"
      )
    } else {
      gene.plot <- NULL
    }
  } else {
    gene.plot <- AnnotationPlot(
      object = object[[assay[[1]]]],
      region = region,
      mode = annotation
    )
  }
  if (!is.null(x = links)) {
    linkplot.list <- list()
    for (i in seq_along(along.with = links)) {
      # generate links plot for each key
      linkplot.list[[i]] <- LinkPlot(
        object = object[[assay[[1]]]],
        region = region,
        key = links[[i]]
      )
      if (length(x = links) > 1) {
        linkplot.list[[i]] <- linkplot.list[[i]] + ylab(links[[i]])
      }
    }
    if (length(x = linkplot.list) > 0) {
      link.plot <- CombineTracks(
        plotlist = linkplot.list,
        heights = rep(1, length(x = linkplot.list))
      )
    } else {
      link.plot <- NULL
    }
  } else {
    link.plot <- NULL
  }
  if (peaks && is.granges) {
    peak.plot <- PeakPlot(
      object = object,
      assay = assay[[1]],
      region = region,
      group.by = peaks.group.by
    )
  } else {
    peak.plot <- NULL
  }
  if (!is.null(x = ranges)) {
    range.plot <- PeakPlot(
      object = object,
      assay = assay[[1]],
      region = region,
      peaks = ranges,
      group.by = ranges.group.by,
      color = "brown3"
    ) +
      ylab(ranges.title)
  } else {
    range.plot <- NULL
  }
  if (tile) {
    # reuse cut matrix
    # TODO implement for multi assay
    tile.df <- ComputeTile(
      cutmatrix = cm.list[[1]],
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
      assay = assay[[1]],
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
      cutmat = cm.list,
      region = region,
      # CoverageTrack indexes these per assay, so recycle across the assays
      group.scale.factors = rep(
        x = list(bulk.scale.factor), length(x = cm.list)
      ),
      # sf.list holds the resolved (non-NULL) scale factor for each assay
      scale.factor = sf.list,
      assay.scale = assay.scale,
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

  if (!is.null(x = gwas)) {
    # Convert to list if needed (following bigwig pattern)
    if (!inherits(x = gwas, what = "list")) {
      gwas <- list(gwas)
      names(gwas) <- "GWAS"
    }

    # Handle associated parameters - convert to lists
    if (
      length(x = gwas.ld.file) == 1 ||
        !inherits(x = gwas.ld.file, what = "list")
    ) {
      gwas.ld.file <- rep(list(gwas.ld.file), length(x = gwas))
    }
    if (
      length(x = gwas.ld.lead.snp) == 1 ||
        !inherits(x = gwas.ld.lead.snp, what = "list")
    ) {
      gwas.ld.lead.snp <- rep(list(gwas.ld.lead.snp), length(x = gwas))
    }
    if (
      length(x = gwas.credset.file) == 1 ||
        !inherits(x = gwas.credset.file, what = "list")
    ) {
      gwas.credset.file <- rep(list(gwas.credset.file), length(x = gwas))
    }

    # Create tracks
    gwas.all <- list()
    for (i in seq_along(along.with = gwas)) {
      gwas.all[[i]] <- GWASTrack(
        gwas = gwas[[i]],
        region = region,
        ld.file = gwas.ld.file[[i]],
        ld.lead.snp = gwas.ld.lead.snp[[i]],
        credset.file = gwas.credset.file[[i]],
        credset.threshold = gwas.credset.threshold,
        show.axis = FALSE
      )
      if (length(x = gwas) > 1) {
        gwas.all[[i]] <- gwas.all[[i]] + ylab(label = names(x = gwas)[[i]])
      }
    }

    # Combine tracks (following bigwig pattern)
    gwas.tracks <- CombineTracks(
      plotlist = gwas.all,
      heights = rep(10, length(x = gwas))
    )
  } else {
    gwas.tracks <- NULL
  }

  # variants
  if (!is.null(x = variants)) {
    variant.track <- VariantTrack(variants = variants, region = region)
  } else {
    variant.track <- NULL
  }

  # AlphaGenome Atlas AVI SHAP feature attributions
  if (!is.null(x = avi)) {
    avi.track <- AVITrack(
      avi = avi,
      region = region,
      features = avi.features,
      type = avi.type,
      bins = avi.bins,
      show.axis = FALSE,
      api.key = avi.api.key,
      max.workers = avi.max.workers,
      verbose = avi.verbose
    )
  } else {
    avi.track <- NULL
  }

  nident <- length(x = unique(x = obj.groups))
  if (split.assays) {
    nident <- nident * length(x = assay)
  }
  bulk.height <- (1 / nident) * 10
  bw.height <- 10
  gwas.height <- 3
  variants.height <- 1
  avi.height <- 3
  plot.list <- list(
    gwas.tracks,
    variant.track,
    avi.track,
    p,
    bulk.plot,
    bigwig.tracks,
    tile.plot,
    gene.plot,
    peak.plot,
    range.plot,
    link.plot
  )
  if (is.null(x = heights)) {
    # one entry per slot in plot.list; CombineTracks drops the entries for
    # tracks that were not requested along with the NULL plots
    heights <- c(
      gwas.height,
      variants.height,
      avi.height,
      10,
      bulk.height,
      bw.height,
      10, 3, 1, 1, 3
    )
  } else {
    # a user-supplied vector already refers to the displayed tracks, so drop
    # the unused slots before passing it on
    plot.list <- Filter(f = Negate(f = is.null), x = plot.list)
    if (length(x = heights) != length(x = plot.list)) {
      stop(
        "heights must contain one value for each displayed track (",
        length(x = plot.list), " tracks are shown), or be NULL"
      )
    }
  }
  # the expression plot is aligned with the coverage track, which follows the
  # optional GWAS, variant and AVI tracks
  coverage.index <- 1 + sum(!vapply(
    X = list(gwas.tracks, variant.track, avi.track),
    FUN = is.null,
    FUN.VALUE = logical(length = 1)
  ))
  p <- CombineTracks(
    plotlist = plot.list,
    expression.plot = ex.plot,
    heights = heights,
    widths = widths,
    expression.index = coverage.index
  ) & theme(
    legend.key.size = unit(x = 1 / 2, units = "lines"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
  return(p)
}

# Coverage Track
#
#' @importFrom ggplot2 geom_area geom_hline facet_wrap xlab ylab theme_classic aes ylim theme element_blank element_text geom_segment scale_color_identity scale_fill_manual geom_rect aes
#' @importFrom IRanges IRanges width
#' @importFrom Seqinfo seqnames
#' @importFrom Matrix colSums
#' @importFrom stats median
#' @importFrom dplyr mutate group_by ungroup slice_sample
#' @importFrom RcppRoll roll_sum
#' @importFrom methods is
#' @importFrom scales hue_pal
#' @importFrom S4Vectors mcols
#' @importFrom withr with_seed
#' @importMethodsFrom GenomicRanges start end
CoverageTrack <- function(
  cutmat,
  region,
  group.scale.factors,
  scale.factor,
  assay.scale,
  obj.groups,
  ymax,
  downsample.rate,
  split.assays = FALSE,
  region.highlight = NULL,
  window = 100,
  max.downsample = 3000
) {
  window.size <- width(x = region)
  levels.use <- levels(x = obj.groups)
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  multicov <- length(x = cutmat) > 1

  cov.df <- data.frame()
  for (i in seq_along(along.with = cutmat)) {
    coverages <- ApplyMatrixByGroup(
      mat = cutmat[[i]],
      fun = colSums,
      groups = obj.groups,
      group.scale.factors = group.scale.factors[[i]],
      scale.factor = scale.factor[[i]],
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

    coverages <- coverages[!is.na(x = coverages$coverage), ]
    coverages <- group_by(.data = coverages, group)
    sampling <- min(max.downsample, window.size * downsample.rate)
    coverages <- with_seed(
      seed = 1234,
      code = slice_sample(.data = coverages, n = as.integer(x = sampling))
    )
    coverages$Assay <- names(x = cutmat)[[i]]
    if (multicov) {
      if (assay.scale == "separate") {
        # scale to fraction of max for each separately
        assay.max <- max(coverages$coverage, na.rm = TRUE)
        coverages$coverage <- coverages$coverage / assay.max
      }
    }
    cov.df <- rbind(cov.df, coverages)
  }
  coverages <- cov.df
  coverages$Assay <- factor(x = coverages$Assay, levels = names(x = cutmat))
  coverages$assay_group <- paste(coverages$group, coverages$Assay, sep = "_")

  # restore factor levels
  if (!is.null(x = levels.use)) {
    colors_all <- hue_pal()(length(x = levels.use))
    names(x = colors_all) <- levels.use
    coverages$group <- factor(x = coverages$group, levels = levels.use)
  }
  covmax <- signif(x = max(coverages$coverage, na.rm = TRUE), digits = 2)
  if (is.null(x = ymax)) {
    ymax <- covmax
  } else if (is.character(x = ymax)) {
    if (!startsWith(x = ymax, prefix = "q")) {
      stop("Unknown ymax requested. Must be NULL, a numeric value, or
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
    }
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
    ) / 100
    ymax <- covmax * percentile.use
  }
  ymin <- 0

  # perform clipping
  coverages$coverage[coverages$coverage > ymax] <- ymax

  if (multicov) {
    p <- ggplot(
      data = coverages,
      mapping = aes(x = .data[["position"]], y = .data[["coverage"]], fill = .data[["Assay"]])
    )
  } else {
    p <- ggplot(
      data = coverages,
      mapping = aes(x = .data[["position"]], y = .data[["coverage"]], fill = .data[["group"]])
    )
  }
  p <- p +
    geom_area(
      stat = "identity",
      alpha = ifelse(test = !split.assays & multicov, yes = 0.5, no = 1)
    ) +
    geom_hline(yintercept = 0, linewidth = 0.1)
  if (split.assays) {
    p <- p +
      facet_wrap(facets = ~assay_group, strip.position = "left", ncol = 1)
  } else {
    p <- p + facet_wrap(facets = ~group, strip.position = "left", ncol = 1)
  }
  p <- p +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = paste0(
      "Normalized signal \n(range ",
      as.character(x = ymin), " - ",
      as.character(x = ymax), ")"
    )) +
    ylim(c(ymin, ymax)) +
    theme_browser(legend = multicov) +
    theme(panel.spacing.y = unit(x = 0, units = "line"))
  if (!is.null(x = levels.use) && !multicov) {
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
      n.panels <- if (split.assays) {
        length(x = unique(x = coverages$assay_group))
      } else {
        length(x = unique(x = coverages$group))
      }
      p <- p +
        geom_rect(
          data = df,
          inherit.aes = FALSE,
          aes(
            xmin = .data[["start"]],
            xmax = .data[["end"]],
            ymin = 0,
            ymax = ymax
          ),
          fill = rep(x = df$color, n.panels),
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
#' given regions of the genome. Tracks are normalized using a per-group scaling
#' factor computed as the number of cells in the group multiplied by the mean
#' sequencing depth for that group of cells. This accounts for differences in
#' number of cells and potential differences in sequencing depth between groups.
#'
#' Additional information can be layered on the coverage plot by setting several
#' different options in the CoveragePlot function. This includes showing:
#' \itemize{
#' \item{gene annotations}
#' \item{peak positions}
#' \item{additional genomic ranges}
#' \item{additional data stored in a bigWig file, which may be hosted remotely}
#' \item{gene or protein expression data alongside coverage tracks}
#' \item{peak-gene links}
#' \item{the position of individual sequenced fragments as a heatmap}
#' \item{data for multiple chromatin assays simultaneously}
#' \item{a pseudobulk for all cells combined}
#' }
#'
#' @param object A Seurat object
#' @param region A set of genomic coordinates to show. Can be a GRanges object,
#' a string encoding a genomic position, a gene name, or a vector of strings
#' describing the genomic coordinates or gene names to plot. If a gene name is
#' supplied, annotations must be present in the assay. If a GRanges object
#' containing more than one range is supplied, a separate plot is produced for
#' each range.
#' @param features A vector of features present in another assay to plot
#' alongside accessibility tracks (for example, gene names).
#' @param assay Name of the assay to plot. If a list of assays is provided,
#' data from each assay will be shown overlaid on each track. The first assay in
#' the list will define the assay used for gene annotations, links, and peaks
#' (if shown). The order of assays given defines the plotting order.
#' @param split.assays When plotting data from multiple assays, display each
#' assay as a separate track. If FALSE, data from different assays are overlaid
#' on a single track with transparancy applied.
#' @param assay.scale Scaling to apply to data from different assays. Can be:
#' \itemize{
#' \item{common: plot all assays on a common scale (default)}
#' \item{separate: plot each assay on a separate scale ranging from zero to the
#' maximum value for that assay within the plotted region}
#' }
#' @param show.bulk Include coverage track for all cells combined (pseudo-bulk).
#' Note that this will plot the combined accessibility for all cells included in
#' the plot (rather than all cells in the object).
#' @param expression.assay Name of the assay containing expression data to plot
#' alongside accessibility tracks. Only needed if supplying `features`
#' argument.
#' @param expression.slot Name of slot to pull expression data from. Only needed
#' if supplying the `features` argument.
#' @param annotation Display gene annotations. Set to TRUE or FALSE to control
#' whether genes models are displayed, or choose "transcript" to display all
#' transcript isoforms, or "gene" to display gene models only (same as setting
#' TRUE).
#' @param peaks Display peaks
#' @param peaks.group.by Grouping variable to color peaks by. Must be a variable
#' present in the feature metadata. If NULL, do not color peaks by any variable.
#' @param ranges Additional genomic ranges to plot
#' @param ranges.group.by Grouping variable to color ranges by. Must be a
#' variable present in the metadata stored in the `ranges` genomic ranges.
#' If NULL, do not color by any variable.
#' @param ranges.title Y-axis title for ranges track. Only relevant if
#' `ranges` parameter is set.
#' @param region.highlight Region to highlight on the plot. Should be a GRanges
#' object containing the coordinates to highlight. By default, regions will be
#' highlighted in grey. To change the color of the highlighting, include a
#' metadata column in the GRanges object named "color" containing the color to
#' use for each region.
#' @param links Character vector containing the keys of link information present
#' in the assay to display. Default is "linkpeaks" which is the default key for
#' peak-gene links stored using the [LinkPeaks()] function. If `NULL`, links
#' will not be displayed.
#' @param tile Display per-cell fragment information in sliding windows. If
#' plotting multi-assay data, only the first assay is shown in the tile plot.
#' @param tile.size Size of the sliding window for per-cell fragment tile plot
#' @param tile.cells Number of cells to display fragment information for in tile
#' plot.
#' @param bigwig A named list or named character vector of bigWig file paths to
#' plot data from. Files can be remotely hosted. The name of each element
#' determines the y-axis label given to the track.
#' @param bigwig.type Type of track to use for bigWig files ("line", "heatmap",
#' or "coverage"). Should either be a single value, or a list of values giving
#' the type for each individual track in the provided list of bigwig files.
#' @param bigwig.scale Same as `assay.scale` parameter, except for bigWig
#' files when plotted with `bigwig.type="coverage"`
#' @param cells Which cells to plot. Default all cells
#' @param idents Which identities to include in the plot. Default is all
#' identities.
#' @param window Smoothing window size
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' @param ymax Maximum value for Y axis. Can be one of:
#'  - `NULL`: set to the highest value among all the tracks (default)
#'  - qXX: clip the maximum value to the XX quantile (for example, q95 will
#' set the maximum value to 95% of the maximum value in the data). This can
#' help remove the effect of extreme values that may otherwise distort the
#' scale.
#'  - numeric: manually define a Y-axis limit
#' @param scale.factor Scaling factor for track height. If NULL (default),
#' use the median group scaling factor determined by total number of fragments
#' sequences in each group.
#' @param group.by Name of one or more metadata columns to group (color) the
#' cells by. Default is the current cell identities
#' @param split.by A metadata variable to split the tracks by. For example,
#' grouping by "celltype" and splitting by "batch" will create separate tracks
#' for each combination of celltype and batch.
#' @param heights Relative heights for each displayed track, given in plotting
#' order (top to bottom): GWAS, variants, AVI SHAP, accessibility, pseudobulk,
#' bigWig, tile, gene annotations, peaks, ranges, links. Only the tracks actually shown
#' are counted, so the vector must have one value per displayed track. If
#' `NULL` (default), heights are chosen automatically.
#' @param max.downsample Maximum amount of downsampling to apply. Corresponds to
#' the minimum number of positions to be kept when downsampling. The
#' downsampling rate is adaptive to the window size, but this parameter will set
#' the minimum possible number of positions to include so that plots do not
#' become too sparse when the window size is small.
#' @param downsample.rate Fraction of positions to retain when downsampling.
#' Retaining more positions can give a higher-resolution plot but can make the
#' number of points large, resulting in larger file sizes when saving the plot
#' and a longer period of time needed to draw the plot.
#' @param gwas GWAS summary statistics to display on the plot. Can be the path
#' to a GWAS-SSF file on-disk or a dataframe in the GWAS-SSF format.
#' @param gwas.ld.file Path to LD data file for coloring GWAS points by r².
#' Optional.
#' @param gwas.ld.lead.snp Lead SNP for LD calculations.
#' Required if gwas.ld.file provided.
#' @param gwas.credset.file Path to fine-mapping credible sets file. Optional.
#' @param gwas.credset.threshold Posterior probability threshold for credible
#' sets (default: 0.01)
#' @param variants Dataframe containing variants to display
#' (see [VariantTrack()])
#' @param avi AlphaGenome Atlas AVI SHAP feature attributions to display on
#' the plot. Either a dataframe with one row per variant (see [AVITrack()]
#' and [LoadAVIFromAtlas()]) or `TRUE` to query the AlphaGenome Atlas API for
#' the plotted region. `NULL` or `FALSE` disables the track. The AlphaGenome
#' Atlas covers hg38 only: if the genome of the gene annotations stored in
#' the assay (`genome(Annotation(object))`) is set to a different build, the
#' track is skipped with a warning.
#' @param avi.features AVI SHAP features to plot, given as feature column
#' names and/or modality group labels (see [AVITrack()]), or "Total"
#' for the total AVI score. Groups without any feature in the data are skipped.
#' Set to `NULL` to plot all features.
#' @param avi.type Plot type for the AVI SHAP track. One of "stacked" or
#' "heatmap" (see [AVITrack()]).
#' @param avi.bins Number of bins (bars) across the region for the AVI SHAP
#' track.
#' @param avi.api.key AlphaGenome API key, used when `avi = TRUE`. If
#' NULL, the `ALPHA_GENOME_API_KEY` environment variable is used.
#' @param avi.max.workers Number of concurrent Atlas API requests.
#' @param avi.verbose Show Atlas query progress and fallback messages.
#' @param ... Additional arguments passed to [patchwork::wrap_plots()]
#'
#' @importFrom patchwork wrap_plots
#' @importFrom SeuratObject DefaultAssay
#' @importFrom Seqinfo genome
#' @export
#' @concept visualization
#' @return Returns a [patchwork::patchwork()] object
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#'
#' # Basic coverage plot
#' CoveragePlot(object = atac_small, region = c("chr1:713500-714500"))
#'
#' # Show additional ranges
#' ranges.show <- GenomicRanges::GRanges("chr1:713750-714000")
#' CoveragePlot(
#'   object = atac_small,
#'   region = c("chr1:713500-714500"),
#'   ranges = ranges.show
#' )
#'
#' # Highlight region
#' CoveragePlot(
#'   object = atac_small,
#'   region = c("chr1:713500-714500"),
#'   region.highlight = ranges.show
#' )
#'
#' # Change highlight color
#' ranges.show$color <- "orange"
#' CoveragePlot(
#'   object = atac_small,
#'   region = c("chr1:713500-714500"),
#'   region.highlight = ranges.show
#' )
#'
#' # Show expression data
#' CoveragePlot(
#'   object = atac_small,
#'   region = c("chr1:713500-714500"),
#'   features = "GYG2"
#' )
#' }
CoveragePlot <- function(
  object,
  region,
  features = NULL,
  assay = NULL,
  split.assays = FALSE,
  assay.scale = "common",
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
  links = "linkpeaks",
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  bigwig.scale = "common",
  heights = NULL,
  group.by = NULL,
  split.by = NULL,
  window = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  scale.factor = NULL,
  ymax = NULL,
  cells = NULL,
  idents = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1,
  gwas = NULL,
  gwas.ld.file = NULL,
  gwas.ld.lead.snp = NULL,
  gwas.credset.file = NULL,
  gwas.credset.threshold = 0.01,
  variants = NULL,
  avi = NULL,
  avi.features = c("Expression", "Accessibility", "Total"),
  avi.type = "stacked",
  avi.bins = 200,
  avi.api.key = NULL,
  avi.max.workers = 4,
  avi.verbose = FALSE,
  ...
) {
  if (inherits(x = region, what = "GRanges") && length(x = region) > 1) {
    region <- as.list(x = split(x = region, f = seq_along(region)))
  } else if (length(x = region) == 1) {
    region <- list(region)
  }
  if (isFALSE(x = avi)) {
    avi <- NULL
  }
  if (!is.null(x = avi)) {
    # the AlphaGenome Atlas covers hg38 only
    avi.assay <- if (inherits(x = assay, what = "list")) assay[[1]] else assay
    avi.assay <- avi.assay %||% DefaultAssay(object = object)
    build <- unique(x = genome(x = Annotation(object = object[[avi.assay]])))
    build <- build[!is.na(x = build)]
    is.hg38 <- grepl(pattern = "^(hg38|grch38)", x = tolower(x = build))
    if (length(x = build) > 0 && !any(is.hg38)) {
      warning(
        "The AlphaGenome Atlas provides scores for hg38 only, but the genome ",
        "of the gene annotations in the '", avi.assay, "' assay is ",
        paste(build, collapse = ", "), ". The AVI SHAP track will not be shown.",
        immediate. = TRUE, call. = FALSE
      )
      avi <- NULL
    }
  }
  plot.list <- lapply(
    X = seq_along(region),
    FUN = function(x) {
      SingleCoveragePlot(
        object = object,
        region = region[[x]],
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
        split.assays = split.assays,
        assay.scale = assay.scale,
        links = links,
        tile = tile,
        tile.size = tile.size,
        tile.cells = tile.cells,
        bigwig = bigwig,
        bigwig.type = bigwig.type,
        bigwig.scale = bigwig.scale,
        group.by = group.by,
        split.by = split.by,
        window = window,
        ymax = ymax,
        scale.factor = scale.factor,
        extend.upstream = extend.upstream,
        extend.downstream = extend.downstream,
        cells = cells,
        idents = idents,
        heights = heights,
        max.downsample = max.downsample,
        downsample.rate = downsample.rate,
        gwas = gwas,
        gwas.ld.file = gwas.ld.file,
        gwas.ld.lead.snp = gwas.ld.lead.snp,
        gwas.credset.file = gwas.credset.file,
        gwas.credset.threshold = gwas.credset.threshold,
        variants = variants,
        avi = avi,
        avi.features = avi.features,
        avi.type = avi.type,
        avi.bins = avi.bins,
        avi.api.key = avi.api.key,
        avi.max.workers = avi.max.workers,
        avi.verbose = avi.verbose
      )
    }
  )
  return(wrap_plots(plot.list, ...))
}

#' Plot DNA sequence motif
#'
#' Plot position weight matrix or position frequency matrix for different DNA
#' sequence motifs.
#'
#' @param object A Seurat object
#' @param motifs A list of motif IDs or motif names to plot
#' @param assay Name of the assay to use
#' @param use.names Use motif names stored in the motif object
#' @param ... Additional parameters passed to [ggseqlogo::ggseqlogo()]
#'
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept visualization
#' @concept motifs
#' @return Returns a [ggplot2::ggplot()] object
#' @examples
#' \donttest{
#' motif.obj <- Motifs(atac_small)
#' MotifPlot(atac_small, motifs = head(colnames(motif.obj)))
#' }
MotifPlot <- function(
  object,
  motifs,
  assay = NULL,
  use.names = TRUE,
  ...
) {
  if (!inherits(x = motifs, what = "character")) {
    stop("Please provide motif names, not a ", class(x = motifs), " vector")
  }
  if (!requireNamespace(package = "ggseqlogo", quietly = TRUE)) {
    stop("Please install ggseqlogo: install.packages('ggseqlogo')")
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "GRangesAssay")) {
    stop("The requested assay is not a GRangesAssay.")
  }
  data.use <- GetMotifData(object = object, assay = assay, slot = "pwm")
  if (length(x = data.use) == 0) {
    stop("Position weight matrix list for the requested assay is empty")
  }
  requested.motifs <- motifs
  missing.motifs <- !(motifs %in% names(x = data.use))
  for (i in seq_along(along.with = motifs)) {
    if (missing.motifs[i]) {
      # try looking up ID
      motifs[i] <- ConvertMotifID(object = object, name = motifs[i])
    }
  }
  # motifs not found by name or ID are NA after conversion
  not.found <- is.na(x = motifs) | !(motifs %in% names(x = data.use))
  if (any(not.found)) {
    if (all(not.found)) {
      stop(
        "None of the requested motifs were found in the ", assay,
        " assay: ", paste(requested.motifs, collapse = ", "),
        call. = FALSE
      )
    }
    warning(
      "The following requested motifs were not found in the ", assay,
      " assay and will be skipped: ",
      paste(requested.motifs[not.found], collapse = ", "),
      call. = FALSE
    )
    motifs <- motifs[!not.found]
  }
  data.use <- data.use[motifs]
  if (use.names) {
    names(x = data.use) <- GetMotifData(
      object = object, assay = assay, slot = "motif.names"
    )[motifs]
  }
  p <- ggseqlogo::ggseqlogo(data = data.use, ...)
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
#' @param max.length Maximum fragment length to display. Fragments longer than
#' this are not shown. If `NULL`, show all fragment lengths present.
#' @param ... Arguments passed to other functions
#'
#' @importFrom ggplot2 ggplot geom_histogram theme_classic aes facet_wrap scale_y_log10 theme element_blank xlim
#' @importFrom SeuratObject DefaultAssay
#'
#' @export
#' @concept visualization
#' @concept qc
#' @return Returns a [ggplot2::ggplot()] object
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' FragmentHistogram(object = atac_small, region = "chr1:10245-780007")
#' }
FragmentHistogram <- function(
  object,
  assay = NULL,
  region = "chr1:1-2000000",
  group.by = NULL,
  cells = NULL,
  log.scale = FALSE,
  max.length = 800,
  ...
) {
  cells <- cells %||% colnames(x = object)
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("The requested assay is not a ChromatinAssay5.")
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
  if (!is.null(x = max.length)) {
    n.dropped <- sum(reads$length > max.length)
    if (n.dropped > 0) {
      message(
        n.dropped, " fragments longer than max.length (", max.length,
        " bp) are not shown"
      )
      reads <- reads[reads$length <= max.length, ]
    }
  }
  if (length(x = unique(x = reads$group)) == 1) {
    p <- ggplot(data = reads, mapping = aes(x = .data[["length"]])) +
      geom_histogram(bins = 200)
  } else {
    p <- ggplot(
      data = reads,
      mapping = aes(x = .data[["length"]], fill = .data[["group"]])
    ) +
      geom_histogram(bins = 200) +
      facet_wrap(~group, scales = "free_y")
  }
  if (!is.null(x = max.length)) {
    p <- p + xlim(c(0, max.length))
  }
  p <- p +
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
#' @param expression.index Position in `plotlist` (after removing `NULL`
#' entries) of the plot that the expression plot is aligned with. Default: the
#' first plot.
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
#' p1 <- PeakPlot(atac_small, region = "chr1:29554-39554")
#' p2 <- AnnotationPlot(atac_small, region = "chr1:29554-39554")
#' CombineTracks(plotlist = list(p1, p2), heights = c(1, 1))
#' }
CombineTracks <- function(
  plotlist,
  expression.plot = NULL,
  heights = NULL,
  widths = NULL,
  expression.index = 1
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
    # two-column grid: the expression plot sits beside the track at
    # expression.index; tracks above and below are stacked in nested
    # patchworks with an empty cell (or the legend) beside them
    n <- length(x = plotlist)
    i <- max(1, min(expression.index, n))
    cells <- list()
    row.heights <- c()
    if (i > 1) {
      above <- seq_len(length.out = i - 1)
      cells <- c(cells, list(
        wrap_plots(plotlist[above], ncol = 1, heights = heights[above]),
        plot_spacer()
      ))
      row.heights <- c(row.heights, sum(heights[above]))
    }
    cells <- c(cells, list(plotlist[[i]], expression.plot))
    row.heights <- c(row.heights, heights[[i]])
    if (i < n) {
      below <- (i + 1):n
      cells <- c(cells, list(
        wrap_plots(plotlist[below], ncol = 1, heights = heights[below]),
        guide_area()
      ))
      row.heights <- c(row.heights, sum(heights[below]))
    } else {
      cells <- c(cells, list(plot_spacer(), guide_area()))
      row.heights <- c(row.heights, 0.01)
    }
    p <- wrap_plots(
      cells, ncol = 2, heights = row.heights, widths = widths,
      guides = "collect"
    )
  } else {
    p <- wrap_plots(plotlist, ncol = 1, heights = heights)
  }
  return(p)
}

#' Plot peaks in a genomic region
#'
#' Display the genomic ranges in a [GRangesAssay-class] object that fall
#' in a given genomic region
#'
#' @param object A [SeuratObject::Seurat()] object
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param region A genomic region to plot
#' @param peaks A GRanges object containing peak coordinates. If NULL, use
#' coordinates stored in the Seurat object.
#' @param group.by Name of variable in feature metadata (if using ranges in the
#' Seurat object) or genomic ranges metadata (if using supplied ranges) to color
#' ranges by. If NULL, do not color by any metadata variable.
#' @param color Color to use. If `group.by` is not NULL, this can be a
#' custom color scale (see examples).
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#'
#' @return Returns a [ggplot2::ggplot()] object
#' @export
#' @concept visualization
#' @importFrom SeuratObject DefaultAssay
#' @importFrom S4Vectors mcols<-
#' @importFrom GenomicRanges start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom Seqinfo seqnames
#' @importFrom ggplot2 ggplot aes geom_segment theme_classic theme xlab ylab scale_color_manual element_blank
#' @examples
#' \donttest{
#' # plot peaks in assay
#' PeakPlot(atac_small, region = "chr1:710000-715000")
#'
#' # manually set color
#' PeakPlot(atac_small, region = "chr1:710000-715000", color = "red")
#'
#' # color by a variable in the feature metadata
#' PeakPlot(atac_small, region = "chr1:710000-715000", group.by = "count")
#' }
PeakPlot <- function(
  object,
  region,
  assay = NULL,
  peaks = NULL,
  group.by = NULL,
  color = "dimgrey",
  extend.upstream = 0,
  extend.downstream = 0
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("The requested assay is not a ChromatinAssay5.")
  }

  if (!inherits(x = region, what = "GRanges")) {
    region <- GRanges(region)
  }
  if (is.null(x = peaks)) {
    if (!inherits(x = object[[assay]], what = "GRangesAssay")) {
      stop("The requested assay is not a GRangesAssay.")
    } else {
      peaks <- granges(x = object[[assay]])
      md <- object[[assay]][[]]
      mcols(x = peaks) <- md
    }
  }
  region <- FindRegion(
    object = object,
    region = region,
    assay = assay[[1]],
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
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
    if (!is.null(x = group.by)) {
      peak.plot <- ggplot(
        data = peak.df,
        aes(color = .data[[group.by]])
      )
    } else {
      peak.plot <- ggplot(data = peak.df)
    }
    if (!is.null(x = group.by)) {
      peak.plot <- peak.plot +
        geom_segment(aes(x = .data[["start"]], y = 0, xend = .data[["end"]], yend = 0),
          linewidth = 2,
          data = peak.df
        )
    } else {
      peak.plot <- peak.plot +
        geom_segment(aes(x = .data[["start"]], y = 0, xend = .data[["end"]], yend = 0),
          linewidth = 2,
          color = color,
          data = peak.df
        )
    }
  } else {
    # no peaks present in region, make empty panel
    peak.plot <- ggplot(data = peak.df)
  }
  peak.plot <- peak.plot + theme_classic() +
    ylab(label = "Peaks") +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start.pos, end.pos))
  return(peak.plot)
}

#' Plot linked genomic elements
#'
#' Display links between pairs of genomic elements within a given region of the
#' genome.
#'
#' @param object A [SeuratObject::Seurat()] object
#' @param region A genomic region to plot
#' @param key Key to use when extracting link information from the Seurat object
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param min.cutoff Minimum absolute score for link to be plotted.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#' @param scale.linewidth Scale thickness of the line according to link score.
#'
#'
#' @return Returns a [ggplot2::ggplot()] object
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom Seqinfo seqnames
#' @importFrom ggplot2 ggplot geom_hline theme_classic xlim ylab theme element_blank scale_color_gradient2 aes
#' @concept visualization
#' @concept links
LinkPlot <- function(
  object,
  region,
  key,
  assay = NULL,
  min.cutoff = 0,
  extend.upstream = 0,
  extend.downstream = 0,
  scale.linewidth = FALSE
) {
  region <- FindRegion(
    object = object,
    region = region,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  chromosome <- seqnames(x = region)

  # extract link information
  links <- Links(object = object)
  links <- links[[key]]

  # if links not set, return NULL
  if (length(x = links) == 0) {
    return(NULL)
  }

  # subset to those in region
  links.keep <- subsetByOverlaps(x = links, ranges = region)

  # filter out links below threshold
  link.df <- as.data.frame(x = links.keep)
  link.df <- link.df[abs(x = link.df$score) > min.cutoff, ]

  # convert to single start and end position
  # take start of gene and midpoint of peak
  link.df$peak_midpoint <- (link.df$start1 + link.df$end1) / 2
  link.df$genestart <- ifelse(
    test = link.df$strand2 == "-",
    yes = link.df$end2,
    no = link.df$start2
  )
  link.df$start <- ifelse(
    test = link.df$peak_midpoint < link.df$genestart,
    yes = link.df$peak_midpoint,
    no = link.df$genestart
  )
  link.df$end <- ifelse(
    test = link.df$peak_midpoint < link.df$genestart,
    yes = link.df$genestart,
    no = link.df$peak_midpoint
  )

  # remove links outside region
  link.df <- link.df[
    link.df$start >= start(x = region) & link.df$end <= end(x = region),
  ]

  # plot
  if (nrow(x = link.df) > 0) {
    if (!requireNamespace(package = "ggforce", quietly = TRUE)) {
      warning(
        "Please install ggforce to enable LinkPlot plotting: ",
        "install.packages('ggforce')"
      )
      p <- ggplot(data = link.df)
    } else {
      # convert to format for geom_bezier
      link.df$group <- seq_len(length.out = nrow(x = link.df))
      df <- data.frame(
        x = c(
          link.df$start,
          (link.df$start + link.df$end) / 2,
          link.df$end
        ),
        y = c(
          rep(x = 0, nrow(x = link.df)),
          rep(x = -1, nrow(x = link.df)),
          rep(x = 0, nrow(x = link.df))
        ),
        group = rep(x = link.df$group, 3),
        score = rep(link.df$score, 3)
      )
      min.color <- min(0, min(df$score))
      if (scale.linewidth) {
        p <- ggplot(data = df) +
          ggforce::geom_bezier(
            mapping = aes(
              x = .data[["x"]],
              y = .data[["y"]],
              group = .data[["group"]],
              color = .data[["score"]],
              linewidth = .data[["score"]]
            )
          )
      } else {
        p <- ggplot(data = df) +
          ggforce::geom_bezier(
            mapping = aes(
              x = .data[["x"]],
              y = .data[["y"]],
              group = .data[["group"]],
              color = .data[["score"]]
            )
          )
      }
      p <- p +
        geom_hline(yintercept = 0, color = "grey") +
        scale_color_gradient2(
          low = "red", mid = "grey", high = "blue",
          limits = c(min.color, max(df$score)),
          n.breaks = 3
        )
    }
  } else {
    p <- ggplot(data = link.df)
  }
  p <- p +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    ylab("Links") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
  return(p)
}

#' Plot gene annotations
#'
#' Display gene annotations in a given region of the genome.
#'
#' @param object A [SeuratObject::Seurat()] object
#' @param region A genomic region to plot
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param mode Display mode. Choose either "gene" or "transcript" to determine
#' whether genes or transcripts are plotted.
#' @param extend.upstream Number of bases to extend the region upstream.
#' @param extend.downstream Number of bases to extend the region downstream.
#'
#' @return Returns a [ggplot2::ggplot()] object
#' @export
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom Seqinfo seqnames
#' @importFrom ggplot2 theme_classic ylim xlim ylab xlab geom_segment geom_text aes scale_color_manual
#' @importFrom grid arrow
#' @importFrom S4Vectors split
#' @importFrom fastmatch fmatch
#' @concept visualization
#' @examples
#' \donttest{
#' AnnotationPlot(object = atac_small, region = c("chr1:29554-39554"))
#' }
AnnotationPlot <- function(
  object,
  region,
  assay = NULL,
  mode = "gene",
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (mode == "gene") {
    collapse_transcript <- TRUE
    label <- "gene_name"
  } else if (mode == "transcript") {
    collapse_transcript <- FALSE
    label <- "tx_id"
  } else {
    stop("Unknown mode requested, choose either 'gene' or 'transcript'")
  }
  annotation <- Annotation(object = object)
  if (is.null(x = annotation)) {
    return(NULL)
  }
  region <- FindRegion(
    object = object,
    region = region,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)

  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  if (mode == "gene") {
    genes.keep <- unique(x = annotation.subset$gene_name)
    annotation.subset <- annotation[
      fmatch(x = annotation$gene_name, table = genes.keep, nomatch = 0L) > 0L
    ]
  } else {
    tx.keep <- unique(x = annotation.subset$tx_id)
    annotation.subset <- annotation[
      fmatch(x = annotation$tx_id, table = tx.keep, nomatch = 0L) > 0L
    ]
  }

  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      collapse_transcript = collapse_transcript
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes(
          x = .data[["start"]],
          y = annotation_df_list$exons$dodge,
          xend = .data[["end"]],
          yend = annotation_df_list$exons$dodge,
          color = .data[["strand"]]
        ),
        show.legend = FALSE,
        linewidth = 3
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes(
          x = .data[["start"]],
          y = .data[["dodge"]],
          xend = .data[["end"]],
          yend = .data[["dodge"]],
          color = .data[["strand"]]
        ),
        show.legend = FALSE,
        linewidth = 1 / 2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes(
          x = .data[["start"]],
          y = annotation_df_list$plus$dodge,
          xend = .data[["end"]],
          yend = annotation_df_list$plus$dodge,
          color = .data[["strand"]]
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        linewidth = 1 / 2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes(
          x = .data[["start"]],
          y = annotation_df_list$minus$dodge,
          xend = .data[["end"]],
          yend = annotation_df_list$minus$dodge,
          color = .data[["strand"]]
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.04, units = "inches")
        ),
        show.legend = FALSE,
        linewidth = 1 / 2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + 0.2
    p <- p + geom_text(
      data = annotation_df_list$labels,
      mapping = aes(x = .data[["position"]], y = .data[["dodge"]], label = .data[[label]]),
      size = 2.5
    )
    y_limit <- c(0.9, n_stack + 0.4)
  }
  p <- p +
    theme_classic() +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    xlim(start.pos, end.pos) +
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    # name the values so that a strand keeps the same colour whether or not
    # the other strands are present in the region
    scale_color_manual(
      values = c("+" = "darkblue", "-" = "darkgreen", "*" = "grey40")
    )
  return(p)
}

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
#' @param layer Name of layer to pull expression data from
#' @param group.by A grouping variable to group cells by. If NULL, use the
#' current cell identities
#' @param idents A list of identities to include in the plot. If NULL, include
#' all identities
#' @param slot Which slot to pull expression data from
#'
#' @importFrom SeuratObject LayerData DefaultAssay as.sparse
#' @importFrom ggplot2 ggplot geom_violin facet_wrap aes theme_classic element_blank scale_y_discrete scale_x_continuous scale_fill_manual theme
#' @importFrom scales hue_pal
#' @importFrom Seqinfo seqnames
#' @importFrom IRanges start end
#' @importFrom patchwork wrap_plots
#' @importFrom fastmatch fmatch
#' @importFrom lifecycle deprecated deprecate_soft
#'
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' ExpressionPlot(atac_small, features = "ASMTL", assay = "RNA")
#' }
ExpressionPlot <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  layer = "data",
  slot = deprecated()
) {
  if (is_present(arg = slot)) {
    deprecate_soft(
      when = "2.0.0",
      what = "ExpressionPlot(slot = )",
      with = "ExpressionPlot(layer = )"
    )
    layer <- slot %||% layer
  }
  assay <- assay %||% DefaultAssay(object = object)
  data.plot <- LayerData(object = object, assay = assay, layer = layer)
  common.features <- intersect(features, rownames(x = data.plot))
  if (length(x = common.features) == 0) {
    stop("None of the requested features were found in the assay")
  } else if (length(x = common.features) != length(x = features)) {
    warning(
      "Some features not found: ",
      setdiff(x = features, y = rownames(x = data.plot))
    )
    features <- common.features
  }
  data.plot <- data.plot[features, , drop = FALSE]
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = NULL
  )
  obj.groups <- obj.groups[colnames(object[[assay]])]
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
    data.plot <- data.plot[, cells.keep, drop = FALSE]
    obj.groups <- obj.groups[cells.keep]
  }
  # construct data frame
  if (inherits(x = data.plot, what = "IterableMatrix")) {
    data.plot <- as.sparse(x = data.plot)
  }
  if (length(x = features) == 1) {
    df <- data.frame(
      gene = features,
      expression = as.vector(x = data.plot),
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
  missing.levels <- setdiff(x = levels(x = df$group), y = unique(x = df$group))
  if (!is.null(x = idents)) {
    missing.levels <- intersect(x = missing.levels, y = idents)
  }
  if (length(x = missing.levels) > 0) {
    # fill missing idents with zero
    for (i in features) {
      df.1 <- data.frame(
        gene = i,
        expression = 0,
        group = missing.levels
      )
      df <- rbind(df, df.1)
    }
  }
  p.list <- list()
  lower.limit <- ifelse(test = layer == "scale.data", yes = NA, no = 0)
  for (i in seq_along(along.with = features)) {
    df.use <- df[df$gene == features[[i]], ]
    p <- ggplot(data = df.use, aes(
      x = .data[["expression"]], y = .data[["gene"]], fill = .data[["group"]]
    )) +
      geom_violin(linewidth = 1 / 4) +
      facet_wrap(~group, ncol = 1, strip.position = "right") +
      theme_classic() +
      scale_y_discrete(position = "top") +
      scale_x_continuous(position = "bottom", limits = c(lower.limit, NA)) +
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

#' Plot strand concordance vs. VMR
#'
#' Plot the Pearson correlation between allele frequencies on each strand
#' versus the log10 mean-variance ratio for the allele.
#'
#' @param variants A dataframe containing variant information. This should be
#' computed using [IdentifyVariants()]
#' @param min.cells Minimum number of high-confidence cells detected with the
#' variant for the variant to be displayed.
#' @param concordance.threshold Strand concordance threshold
#' @param vmr.threshold Mean-variance ratio threshold
#'
#' @concept mito
#' @concept visualization
#' @export
#' @importFrom ggplot2 ggplot aes geom_point labs scale_y_log10 geom_vline geom_hline theme_classic scale_color_manual theme
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
    mapping = aes(x = .data[["strand_correlation"]], y = .data[["vmr"]], color = .data[["pos"]])
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

#' Plot variant positions
#'
#' Plot variant positions within a genomic region.
#'
#' @param variants Data frame with columns: position (numeric), rsid
#' (character), color (character). Each row defines one SNP marker to display.
#' @param region Genomic region ([GenomicRanges::GRanges] or a string like
#' "chr10:112900000-113100000" that can be converted to `GRanges`)
#'
#' @return Returns a ggplot2 object
#' @export
#' @concept visualization
#' @importFrom ggplot2 geom_segment geom_text ylim margin labs aes
#' @examples
#' # Define SNPs to mark
#' variants <- data.frame(
#'   position = c(112951996, 112952395),
#'   rsid = c("rs10885396", "rs7094871"),
#'   color = c("steelblue", "darkred")
#' )
#'
#' # Create stacked plot
#' VariantTrack(
#'   variants = variants,
#'   region = "chr10:112990000-113010000"
#' )
VariantTrack <- function(
  variants,
  region
) {
  if (!inherits(x = region, what = "GRanges")) {
    region <- GRanges(region)
  }

  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- start(x = region)
  end.pos <- end(x = region)

  snp_plot <- ggplot(data = variants) +
    geom_segment(
      aes(
        x = .data[["position"]],
        xend = .data[["position"]],
        y = 0,
        yend = 1,
        color = .data[["color"]]
      ),
      linewidth = 1,
    ) +
    geom_text(
      aes(x = .data[["position"]], y = 1.2, label = .data[["rsid"]]),
      size = 3.5, fontface = "italic"
    ) +
    scale_color_identity() +
    theme_browser() +
    xlim(start.pos, end.pos) +
    ylim(0, 1.5) +
    labs(
      x = paste0(chromosome, " position (bp)"),
      y = "Variants"
    ) +
    theme(plot.margin = margin(t = 5, r = 5, b = 0, l = 5))

  return(snp_plot)
}

#' Plot AlphaGenome Atlas AVI SHAP feature attributions
#'
#' Plot AlphaGenome Variant Impact (AVI) SHAP feature attributions from the
#' AlphaGenome Atlas across a genomic region. For each variant, the AlphaGenome
#' Atlas decomposes the raw AVI score into 18 additive SHAP contributions, one
#' per input feature (chromatin accessibility, TF binding, histone marks,
#' expression, splicing, protein impact, conservation, etc.). This function
#' summarizes those attributions along the genome and draws them as a track
#' that can be combined with other Signac genome browser tracks (see
#' [CoveragePlot()]).
#'
#' The 18 AVI input features (named as in the AlphaGenome Atlas publication
#' and API) are grouped by modality for display:
#'
#' * **Accessibility**: `MAX_ABS_ATAC`, `MAX_ABS_DNASE`
#' * **TF binding**: `MAX_ABS_CHIP_TF`
#' * **Histone**: `MAX_ABS_CHIP_HISTONE`
#' * **Expression**: `MAX_ABS_CAGE`, `MAX_ABS_PROCAP`, `MAX_ABS_RNA_SEQ`
#' * **Polyadenylation**: `MAX_ABS_POLYADENYLATION`
#' * **Splicing**: `MERGED_SPLICING`
#' * **Contact maps**: `MAX_ABS_CONTACT_MAPS`
#' * **Protein**: `ALPHAMISSENSE`, `PROTEIN_TERMINATION`, `START_LOST`,
#'   `STOP_LOST`
#' * **Conservation**: `PHASTCONS_470_WAY`, `CACTUS_241_WAY`
#' * **Indel**: `IS_INSERTION`, `IS_DELETION`
#'
#' By default the Expression and Accessibility groups are shown. Other
#' features can be selected by column name or by group label, for example
#' `features = c("Accessibility", "Splicing")`, and `features = NULL`
#' shows all features. Only these 18 columns are treated as feature
#' attributions; other numeric columns, such as `AVI_SCORE`, are ignored
#' unless requested explicitly by name. The total impact score
#' (the raw AVI score, which equals the sum of all feature attributions) can
#' be requested with `features = "Total"`. In the "stacked" plot type the
#' total is drawn as points over the stacked feature contributions (as on the
#' AlphaGenome Atlas website) rather than being stacked itself. Note that the raw AVI score is on the model logit
#' scale; the PHRED-scaled AVI score shown on the AlphaGenome Atlas website can
#' be retrieved separately with `scorer = "AVI_SCORE"` in
#' [LoadAVIFromAtlas()].
#'
#' The Atlas contains three SNVs at each genomic position (one per alternate
#' allele). Following the AlphaGenome Atlas publication, the variant with the
#' largest absolute total AVI score is kept at each position. The region is
#' then divided into `bins` equal-width bins and each bin shows the value of
#' largest magnitude among its positions, so that sharp SHAP peaks are
#' preserved while large regions remain fast to plot. Positions without a
#' variant contribute a value of zero.
#'
#' @param avi AVI SHAP feature attributions. Either a `data.frame` with one
#' row per variant containing chromosome, position, ref and alt columns along
#' with one numeric column per feature (for example as returned by
#' [LoadAVIFromAtlas()]), or `TRUE` to query the AlphaGenome Atlas API for the
#' region (see [LoadAVIFromAtlas()]; requires the **reticulate** package, the
#' `alphagenome` Python package, and an API key). If no variants fall in the
#' region an empty track is drawn.
#' @param region Genomic region ([GenomicRanges::GRanges] or a string that can
#' be converted to `GRanges` like "chr11:5225000-5230000")
#' @param features Features to plot, given as feature column names and/or
#' modality group labels (see Details), or "Total" for the total
#' AVI score. Groups without any feature in the data are skipped. Set to `NULL`
#' to plot all features. Features are stacked and listed in the legend in the
#' order given.
#' @param group.features Combine features into modality groups (see Details)
#' before plotting. If `FALSE`, each feature is
#' shown separately.
#' @param type Plot type. One of "stacked" (stacked bars showing the signed
#' contribution of each feature group) or "heatmap" (one row per feature
#' group, colored by SHAP value).
#' @param bins Number of bars (bins) to display across the region. Each bin
#' is `width(region) / bins` bp wide.
#' @param ymax Maximum absolute y-axis value for the "stacked" plot type. If
#' `NULL`, determined from the data.
#' @param colors Named vector of colors for feature groups (or features,
#' if `group.features = FALSE`). Groups without an entry use the default
#' palette.
#' @param y_label Y-axis label. If `NULL`, "AVI SHAP" is used, or "AVI score
#' (raw)" when only the total AVI score is plotted.
#' @param show.axis Show x-axis (default: TRUE)
#' @param api.key AlphaGenome API key, used when `avi = TRUE`. If
#' `NULL`, the `ALPHA_GENOME_API_KEY` environment variable is used.
#' @param max.workers Number of concurrent requests to the Atlas API when
#' `avi = TRUE` (see [LoadAVIFromAtlas()]).
#' @param max.retries Maximum number of retries per throttled Atlas API
#' request when `avi = TRUE` (see [LoadAVIFromAtlas()]).
#' @param verbose Show Atlas query progress messages when `avi = TRUE`.
#'
#' @return Returns a [ggplot2::ggplot()] object
#'
#' @references Cheng et al. (2026). AlphaGenome Atlas: in silico mutagenesis of
#' the entire human genome improves prioritization and interpretation of
#' non-coding variants
#' <https://deepmind.google.com/science/alphagenome/atlas>
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_hline geom_tile scale_fill_manual scale_fill_gradient2 scale_shape_manual guides guide_legend coord_cartesian xlab ylab theme element_blank
#' @importFrom GenomicRanges GRanges start end width
#' @importFrom Seqinfo seqnames
#' @importFrom scales hue_pal
#' @importFrom stats setNames
#' @export
#' @concept visualization
#' @concept alphagenome
#' @examples
#' # simulated AVI SHAP feature attributions for three positions
#' avi <- data.frame(
#'   chromosome = "chr1",
#'   position = rep(c(713600, 713800, 714200), each = 3),
#'   ref = "A",
#'   alt = rep(c("C", "G", "T"), times = 3),
#'   MAX_ABS_ATAC = runif(9, 0, 1),
#'   MERGED_SPLICING = runif(9, 0, 2),
#'   CACTUS_241_WAY = runif(9, -0.5, 0.5)
#' )
#' # expression and accessibility features (default)
#' AVITrack(avi = avi, region = "chr1:713500-714500")
#' # all features, one row per modality group
#' AVITrack(avi = avi, region = "chr1:713500-714500", features = NULL,
#'          type = "heatmap")
#' # feature contributions with the total AVI score overlaid
#' AVITrack(avi = avi, region = "chr1:713500-714500",
#'          features = c("Accessibility", "Splicing", "Total"))
#'
#' \dontrun{
#' # query the AlphaGenome Atlas API (requires ALPHA_GENOME_API_KEY)
#' AVITrack(avi = TRUE, region = "chr11:5225000-5230000")
#' }
AVITrack <- function(
  avi,
  region,
  features = c("Expression", "Accessibility", "Total"),
  group.features = TRUE,
  type = "stacked",
  bins = 200,
  ymax = NULL,
  colors = NULL,
  y_label = NULL,
  show.axis = TRUE,
  api.key = NULL,
  max.workers = 4,
  max.retries = 5,
  verbose = FALSE
) {
  possible.types <- c("stacked", "heatmap")
  if (!(type %in% possible.types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible.types, collapse = ", ")
    )
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- GRanges(region)
  }
  if (length(x = region) != 1) {
    stop("region must describe a single genomic range")
  }
  if (isTRUE(x = avi)) {
    avi <- LoadAVIFromAtlas(
      region = region, api.key = api.key, max.workers = max.workers,
      max.retries = max.retries, verbose = verbose
    )
  } else if (inherits(x = avi, what = "data.frame")) {
    avi <- as.data.frame(x = avi)
    avi <- CheckAVIColumns(avi = avi)
    keep <- avi[["chromosome"]] == as.character(x = seqnames(x = region)) &
      avi[["position"]] >= start(x = region) &
      avi[["position"]] <= end(x = region)
    avi <- avi[keep, , drop = FALSE]
  } else {
    stop("avi must be a data.frame or TRUE to query the Atlas API")
  }
  if (nrow(x = avi) == 0) {
    message("No AVI data found in region; drawing an empty track")
  }
  # the SHAP attributions sum to the raw AVI score
  all.features <- AVIFeatureColumns(avi = avi, allow.empty = TRUE)
  if (length(x = all.features) > 0) {
    avi[["Total"]] <- rowSums(
      x = as.matrix(x = avi[, all.features, drop = FALSE]), na.rm = TRUE
    )
  }
  features <- ResolveAVIFeatures(avi = avi, features = features)
  collapsed <- CollapseAVIAlleles(avi = avi, features = features)
  # combine contributions of features in the same group
  if (group.features) {
    feature.groups <- AVIFeatureGroups(features = features)
    tracks <- unique(x = unname(obj = feature.groups))
    for (track in tracks) {
      members <- features[feature.groups == track]
      collapsed[[track]] <- rowSums(
        x = as.matrix(x = collapsed[, members, drop = FALSE]), na.rm = TRUE
      )
    }
  } else {
    tracks <- features
  }
  # groups are ordered as the features were requested (or as the columns
  # appear in the data when features = NULL)
  group.levels <- unique(x = tracks)
  # the attributions are labelled "AVI SHAP" as in the Atlas publication; on
  # its own the total is the raw AVI score
  only.total <- identical(x = group.levels, y = "Total")
  value.label <- if (only.total) "AVI score (raw)" else "AVI SHAP"
  y_label <- y_label %||% value.label

  # bin, keeping the value of largest magnitude in each bin
  long <- BinAVI(
    values = collapsed,
    region = region,
    tracks = tracks,
    bins = bins,
    rank.by = if (type == "stacked") "Total" else NULL
  )
  bin.size <- attr(x = long, which = "bin.size")
  long[["group"]] <- factor(x = long[["group"]], levels = group.levels)
  feature.table <- AVIFeatureTable()
  defaults <- c(
    setNames(object = feature.table[["color"]], nm = feature.table[["group"]]),
    "Total" = "#000000"
  )
  colors <- c(colors, defaults)
  colors <- colors[!duplicated(x = names(x = colors))]
  missing.groups <- setdiff(x = group.levels, y = names(x = colors))
  if (length(x = missing.groups) > 0) {
    colors <- c(
      colors,
      setNames(
        object = hue_pal()(length(x = missing.groups)),
        nm = missing.groups
      )
    )
  }
  colors <- colors[group.levels]

  chromosome <- as.character(x = seqnames(x = region))
  xlim <- c(start(x = region), end(x = region))
  if (type == "stacked") {
    ylim <- NULL
    if (!is.null(x = ymax)) {
      ylim <- c(min(0, -ymax), ymax)
    }
    # the total is drawn over the stacked components rather than stacked
    has.total <- "Total" %in% group.levels && length(x = group.levels) > 1
    components <- long
    if (has.total) {
      components <- long[long[["group"]] != "Total", , drop = FALSE]
      total <- long[long[["group"]] == "Total", , drop = FALSE]
    }
    p <- ggplot(
      data = components,
      mapping = aes(
        x = .data[["position"]], y = .data[["value"]], fill = .data[["group"]]
      )
    ) +
      geom_col(position = "stack", width = bin.size) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
      scale_fill_manual(values = colors, name = "Feature") +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      theme_browser(axis.text.y = TRUE)
    if (has.total) {
      p <- p +
        geom_point(
          data = total,
          mapping = aes(
            x = .data[["position"]], y = .data[["value"]], shape = "Total"
          ),
          inherit.aes = FALSE,
          color = colors[["Total"]],
          size = 0.4
        ) +
        scale_shape_manual(values = c("Total" = 18), name = NULL) +
        guides(shape = guide_legend(override.aes = list(size = 2)))
    }
  } else {
    long[["group"]] <- factor(
      x = long[["group"]], levels = rev(x = group.levels)
    )
    p <- ggplot(
      data = long,
      mapping = aes(
        x = .data[["position"]], y = .data[["group"]], fill = .data[["value"]]
      )
    ) +
      geom_tile(width = bin.size) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        name = value.label
      ) +
      coord_cartesian(xlim = xlim) +
      theme_browser(axis.text.y = TRUE)
  }
  p <- p +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label)
  if (!show.axis) {
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    )
  }
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
#' @return Returns a [ggplot2::ggplot()] object
#' @importFrom SeuratObject DefaultAssay
#' @importFrom ggplot2 xlab
#' @importFrom Seqinfo seqnames
#'
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package = "Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' TilePlot(object = atac_small, region = c("chr1:713500-714500"))
#' }
TilePlot <- function(
  object,
  region,
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
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("Requested assay is not a ChromatinAssay5.")
  }
  region <- FindRegion(
    object = object,
    region = region,
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
  smoothed$bin <- (smoothed$bin - 1) * window +
    as.numeric(x = colnames(x = cutmatrix)[[1]])
  return(smoothed)
}

#' @importFrom ggplot2 ggplot aes geom_raster ylab scale_fill_gradient scale_y_reverse guides guide_legend geom_hline
CreateTilePlot <- function(df, n, legend = TRUE) {
  # create plot
  p <- ggplot(
    data = df,
    aes(x = .data[["bin"]], y = .data[["idx"]], fill = .data[["value"]])
  ) +
    facet_wrap(
      facets = ~group,
      scales = "free_y",
      ncol = 1,
      strip.position = "left"
    ) +
    geom_raster() +
    theme_browser(legend = legend) +
    geom_hline(yintercept = c(0, n), linewidth = 0.1) +
    ylab(paste0("Fragments (", n, " cells)")) +
    scale_fill_gradient(low = "white", high = "darkred") +
    scale_y_reverse() +
    guides(fill = guide_legend(
      title = "Fragment\ncount",
      keywidth = 1 / 2, keyheight = 1
    )) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    )
  return(p)
}

#' Genome browser theme
#'
#' Theme applied to panels in the [CoveragePlot()] function.
#'
#' @param ... Additional arguments
#' @param legend Display plot legend
#' @param axis.text.y Display y-axis text
#'
#' @importFrom ggplot2 theme theme_classic element_blank element_text
#' @export
#' @concept visualization
#' @examples
#' \donttest{
#' PeakPlot(atac_small, region = "chr1:710000-715000") + theme_browser()
#' }
theme_browser <- function(..., legend = TRUE, axis.text.y = FALSE) {
  browser.theme <- theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0)
    )
  if (!axis.text.y) {
    browser.theme <- browser.theme +
      theme(
        axis.text.y = element_blank()
      )
  }
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
    starts <- starts[starts <= df$end]
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    tx_id = df$tx_id[[1]],
    gene_name = df$gene_name[[1]],
    gene_biotype = df$gene_biotype[[1]],
    type = "arrow"
  )
  return(breaks)
}

reformat_annotations <- function(
  annotation,
  start.pos,
  end.pos,
  collapse_transcript = TRUE
) {
  total.width <- end.pos - start.pos
  tick.freq <- total.width / 8
  annotation <- annotation[annotation$type == "exon"]
  exons <- as.data.frame(x = annotation, row.names = NULL)
  if (collapse_transcript) {
    annotation <- split(
      x = annotation,
      f = annotation$gene_name
    )
  } else {
    annotation <- split(
      x = annotation,
      f = annotation$tx_id
    )
  }
  annotation <- lapply(X = annotation, FUN = as.data.frame, row.names = NULL)

  # add gene total start / end
  gene_bodies <- list()
  for (i in seq_along(annotation)) {
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      tx_id = annotation[[i]]$tx_id[[1]],
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
    breaks <- split_body(df = df, width = tick.freq)
    df <- rbind(df, breaks)
    gene_bodies[[i]] <- df
  }
  gene_bodies <- do.call(what = rbind, args = gene_bodies)

  # record if genes overlap
  overlap_idx <- record_overlapping(
    annotation = gene_bodies,
    min.gapwidth = 1000,
    collapse_transcript = collapse_transcript
  )
  if (collapse_transcript) {
    gene_bodies$dodge <- overlap_idx[gene_bodies$gene_name]
    exons$dodge <- overlap_idx[exons$gene_name]
  } else {
    gene_bodies$dodge <- overlap_idx[gene_bodies$tx_id]
    exons$dodge <- overlap_idx[exons$tx_id]
  }

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
record_overlapping <- function(
  annotation,
  min.gapwidth = 1000,
  collapse_transcript = TRUE
) {
  # convert back to granges
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
  if (collapse_transcript) {
    names(x = idx) <- gr$gene_name
  } else {
    names(x = idx) <- gr$tx_id
  }
  return(idx)
}

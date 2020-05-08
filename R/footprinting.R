#' @include generics.R
#'
NULL

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

#' Extract footprint data for a set of transcription factors
#'
#' @param object A Seurat object
#' @param features A vector of features to plot
#' @param assay Name of assay to use
#' @param group.by A grouping variable
#' @param idents Set of identities to include in the plot
#' @export
#' @importFrom Seurat DefaultAssay
GetFootprintData <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  positionEnrichment <- GetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment"
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  all.groups <- unique(x = obj.groups)
  plot.data <- lapply(X = features, FUN = function(x) {
    if (!(x %in% names(x = positionEnrichment))) {
      warning("Footprint information for ", x, " not found in assay")
      return()
    } else {
      fp <- positionEnrichment[[x]]
      # extract expected
      expected <- fp["expected", ]
      # remove expected and motif position from main matrix
      fp <- fp[1:(nrow(x = fp) - 2), ]
      bg.norm <- lapply(X = all.groups, FUN = function(x) {
        cells.use <- names(x = obj.groups)[obj.groups == x]
        mat.use <- fp[cells.use, ]
        return(BackgroundMeanNorm(x = mat.use, background = 50))
      })
      bg.norm <- do.call(what = rbind, args = bg.norm)
      groupmeans <- ApplyMatrixByGroup(
        mat = bg.norm,
        groups = obj.groups,
        fun = colMeans,
        normalize = FALSE
      )
      # add feature information
      groupmeans$feature <- x
      groupmeans$class <- "Observed"

      # add expected insertions
      expect.df <- data.frame(
        group = NA,
        count = expected,
        norm.value = expected,
        position = as.numeric(x = names(x = expected)),
        feature = x,
        class = "Expected"
      )
      groupmeans <- rbind(groupmeans, expect.df)
      return(groupmeans)
    }
  })
  plot.data <- do.call(what = rbind, args = plot.data)
  return(plot.data)
}

#' @param regions A set of genomic ranges containing the motif instances
#' @param genome A \code{\link[BSgenome]{BSgenome}} object
#' @param motif.name Name of a motif stored in the assay to footprint. If not
#' supplied, must supply a set of regions.
#' @param key Key to store positional enrichment information under.
#' @param upstream Number of bases to extend upstream
#' @param downstream Number of bases to extend downstream
#' @param verbose Display messages
#' @param compute.expected Find the expected number of insertions at each
#' position given the local DNA sequence context and the insertion bias of Tn5
#' @param ... Arguments passed to other methods
#' @importFrom IRanges width
#' @export
#' @rdname Footprint
#' @method Footprint ChromatinAssay
Footprint.ChromatinAssay <- function(
  object,
  genome,
  motif.name = NULL,
  key = motif.name,
  regions = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  compute.expected = TRUE,
  verbose = TRUE,
  ...
) {
  if (is.null(x = motif.name) & is.null(x = regions)) {
    stop("Must supply the name of a motif or a set of regions")
  } else if (!is.null(x = motif.name) & !is.null(x = regions)) {
    stop("Supplied both a motif name and set of regions. Choose one only.")
  } else if (!is.null(x = motif.name)) {
    # pull motif positions from object
    motif.obj <- Motifs(object = object)
    if (length(x = motif.name) != length(x = key)) {
      stop("A Key needs to be supplied for each motif")
    }
    regionlist <- lapply(X = motif.name, FUN = function(x) {
      GetFootprintRegions(motif.obj = motif.obj, motif.name = x)
    })
  } else {
    if (is.null(x = key)) {
      stop("Must set a key to store positional enrichment information")
    }
    # supplied regions, put into list
    regionlist <- list(regions)
  }
  for (i in seq_along(along.with = regionlist)) {
    object <- RunFootprint(
      object = object,
      genome = genome,
      regions = regionlist[[i]],
      key = key[[i]],
      upstream = upstream,
      downstream = downstream,
      compute.expected = compute.expected,
      verbose = verbose
    )
  }
  return(object)
}

# Extract regions for a given TF name
GetFootprintRegions <- function(
  motif.obj,
  motif.name
) {
  motif.positions <- GetMotifData(object = motif.obj, slot = "positions")
  if (is.null(x = motif.positions)) {
    stop("Motif positions not present in Motif object.")
  } else {
    if (motif.name %in% names(x = motif.positions)) {
      regions <- motif.positions[[motif.name]]
    } else {
      # convert to common name and look up
      common.names <- GetMotifData(object = motif.obj, slot = "motif.names")
      if (motif.name %in% common.names) {
        motif.idx <- names(x = which(x = common.names == motif.name))
        regions <- motif.positions[[motif.idx]]
      } else {
        stop("Motif not found")
      }
    }
    return(regions)
  }
}

#' @rdname Footprint
#' @param assay Name of assay to use
#' @method Footprint Seurat
#' @export
#' @importFrom Seurat DefaultAssay
Footprint.Seurat <- function(
  object,
  genome,
  regions = NULL,
  motif.name = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  object[[assay]] <- Footprint(
    object = object[[assay]],
    regions = regions,
    motif.name = motif.name,
    genome = genome,
    group.by = group.by,
    idents = idents,
    upstream = upstream,
    downstream = downstream,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @param genome A BSgenome object
#' @param region Region to use when assessing bias. Default is human chromosome 1.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings oligonucleotideFrequency
#' @export
#' @rdname InsertionBias
#' @method InsertionBias ChromatinAssay
#' @examples
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm10)
#'
#' region.use <- GRanges(
#'   seqnames = c('chr1', 'chr2'),
#'   IRanges(start = c(1,1), end = c(195471971, 182113224))
#' )
#'
#' InsertionBias(
#'  object = atac_small,
#'  genome = BSgenome.Mmusculus.UCSC.mm10,
#'  region = region.use
#' )
#' }
InsertionBias.ChromatinAssay <- function(
  object,
  genome,
  region = 'chr1-1-249250621',
  verbose = TRUE,
  ...
) {
  reads <- MultiGetReadsInRegion(
    object = object,
    region = region,
    ...
  )
  insertions <- GRanges(
    seqnames = c(reads$chr, reads$chr),
    ranges = IRanges(
      start = c(reads$start, reads$end),
      width = 1
    ),
    strand = '+'
  )
  insertions <- Extend(x = insertions, upstream = 3, downstream = 2)
  sequences <- as.vector(x = getSeq(x = genome, insertions))
  seq.freq <- table(sequences)
  # remove sequences containing N
  keep.seq <- !grepl(pattern = "N", x = names(x = seq.freq))
  insertion_hex_freq <- as.matrix(x = seq.freq[keep.seq])
  genome_freq <- oligonucleotideFrequency(
    x = getSeq(x = genome, names = 'chr1'),
    width = 6
  )
  if (nrow(x = insertion_hex_freq) != length(x = genome_freq)) {
    stop("Not all hexamers represented in input region")
  }
  insertion_hex_freq <- insertion_hex_freq[names(x = genome_freq), ]
  bias <- insertion_hex_freq / genome_freq
  object <- SetAssayData(object = object, slot = "bias", new.data = bias)
  return(object)
}

#' @param assay Name of assay to use
#' @rdname InsertionBias
#' @method InsertionBias Seurat
#' @importFrom Seurat DefaultAssay
#' @export
InsertionBias.Seurat <- function(
  object,
  genome,
  assay = NULL,
  region = 'chr1-1-249250621',
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  object[[assay]] <- InsertionBias(
    object = object[[assay]],
    genome = genome,
    region = region,
    verbose = verbose,
    ...
  )
  return(object)
}

####### Not exported #######

# Divide matrix by flanks
#' @importMethodsFrom Matrix mean
BackgroundMeanNorm <- function(x, background = 50) {
  positions.use <- c(1:background, (ncol(x = x) - background):ncol(x = x))
  flanks <- mean(x = x[ ,positions.use])
  x <- x / flanks
  return(x)
}

# Find the expected number of insertions over a set of genomic regions
# given the DNA sequences and the insertion bias of Tn5 for the experiment
# @param dna.sequence A set of DNA sequences
# @param bias A vector describing Tn5 insertion frequency at each hexamer
# @param verbose Display messages
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#' @importFrom IRanges width narrow
FindExpectedInsertions <- function(dna.sequence, bias, verbose = TRUE) {
  if (verbose) {
    message("Computing base composition at motif sites")
  }
  total.hexamer.positions <- width(x = dna.sequence)[[1]] - 6
  hex.key <- seq_along(along.with = bias)
  names(hex.key) <- names(bias)

  # x is the hexamer frequency
  x <- vector(
    mode = "numeric",
    length = length(x = bias) * total.hexamer.positions
  )
  # i is the hexamers sequence
  i <- vector(
    mode = "numeric",
    length = length(x = bias) * total.hexamer.positions
  )
  # j is the base position
  j <- vector(
    mode = "numeric",
    length = length(x = bias) * total.hexamer.positions
  )
  current.pos <- 1

  for (jj in seq_len(length.out = total.hexamer.positions)) {
    # resize dna string set
    resized <- narrow(x = dna.sequence, start = jj, width = 6)
    resized <- as.character(x = resized)
    # need to remove any that contain N
    resized <- resized[!grepl(pattern = "N", x = resized)]
    # count
    frequencies <- table(resized)
    end.pos <- current.pos + length(x = frequencies) - 1
    # append
    x[current.pos:end.pos] <- as.numeric(x = frequencies)
    j[current.pos:end.pos] <- jj
    i[current.pos:end.pos] <- as.vector(x = hex.key[names(x = frequencies)])
    # shift current position
    current.pos <- end.pos + 1
  }
  # trim vectors
  x <- x[1:(current.pos - 1)]
  i <- i[1:(current.pos - 1)]
  j <- j[1:(current.pos - 1)]

  # construct matrix
  hexamer.matrix <- sparseMatrix(i = i, j = j, x = x)
  rownames(hexamer.matrix) <- names(x = hex.key)
  colnames(hexamer.matrix) <- seq_len(length.out = total.hexamer.positions)
  hexamer.matrix <- as.matrix(x = hexamer.matrix)

  if (verbose) {
    message("Computing expected Tn5 insertions per base")
  }
  # ensure correct order
  hexamer.matrix <- hexamer.matrix[names(x = bias), ]
  expected.insertions <- as.vector(
    x = crossprod(x = hexamer.matrix, y = as.matrix(x = bias))
  )

  # normalize expected by dividing by flanks
  # TODO use BackgroundMeanNorm function here
  flanks <- mean(
    x = c(expected.insertions[1:50],
          expected.insertions[
            (total.hexamer.positions - 50):total.hexamer.positions
            ])
  )
  expected.insertions <- expected.insertions / flanks
  return(expected.insertions)
}

# Get size of motif that was footprinted
# @param object A Seurat object
# @param feature A vector of footprinted TFs
# @param assay Name of assay to use
#' @importFrom Seurat DefaultAssay GetAssayData
GetMotifSize <- function(
  object,
  features,
  assay = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  positionEnrichment <- GetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment"
  )
  sizes <- c()
  for (i in features) {
    motif <- positionEnrichment[[i]]["motif", ]
    sizes <- c(sizes, sum(motif))
  }
  names(x = sizes) <- features
  return(sizes)
}

# @param object A ChromatinAssay object
# @param genome a BSgenome object
# @param regions a set of regions of the same width
# @param compute.expected Compute the expected number of insertions given
# DNA sequence bias of Tn5
# @param verbose Display messages
#' @importFrom Seurat GetAssayData
#' @importFrom Biostrings getSeq
Pileup <- function(
  object,
  genome,
  regions,
  compute.expected = TRUE,
  verbose = TRUE
) {
  # add three bases each side here so we can get the hexamer frequencies
  # for every position
  dna.sequence <- getSeq(x = genome, Extend(
    x = regions,
    upstream = 3,
    downstream = 3
  )
  )
  if (compute.expected) {
    bias <- GetAssayData(object = object, slot = "bias")
    if (is.null(x = bias)) {
      stop("Insertion bias not computed")
    } else {
      expected.insertions <- FindExpectedInsertions(
        dna.sequence = dna.sequence,
        bias = bias,
        verbose = verbose
      )
    }
  } else {
    expected.insertions <- rep(1, width(x = dna.sequence)[[1]])
  }

  if (verbose) {
    message("Computing observed Tn5 insertions per base")
  }
  # count insertions at each position for each cell
  insertion.matrix <- CreateRegionPileupMatrix(
    object = object,
    regions = regions
  )

  # store expected as one additional row in the matrix
  expected.insertions <- t(x = as.matrix(x = expected.insertions))
  rownames(x = expected.insertions) <- "expected"
  insertion.matrix <- rbind(insertion.matrix, expected.insertions)
  return(insertion.matrix)
}

# Run Footprint function for a single set of regions
# @param object A ChromatinAssay object
# @param genome A BSgenome object
# @param regions A set of genomic regions
# @param key Key to store positional enrichment matrix under
# @param upstream Number of bases to extend upstream
# @param downstream Number of bases to extend downstream
#' @importFrom IRanges width
#' @importFrom BiocGenerics sort
RunFootprint <- function(
  object,
  genome,
  key,
  regions,
  upstream = 250,
  downstream = 250,
  compute.expected = TRUE,
  verbose = TRUE
) {
  motif.size <- width(x = regions)[[1]]
  regions <- sort(x = regions)
  # extend upstream and downstream
  regions <- Extend(
    x = regions,
    upstream = upstream,
    downstream = downstream
  )
  if (compute.expected) {
    # check that bias is computed
    bias <- GetAssayData(object = object, slot = "bias")
    if (is.null(x = bias)) {
      if (verbose) {
        message("Computing Tn5 insertion bias")
      }
      object <- InsertionBias(
        object = object,
        genome = genome
      )
    }
  }
  # find expected and observed insertions across all regions
  insertion.matrix <- Pileup(
    object = object,
    genome = genome,
    regions = regions,
    compute.expected = compute.expected,
    verbose = verbose
  )
  # encode motif position as additional row in matrix
  motif.vec <- t(x = matrix(
    data = c(
      rep(x = 0, upstream),
      rep(x = 1, motif.size),
      rep(x = 0, downstream)
    )
  )
  )
  rownames(x = motif.vec) <- "motif"
  insertion.matrix <- rbind(insertion.matrix, motif.vec)
  # store results in the assay
  object <- suppressWarnings(expr = SetAssayData(
    object = object,
    slot = "positionEnrichment",
    new.data = insertion.matrix,
    key = key
  ))
  return(object)
}

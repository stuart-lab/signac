#' @include generics.R
#'
NULL

#' Plot footprinting results
#'
#' @param object A Seurat object
#' @param features A vector of features to plot
#' @param assay Name of assay to use
#' @param group.by A grouping variable
#' @param idents Set of identities to include in the plot
#' @export
#' @importFrom Seurat DefaultAssay
#' @importFrom ggplot2 ggplot aes
PlotFootprint <- function(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  # TODO add option to plot expected below plot
  # TODO add option to label top identities
  # TODO add option to show variance among cells
  plot.data <- GetFootprintData(
    object = object,
    features = features,
    assay = assay,
    group.by = group.by,
    idents = idents
  )
  p <- ggplot(
    data = plot.data,
    mapping = aes(x = position, y = norm.value, color = group)
  ) +
    geom_line(stat = "identity", size = 0.2) +
    facet_wrap(facets = ~feature)
  p <- p +
    xlab("Distance from motif") +
    ylab(label = "Normalized Tn5 insertions") +
    theme_minimal()
  return(p)
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
  # TODO add normalization options
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
  plot.data <- lapply(X = features, FUN = function(x) {
    if (!(x %in% names(x = positionEnrichment))) {
      warning("Footprint information for ", x, " not found in assay")
      return()
    } else {
      fp <- positionEnrichment[[x]]
      # remove row containing expected insertions
      expected <- fp["expected", ]
      fp <- fp[1:(nrow(x = fp) - 1), ]
      # average the signal per group per base
      groupmeans <- ApplyMatrixByGroup(
        mat = fp,
        groups = obj.groups,
        fun = colMeans,
        normalize = FALSE
      )
      # add feature information
      groupmeans$feature <- x
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
          regions = motif.positions[[motif.idx]]
        } else {
          stop("Motif not found")
        }
      }
    }
  } else {
    if (is.null(x = key)) {
      stop("Must set a key to store positional enrichment information")
    }
  }
  regions <- sort(x = regions)
  # extend upstream and downstream
  regions <- Extend(
    x = regions,
    upstream = upstream,
    downstream = downstream
  )
  # find expected and observed insertions across all regions
  insertion.matrix <- Pileup(
    object = object,
    genome = genome,
    regions = regions,
    compute.expected = compute.expected,
    verbose = verbose
  )
  # store results in the assay
  object <- suppressWarnings(expr = SetAssayData(
    object = object,
    slot = "positionEnrichment",
    new.data = insertion.matrix,
    key = key
  ))
  return(object)

}

#' @rdname Footprint
#' @param assay Name of assay to use
#' @method Footprint Seurat
#' @importFrom Seurat DefaultAssay
Footprint.Seurat <- function(
  object,
  genome,
  regions = NULL,
  motif.name = NULL,
  group.by = NULL,
  idents = NULL,
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
#' \donttest{
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
  flanks <- mean(
    x = c(expected.insertions[1:50],
          expected.insertions[
            (total.hexamer.positions - 50):total.hexamer.positions
            ])
  )
  expected.insertions <- expected.insertions / flanks
  return(expected.insertions)
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
      if (verbose) {
        message("Computing Tn5 insertion bias")
      }
      object <- InsertionBias(
        object = object,
        genome = genome
      )
    }
    expected.insertions <- FindExpectedInsertions(
      dna.sequence = dna.sequence,
      bias = bias,
      verbose = verbose
    )
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


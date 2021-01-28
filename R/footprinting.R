#' @include generics.R
#'
NULL

#' Get footprinting data
#'
#' Extract footprint data for a set of transcription factors or metafeatures.
#' This function will pull accessibility data for a given feature (eg, a TF),
#' and perform background normalization for each identity class. This is the
#' data that's used to create TF footprinting plots with the
#' \code{PlotFootprint} function.
#'
#' @param object A Seurat object
#' @param features A vector of features to extract data for
#' @param assay Name of assay to use
#' @param group.by A grouping variable
#' @param idents Set of identities to group cells by
#' @export
#' @return Returns a matrix
#' @concept footprinting
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
  levels.stash <- levels(x = obj.groups)
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
  plot.data$group <- factor(x = plot.data$group, levels = levels.stash)
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
#' @param in.peaks Restrict motifs to those that fall in peaks
#' @param ... Arguments passed to other methods
#' @importFrom IRanges width
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @export
#' @concept footprinting
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
  in.peaks = FALSE,
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
  # run in parallel
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  matrices <- mylapply(
    X = seq_along(along.with = regionlist),
    FUN = function(x) {
      RunFootprint(
        object = object,
        genome = genome,
        regions = regionlist[[x]],
        upstream = upstream,
        downstream = downstream,
        compute.expected = compute.expected,
        in.peaks = in.peaks,
        verbose = verbose
      )
    })

  # store in object
  for (i in seq_along(along.with = matrices)) {
    object <- SetAssayData(
      object = object,
      slot = "positionEnrichment",
      new.data = matrices[[i]],
      key = key[[i]]
    )
  }
  return(object)
}

#' @rdname Footprint
#' @param assay Name of assay to use
#' @method Footprint Seurat
#' @export
#' @concept footprinting
#' @importFrom Seurat DefaultAssay
Footprint.Seurat <- function(
  object,
  genome,
  regions = NULL,
  motif.name = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  in.peaks = FALSE,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  object[[assay]] <- Footprint(
    object = object[[assay]],
    regions = regions,
    motif.name = motif.name,
    genome = genome,
    upstream = upstream,
    downstream = downstream,
    in.peaks = in.peaks,
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
#' @concept footprinting
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
  chr.use <- unlist(x = strsplit(x = region, split = "-", fixed = TRUE))[[1]]
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
    x = getSeq(x = genome, names = chr.use),
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
#' @concept footprinting
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

# Extract regions for a given TF name
# @param motif.obj A Motif object
# @param motif.name Name of a motif to pull positional information for
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

# Run Footprint function for a single set of regions
# @param object A ChromatinAssay object
# @param genome A BSgenome object
# @param regions A set of genomic regions
# @param upstream Number of bases to extend upstream
# @param downstream Number of bases to extend downstream
# @param in.peaks Restrict to motifs in peaks
#' @importFrom IRanges width subsetByOverlaps
#' @importFrom Biostrings getSeq
#' @importFrom BiocGenerics sort
RunFootprint <- function(
  object,
  genome,
  regions,
  upstream = 250,
  downstream = 250,
  compute.expected = TRUE,
  in.peaks = FALSE,
  verbose = TRUE
) {
  motif.size <- width(x = regions)[[1]]
  regions <- sort(x = regions)
  if (in.peaks) {
    regions <- subsetByOverlaps(x = regions, ranges = granges(x = object))
  }
  # extend upstream and downstream
  regions <- Extend(
    x = regions,
    upstream = upstream,
    downstream = downstream
  )
  if (verbose) {
    message("Computing observed Tn5 insertions per base")
  }
  if (compute.expected) {
    bias <- GetAssayData(object = object, slot = "bias")
    if (is.null(x = bias)) {
      stop("Insertion bias not computed")
    } else {
      # add three bases each side here so we can get the hexamer frequencies
      # for every position
      dna.sequence <- getSeq(x = genome, Extend(
        x = regions,
        upstream = 3,
        downstream = 3
       )
      )
      expected.insertions <- FindExpectedInsertions(
        dna.sequence = dna.sequence,
        bias = bias,
        verbose = verbose
      )
    }
  } else {
    expected.insertions <- rep(1, width(x = dna.sequence)[[1]])
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
  return(insertion.matrix)
}

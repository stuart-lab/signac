#' @include generics.R
#'
NULL

#' Get footprinting data
#'
#' Extract footprint data for a set of transcription factors or metafeatures.
#' This function will pull accessibility data for a given feature (eg, a TF),
#' and perform background normalization for each identity class. This is the
#' data that's used to create TF footprinting plots with the
#' `PlotFootprint` function.
#'
#' @param object A Seurat object
#' @param features A vector of features to extract data for
#' @param assay Name of assay to use
#' @param group.by A grouping variable
#' @param idents Set of identities to group cells by
#' @export
#' @return Returns a data.frame with the following columns:
#' \itemize{
#'   \item{group: Cell group (determined by group.by parameter}
#'   \item{position: Position relative to motif center}
#'   \item{count: Normalized Tn5 insertion counts at each position}
#'   \item{norm.value: Normalized Tn5 insertion counts at each position (same as count)}
#'   \item{feature: Name of the footprinted motif}
#'   \item{class: observed or expected}
#'  }
#' @concept footprinting
#' @importFrom SeuratObject DefaultAssay
GetFootprintData <- function(
    object,
    features,
    assay = NULL,
    group.by = NULL,
    idents = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("The requested assay is not a ChromatinAssay5")
  }
  
  region.enrichment <- RegionAggr(object[[assay]])
  # get existing features 
  region.enrichment.names <- vapply(region.enrichment, FUN = function(x) x@name, FUN.VALUE = character(1))

  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  levels.stash <- levels(x = obj.groups)
  all.groups <- unique(x = obj.groups)
  
  # Loop through features to plot
  plot.data <- lapply(X = features, FUN = function(feature) {
    if (!(feature %in% region.enrichment.names)) {
      warning("Footprint information for ", feature, " not found in assay")
      return()
    } 
    # V2 edit: 
    agg.idx <- which(region.enrichment.names == feature) # could be a list of more than one 
    region.agg.list <- region.enrichment[agg.idx] # get a list of agg.obj
    
    if (length(region.agg.list)>1){
      # assert same motif size 
      w <- unique(unlist(lapply(region.agg.list, function(x) width(x@regions))))
      if (length(w) != 1) {
        warning("Feature ",feature, " have different motif width! Skipping")
        return(NULL)
      }
      motif.width <- w[[1]]
      # ensure upstream & downstream matches
      matrix.upstream <- vapply(region.agg.list, function(x) x@upstream, integer(1))
      target.upstream <- min(matrix.upstream)
      
      matrix.downstream <- vapply(region.agg.list, function(x) x@downstream, integer(1))
      target.downstream <- min(matrix.downstream)
      
      if (any(matrix.downstream != target.downstream ) || any(matrix.upstream != target.upstream)){
         warning("Truncating matrices to the smallest width to align Footprint for ",feature)
      }
      
      all.mat <- lapply(region.agg.list, function(x){
        mat <- x@matrix
        start_col <- x@upstream - target.upstream + 1
        end_col <- ncol(mat) - x@downstream + target.downstream
        if (start_col <1 || end_col > ncol(mat) || start_col > end_col){
          warning("Skipping incompatible RegionAggregation for ", feature)
          return(NULL)
        }
        mat <- mat[,start_col:end_col]
        rownames(mat) <- x@cells 
        expected <- x@expected[start_col: end_col]
        list(
          matrix = mat,
          expected = expected, 
          weight = nrow(mat) # number of cells
          )
      })
      
      fp <- do.call(rbind, lapply(all.mat, '[[', "matrix"))
      expected.mat <- do.call(rbind, lapply(all.mat, function(x) x$expected * x$weight))
      weights <- vapply(all.mat, '[[', numeric(1), "weight")
      # expected.weighted.mean
      expected <- colSums(expected.mat)/sum(weights) 
      
    } else { # if only one agg.obj
      agg <- region.agg.list[[1]]
      motif.width <- width(x = agg@regions)[[1]]
      target.upstream <- agg@upstream
      target.downstream <- agg@downstream
      # extract cell x position insertion matrix 
      fp <- agg@matrix
      rownames(fp) <- agg@cells
      # extract vector of expected 
      expected <- agg@expected
    }

    ## --- Split into groups ----------------
    # background normalization by group 
    bg.norm <- lapply(X = all.groups, FUN = function(x) {
      cells.use <- names(x = obj.groups)[obj.groups == x]
      cells.present <- intersect(cells.use, rownames(fp))
      if (length(cells.present) == 0) {
        return(NULL)
      }
      mat.use <- fp[cells.use, , drop = FALSE]
      return(BackgroundMeanNorm(x = mat.use, background = 50))
    })
    bg.norm <- Filter(Negate(is.null), bg.norm)
    bg.norm <- do.call(what = rbind, args = bg.norm)
    # add position
    center.offset <- floor(motif.width/2)
    positions <-seq(
      from = -target.upstream - center.offset,
      to = target.downstream + motif.width - center.offset - 1
    )
    
    colnames(bg.norm) <- positions
    groupmeans <- ApplyMatrixByGroup(
      mat = bg.norm,
      groups = obj.groups[rownames(bg.norm)],
      fun = colMeans,
      normalize = FALSE
    )
    # add feature information 
    groupmeans$feature <- feature 
    groupmeans$class <- "Observed"
    groupmeans$position <- positions
    # add expected insertions
    expect.df <- data.frame(
      group = NA,
      count = expected,
      norm.value = expected,
      position = positions,
      feature = feature,
      class = "Expected"
    )
    groupmeans <- rbind(groupmeans, expect.df)
    return(groupmeans)
  })
  
  plot.data <- Filter(Negate(is.null), plot.data)
  plot.data <- do.call(what = rbind, args = plot.data)
  if (!is.null(x = levels.stash)) {
      plot.data$group <- factor(x = plot.data$group, levels = levels.stash)
  }
  return(plot.data)
}

#' @param regions A set of genomic ranges containing the motif instances. These
#' should all be the same width.
#' @param genome A `BSgenome` object or any other object supported by
#' `getSeq`. Do `showMethods("getSeq")` to get the list of all
#' supported object types.
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
#' @importFrom Seqinfo seqlengths
#' @export
#' @concept footprinting
#' @rdname Footprint
#' @method Footprint ChromatinAssay5
Footprint.ChromatinAssay5 <- function(
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
    overwrite = FALSE,
    verbose = TRUE,
    ...
) {
  if (is.null(x = motif.name) && is.null(x = regions)) {
    stop("Must supply the name of a motif or a set of regions")
  } else if (!is.null(x = motif.name) && !is.null(x = regions)) {
    stop("Supplied both a motif name and set of regions. Choose one only.")
  } else if (!is.null(x = motif.name)) {
    if (!inherits(x = object, what = "GRangesAssay")) {
      stop("Must supply motif positions")
    } else {
      # pull motif positions from object
      motif.obj <- Motifs(object = object)
    }
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
    if (!inherits(x = regions, what = "list")) {
      regionlist <- list(regions)
    } else {
      regionlist <- regions
    }
    all.widths <- sapply(X = regionlist, FUN = function(x) {
      length(x = unique(x = width(x = x))) == 1
    })
    if (!all(all.widths)) {
      stop("Manually-supplied regions must all have the same width.")
    }
  }
  if (compute.expected) {
    # check that bias is computed
    bias <- GetAssayData(object = object, layer = "bias")
    if (is.null(x = bias)) {
      if (verbose) {
        message("Computing Tn5 insertion bias")
      }
      region.end <- seqlengths(x = genome)[1]
      object <- InsertionBias(
        object = object,
        genome = genome,
        region = paste0(
          names(x = region.end),
          ":1-",
          as.character(x = region.end)
        )
      )
    }
  }
  # if overwrite is FALSE, skip motifs that already exist for Cells(object)
  if (!overwrite){
    existing <- RegionAggr(object)
    if (length(existing) > 0) {
      existing.names <- vapply(existing, function(x) x@name, character(1))
      
      keep.idx <- logical(length(key)) # which motifs to recompute
      for (i in seq_along(key)){
        same.name.idx <- which(existing.names ==  key[[i]])
        if (length(same.name.idx) == 0){
          keep.idx[i] <- TRUE 
          next
        }
        # union cells 
        existing.cells <- unique(unlist(
          lapply(existing[same.name.idx], function(x) x@cells),
          use.names = FALSE
        ))
        if (all(Cells(object) %in% existing.cells)){
          warning(sprintf(paste0(
            "Footprint for '%s' already exists and will not be recomputed. \n", 
            "Set overwrite=TRUE to replace the existing Footprint, ", 
            "or supply a different name to store it separately"
          ), key[[i]]), call. = FALSE)
          
          # subset key and regionlist
          keep.idx[i] <- FALSE
        } else {
          keep.idx[i] <- TRUE
        }
      }
      key <- key[keep.idx]
      regionlist <- regionlist[keep.idx]
    }
  }
  
  
  # run in parallel
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  
  # v2 edit 
  reg.agg.list <- mylapply(
    X = seq_along(along.with = regionlist),
    FUN = function(x) {
      RunFootprint(
        object = object,
        genome = genome,
        regions = regionlist[[x]],
        name = key[[x]],
        upstream = upstream,
        downstream = downstream,
        compute.expected = compute.expected,
        in.peaks = in.peaks,
        verbose = verbose
      )
    }
  )
  # store in object
  object <- SetAssayData(
    object = object, 
    layer = "region.aggregation", 
    new.data = reg.agg.list,
    overwrite = overwrite
  )

  return(object)
}

#' @rdname Footprint
#' @param assay Name of assay to use
#' @method Footprint Seurat
#' @export
#' @concept footprinting
#' @importFrom SeuratObject DefaultAssay
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
  assay <- assay %||% DefaultAssay(object = object)
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

#' @param genome A `BSgenome` object or any other object supported by
#' `getSeq`. Do `showMethods("getSeq")` to get the list of all
#' supported object types.
#' @param region Genomic region to use when assessing bias.
#' @param verbose Display messages
#'
#' @importFrom IRanges IRanges
#' @export
#' @concept footprinting
#' @rdname InsertionBias
#' @method InsertionBias ChromatinAssay5
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' region.use <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr2"),
#'   IRanges(start = c(1, 1), end = c(195471971, 182113224))
#' )
#'
#' InsertionBias(
#'   object = atac_small,
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   region = region.use
#' )
#' }
InsertionBias.ChromatinAssay5 <- function(
  object,
  genome,
  region = "chr1:1-249250621",
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- GRanges(region)
  }
  chr.use <- as.character(x = seqnames(x = region))
  reads <- MultiGetReadsInRegion(object = object, region = region)
  insertions <- GRanges(
    seqnames = c(reads$chr, reads$chr),
    ranges = IRanges(
      start = c(reads$start, reads$end),
      width = 1
    ),
    strand = "+"
  )
  insertions <- Extend(x = insertions, upstream = 3, downstream = 2)
  sequences <- as.vector(x = Biostrings::getSeq(x = genome, insertions))
  seq.freq <- table(sequences)
  # remove sequences containing non-ATCG characters
  keep.seq <- !grepl(pattern = "[^ATCGatcg]", x = names(x = seq.freq))
  insertion_hex_freq <- as.matrix(x = seq.freq[keep.seq])
  genome_freq <- Biostrings::oligonucleotideFrequency(
    x = Biostrings::getSeq(x = genome, chr.use),
    width = 6
  )
  if (inherits(x = genome_freq, what = "matrix")) {
    genome_freq <- genome_freq[1, ]
  }
  if (nrow(x = insertion_hex_freq) != length(x = genome_freq)) {
    stop("Not all hexamers represented in input region")
  }
  insertion_hex_freq <- insertion_hex_freq[names(x = genome_freq), ]
  bias <- insertion_hex_freq / genome_freq
  object <- SetAssayData(object = object, layer = "bias", new.data = bias)
  return(object)
}

#' @param assay Name of assay to use
#' @rdname InsertionBias
#' @method InsertionBias Seurat
#' @importFrom SeuratObject DefaultAssay
#' @export
#' @concept footprinting
InsertionBias.Seurat <- function(
  object,
  genome,
  assay = NULL,
  region = "chr1:1-249250621",
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
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
  flanks <- mean(x = x[, positions.use])
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

    # remove frequencies not present in hex.key
    frequencies <- frequencies[names(x = frequencies) %in% names(x = hex.key)]

    i[current.pos:end.pos] <- as.vector(x = hex.key[names(x = frequencies)])
    # shift current position
    current.pos <- end.pos + 1
  }
  # trim vectors
  x <- x[1:(current.pos - 1)]
  i <- i[1:(current.pos - 1)]
  j <- j[1:(current.pos - 1)]

  # construct matrix
  hexamer.matrix <- sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(x = hex.key), total.hexamer.positions)
  )
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
    x = c(
      expected.insertions[1:50],
      expected.insertions[
        (total.hexamer.positions - 50):total.hexamer.positions
      ]
    )
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
#' @importFrom SeuratObject DefaultAssay GetAssayData
GetMotifSize <- function(
    object,
    features,
    assay = NULL
) {
    assay <- assay %||% DefaultAssay(object = object)
    agg.list <- RegionAggr(object[[assay]])
    
    sizes <- c()
    agg.list.names <- vapply(agg.list, FUN = function(x) x@name, FUN.VALUE = character(1))
    for (i in features) {
        # motif <- positionEnrichment[[i]]["motif", ]
        if (!i %in% agg.list.names){
            stop("No footprinting data found for feature: ", i)
        }
        motif.agg <- agg.list[which(agg.list.names==i)][[1]]
        w <- unique(width(motif.agg@regions))
        if (length(w) != 1){
            stop("Regions for feature ", i, " have inconsistent widths")
        }
        sizes <- c(sizes, w)
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
#' @importFrom BiocGenerics sort
RunFootprint <- function(
    object,
    genome,
    regions,   
    name = NULL , 
    upstream = 250,
    downstream = 250,
    compute.expected = TRUE,
    in.peaks = FALSE,
    verbose = TRUE
) {
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
        stop("Please install Biostrings: BiocManager::install('Biostrings')")
    }
    # resolve footprint name 
    if (is.null(name)){
        # try to infer from GRanges names
        gr.names <- unique(names(regions))
        gr.names <- gr.names[!is.na(gr.names)]
        if (length(gr.names) ==1){
          name <- gr.names
        } else {
          stop("Footprint name not provided and could not be inferred from regions (",
               gr.names, "), please supply name explicitly")
        }
    }
    motif.size <- width(x = regions)[[1]]
    regions <- sort(x = regions)
    if (in.peaks) {
        regions <- subsetByOverlaps(x = regions, ranges = granges(x = object))
    }
    motif.regions <- regions
    # extend upstream and downstream
    regions <- Extend(
        x = regions,
        upstream = upstream,
        downstream = downstream
    )
    if (verbose) {
        message("Computing observed Tn5 insertions per base")
    }
    dna.sequence <- Biostrings::getSeq(x = genome, Extend(
        x = regions,
        upstream = 3,
        downstream = 3
    ))
    if (compute.expected) {
        bias <- GetAssayData(object = object, layer = "bias")
        if (is.null(x = bias)) {
            stop("Insertion bias not computed")
        } else {
            # add three bases each side here so we can get the hexamer frequencies
            # for every position
            expected.insertions <- FindExpectedInsertions(
                dna.sequence = dna.sequence,
                bias = bias,
                verbose = verbose
            )
        }
    } else {
        expected.insertions <- rep(1, width(x = dna.sequence)[[1]] - 6)
    }
    
    # count insertions at each position for each cell
    insertion.matrix <- CreateRegionPileupMatrix(
        object = object,
        regions = regions 
    ) # returns a sparse dgcMatrix that does not have nrows()
    # get expected insertions 
    expected.insertions <- as.numeric(x = expected.insertions)
    agg.obj <- CreateRegionAggregationObject(
        mat = insertion.matrix, 
        regions = motif.regions, 
        upstream = upstream, 
        downstream = downstream, 
        name = name, 
        expected = expected.insertions, 
        cells = Cells(object)
    )

    return(agg.obj)
} 

#' @include generics.R
#'
NULL

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' AverageCounts
#'
#' Compute the mean counts per group of cells for a given assay
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. Default is the active assay
#' @param group.by Grouping variable to use. Default is the active identities
#' @param verbose Display messages
#' @importFrom Seurat DefaultAssay Idents GetAssayData
#' @importFrom Matrix colSums
#' @importFrom dplyr group_by summarize
#' @export
AverageCounts <- function(
  object,
  assay = NULL,
  group.by = NULL,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = group.by)) {
    group.by <- Idents(object = object)
  } else {
    group.by <- object[[group.by, drop = TRUE]]
  }
  counts <- GetAssayData(object = object, assay = assay, slot = 'counts')
  if (verbose) {
    message('Summing counts per cell')
  }
  totals <- colSums(x = counts)
  total.df <- data.frame(cell = names(x = totals), readcount = totals, stringsAsFactors = FALSE)
  total.df$group <- group.by[total.df$cell]
  total.df <- group_by(total.df, group)
  if (verbose) {
    message('Computing average counts per group')
  }
  group.means <- summarize(.data = total.df, mn = mean(x = readcount))
  results <- group.means$mn
  names(x = results) <- group.means$group
  return(results)
}

#' CellsPerGroup
#'
#' Count the number of cells in each group
#'
#' @param object A Seurat object
#' @param group.by A grouping variable. Default is the active identities
#' @return Returns a vector
#' @importFrom Seurat Idents
#' @export
CellsPerGroup <- function(
  object,
  group.by = NULL
) {
  if (is.null(x = group.by)) {
    cellgroups <- Idents(object = object)
  } else {
    meta.data <- object[[]]
    cellgroups <- meta.data[[group.by]]
  }
  cells.per.group <- table(cellgroups)
  lut <- as.vector(x = cells.per.group)
  names(x = lut) <- names(x = cells.per.group)
  return(lut)
}

#' ClosestFeature
#'
#' Find the closest feature to a given set of genomic regions
#'
#' @param regions A set of genomic regions to query
#' @param annotation Annotation information. Can be a GRanges object or an EnsDb object
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors subjectHits mcols
#' @importFrom GenomicFeatures genes
#' @importFrom GenomeInfoDb seqlevelsStyle "seqlevelsStyle<-"
#' @importFrom methods is
#'
#' @return Returns a dataframe with the name of each region, the closest feature in the annotation,
#' and the distance to the feature.
#'
#' @export
ClosestFeature <- function(
  regions,
  annotation,
  ...
) {
  if (!is(object = regions, class2 = 'GRanges')) {
    regions <- StringToGRanges(regions = regions, ...)
  }
  if (is(object = annotation, class2 = 'EnsDb')) {
    annotation <- genes(x = annotation, filter = ~ gene_biotype == "protein_coding")
    if (seqlevelsStyle(x = regions) != seqlevelsStyle(x = annotation)) {
      seqlevelsStyle(x = annotation) <- seqlevelsStyle(x = regions)
    }
  }
  nearest_feature <- distanceToNearest(x = regions, subject = annotation)
  feature_hits <- annotation[subjectHits(x = nearest_feature)]
  df <- as.data.frame(x = mcols(x = feature_hits))
  df$region <- GRangesToString(grange = feature_hits)
  df$distance <- mcols(x = nearest_feature)$distance
  return(df)
}

#' Set the fragments file path for creating plots
#'
#' Give path of indexed fragments file that goes with data in the object.
#' Checks for a valid path and an index file with the same name (.tbi) at the same path.
#' Stores the path under the tools slot for access by visualization functions.
#' One fragments file can be stored for each assay.
#'
#' @param object A Seurat object
#' @param file Path to indexed fragment file.
#' See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @param assay Assay used to generate the fragments. If NULL, use the active assay.
#'
#' @importFrom methods "slot<-" slot is
#'
#' @export
#'
SetFragments <- function(
  object,
  file,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!(assay %in% names(x = slot(object = object, name = 'assays')))) {
    stop('Requested assay not present in object')
  }
  index.file <- paste0(file, '.tbi')
  if (all(file.exists(file, index.file))) {
    file <- normalizePath(path = file)
    current.tools <- slot(object = object, name = 'tools')
    current.tools$fragments[[assay]] <- file
    slot(object = object, name = 'tools') <- current.tools
    return(object)
  } else {
    stop('Requested file does not exist or is not indexed')
  }
}

#' StringToGRanges
#'
#' Convert a genomic coordinate string to a GRanges object
#'
#' @param regions Vector of genomic region strings
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome
#' and coordinates, second separator is used to separate start and end coordinates.
#' @return Returns a GRanges object
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate
#' @examples
#' regions <- c('chr1-1-10', 'chr2-12-3121')
#' StringToGRanges(regions = regions)
#'
#' @export
StringToGRanges <- function(regions, sep = c("-", "-")) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = 'ranges',
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c('chr', 'start', 'end')
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df)
  return(granges)
}

#' GRangesToString
#'
#' Convert GRanges object to a vector of strings
#'
#' @param grange A GRanges object
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome
#' and coordinates, second separator is used to separate start and end coordinates.
#' @importFrom GenomicRanges seqnames start end
#' @examples
#' GRangesToString(grange = blacklist_hg19)
#'
#' @export
GRangesToString <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange)
  )
  return(regions)
}

#' CalculateCoverages
#'
#' Calculate normalized read coverage per base per cell group
#'
#' @param reads Dataframe containing reads for a region, generated by \code{\link{GetReadsInRegion}}
#' @param cells.per.group Number of cells in each group. Must be a named vector. Used for normalization of track heights.
#' @param reads.per.group Average number of reads per group.
#' @param scale.factor Scaling factor for track height. Default is 10^4.
#' @param window Smoothing window to use
#' @param verbose Display messages
#'
#' @importFrom zoo rollapply
#' @importFrom dplyr group_by mutate summarize arrange ungroup
#' @export
CalculateCoverages <- function(
  reads,
  cells.per.group,
  reads.per.group,
  scale.factor = 10^4,
  window = 100,
  verbose = TRUE
) {
  if (verbose) {
    message("Computing coverage per base")
  }
  templist <- list()
  y <- 1
  n_allocate <- sum(reads$length)
  templist$position <- rep(x = 0, n_allocate)
  templist$value <- rep(x = 0, n_allocate)
  templist$cell <- rep(x = '0', n_allocate)
  templist$group <- rep(x = '0', n_allocate)
  startpos <- reads[,][['start']]
  endpos <- reads[,][['end']]
  cellnames <- as.character(x = reads[,][['cell']])
  groupnames <- as.character(x = reads[,][['group']])
  for (i in seq_len(length.out = nrow(x = reads))) {
    interval <- startpos[[i]]:endpos[[i]]
    for (j in interval) {
      templist$position[y] <- j
      templist$value[y] <- 1
      templist$cell[[y]] <- cellnames[[i]]
      templist$group[[y]] <- groupnames[[i]]
      y <- y + 1
    }
  }
  expanded <- as.data.frame(x = t(x = do.call(what = rbind, args = templist)), stringsAsFactors = FALSE)
  expanded$position <- as.numeric(x = expanded$position)
  expanded$value <- as.numeric(x = expanded$value)
  if (is(object = reads$group, class2 = 'factor')) {
    expanded$group <- factor(x = expanded$group, levels = levels(reads$group))
  }
  expanded$norm.value <- expanded$value / reads.per.group[expanded$group] / cells.per.group[expanded$group] * scale.factor
  expanded <- group_by(.data = expanded, position, group)
  coverages <- summarize(.data = expanded, total = sum(norm.value))
  coverages <- group_by(.data = coverages, group)
  coverages <- arrange(.data = coverages, position)
  coverages <- mutate(.data = coverages, coverage = rollapply(
    data = total,
    width = window,
    FUN = mean,
    align = 'center',
    fill = NA
  ))
  coverages <- ungroup(x = coverages)
  return(coverages)
}

#' ChunkGRanges
#'
#' Split a genomic ranges object into evenly sized chunks
#'
#' @param granges A GRanges object
#' @param nchunk Number of chunks to split into
#'
#' @return Returns a list of GRanges objects
#' @export
ChunkGRanges <- function(granges, nchunk) {
  chunksize <- as.integer(x = (length(granges) / nchunk))
  range.list <- sapply(X = 1:nchunk, FUN = function(x) {
    chunkupper <- (x * chunksize)
    if (x == 1) {
      chunklower <- 1
    } else {
      chunklower <- ((x-1) * chunksize) + 1
    }
    if (x == nchunk) {
      chunkupper <- length(x = granges)
    }
    return(granges[chunklower:chunkupper])
  })
  return(range.list)
}

#' Extend
#'
#' Resize GenomicRanges upstream and or downstream.
#' From \url{https://support.bioconductor.org/p/78652/}
#'
#' @param x A range
#' @param upstream Length to extend upstream
#' @param downstream Length to extend downstream
#'
#' @importFrom GenomicRanges strand start end trim
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
Extend <- function(x, upstream = 0, downstream = 0) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  new_start <- start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
  new_end <- end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}

#' GetCellsInRegion
#'
#' Extract cell names containing reads mapped within a given genomic region
#'
#' @param tabix Tabix object
#' @param region A string giving the region to extract from the fragments file
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome
#' and coordinates, second separator is used to separate start and end coordinates.
#' @param cells Vector of cells to include in output. If NULL, include all cells
#'
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom methods is
#' @export
GetCellsInRegion <- function(tabix, region, sep = c("-", "-"), cells = NULL) {
  if (!is(object = region, class2 = 'GRanges')) {
    region <- StringToGRanges(regions = region)
  }
  bin.reads <- scanTabix(file = tabix, param = region)
  reads <- sapply(X = bin.reads, FUN = ExtractCell)
  if (!is.null(x = cells)) {
    reads <- sapply(X = reads, FUN = function(x) {
      x <- x[x %in% cells]
      if (length(x = x) == 0) {
        return(NULL)
      } else {
        return(x)
      }
    })
  }
  nrep <- sapply(X = reads, FUN = length)
  regions <- rep(x = names(x = reads), nrep)
  regions <- gsub(pattern = ":", replacement = sep[[1]], x = regions)
  regions <- gsub(pattern = "-", replacement = sep[[2]], x = regions)
  cellnames <- as.vector(x = unlist(x = reads))
  return(list(cells = cellnames, region = regions))
}

#' GetReadsInRegion
#'
#' Extract reads for each cell within a given genomic region or set of regions
#'
#' @param object A Seurat object
#' @param region A genomic region, specified as a string in the format
#' 'chr:start-end'. Can be a vector of regions.
#' @param assay Name of assay to use
#' @param fragment.path Path to indexed fragment file
#' @param group.by Cell grouping information to add
#' @param cells Cells to include. Default is all cells present in the object.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom Seurat Idents DefaultAssay
#'
#' @return Returns a data frame
#' @export
GetReadsInRegion <- function(
  object,
  region,
  assay = NULL,
  fragment.path = NULL,
  group.by = NULL,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = group.by)) {
    group.by <- Idents(object = object)
  } else {
    meta.data <- object[[]]
    group.by <- meta.data[[group.by]]
    names(x = group.by) <- rownames(x = meta.data)
  }
  fragment.path <- fragment.path %||% GetFragments(object = object, assay = assay)
  if (verbose) {
    message('Extracting reads in requested region')
  }
  if (!is(object = region, class2 = 'GRanges')) {
    region <- StringToGRanges(regions = region, ...)
  }
  tbx <- TabixFile(file = fragment.path)
  reads <- scanTabix(file = tbx, param = region)
  reads <- TabixOutputToDataFrame(reads = reads)
  reads <- reads[reads$cell %in% names(group.by), ]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    stop('No cells present in the requested region')
  }
  reads$length <- reads$end - reads$start
  reads$group <- group.by[reads$cell]
  return(reads)
}

#' GetFragments
#'
#' Retrieve path to fragments file from assay object, and checks that the file exists and
#' is indexed before returning the file path.
#'
#' @param object A Seurat object
#' @param assay Name of the assay use to store the fragments file path
#'
#' @importFrom methods slot
#'
#' @return Returns the path to a fragments file stored in the Assay if present
#' @export
GetFragments <- function(
  object,
  assay
) {
    tools <- slot(object = object, name = 'tools')
    if ('fragments' %in% names(x = tools)) {
      if (assay %in% names(x = tools$fragments)) {
        fragment.path <- tools$fragments[[assay]]
      } else {
        stop('Fragment file not supplied for the requested assay')
      }
    } else {
      stop('Fragment file not set. Run SetFragments to set the fragment file path.')
    }
  if (!(all(file.exists(fragment.path, paste0(fragment.path, '.tbi'))))) {
    stop('Requested file does not exist or is not indexed')
  } else {
    return(fragment.path)
  }
}

#' CountsInRegion
#'
#' Count reads per cell overlapping a given set of regions
#'
#' @param object A Seurat object
#' @param assay Name of assay in the object to use
#' @param regions A GRanges object
#' @param sep Separator to use when extracting genomic coordinates from the Seurat object
#' @param ... Additional arguments passed to \code{\link[IRanges]{findOverlaps}}
#'
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom Matrix colSums
#'
#' @export
CountsInRegion <- function(
  object,
  assay,
  regions,
  sep = c("-", '-'),
  ...
) {
  obj.regions <- rownames(x = object[[assay]])
  obj.granges <- StringToGRanges(regions = obj.regions, sep = sep)
  overlaps <- findOverlaps(query = obj.granges, subject = regions, ...)
  hit.regions <- GRangesToString(grange = obj.granges[queryHits(x = overlaps)], sep = sep)
  data.matrix <- GetAssayData(object = object, assay = assay, slot = 'counts')[hit.regions, ]
  return(colSums(data.matrix))
}

#' ExtractCell
#'
#' Extract cell barcode from list of tab delimited character
#' vectors (output of \code{\link{scanTabix}})
#'
#' @param x List of character vectors
#' @export
ExtractCell <- function(x) {
  if (length(x = x) == 0) {
    return(NULL)
  } else {
    tmp <- strsplit(x = x, split = "\t")
    return(unlist(x = tmp)[5*(1:length(x = tmp))-1])
  }
}

#' FractionCountsInRegion
#'
#' Find the fraction of counts per cell that overlap a given set of genomic ranges
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param regions A GRanges object containing a set of genomic regions
#' @param sep The separator used to separate genomic coordinate information in the assay feature names
#' @param ... Additional arguments passed to \code{\link{CountsInRegion}}
#'
#' @export
FractionCountsInRegion <- function(
  object,
  assay,
  regions,
  sep = c("-", "-"),
  ...
) {
  reads.in.region <- CountsInRegion(
    object = object,
    regions = regions,
    assay = assay,
    sep = sep,
    ...
  )
  total.reads <- colSums(x = GetAssayData(object = object, assay = assay, slot = 'counts'))
  return(reads.in.region / total.reads)
}

#' TabixOutputToDataFrame
#'
#' Create a single dataframe from list of character vectors
#'
#' @param reads List of character vectors (the output of \code{\link{scanTabix}})
#' @param record.ident Add a column recording which region the reads overlapped with (default TRUE)
#' @importFrom data.table rbindlist
#' @importFrom utils read.table
#' @return Returns a data.frame
#' @export
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  df.list <- lapply(X = 1:length(reads), FUN = function(x) {
    if (length(x = reads[[x]]) == 0) {
      return(NULL)
    }
    df <- read.table(
      file = textConnection(object = reads[[x]]),
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    colnames(x = df) <- c('chr', 'start', 'end', 'cell', 'count')
    if (record.ident) {
      df$ident <- x
    }
    return(df)
  })
  return(rbindlist(l = df.list))
}

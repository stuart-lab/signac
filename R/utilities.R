#' @include generics.R
#' @importFrom utils globalVariables
NULL

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param y The value to provide if x is null
# @return Returns y if x is null, otherwise returns x.
SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}

#' @importFrom Seurat DefaultAssay Misc
AddToMisc <- function(
  object,
  new.data,
  save.as,
  assay = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  misc.slot <- SetIfNull(x = Misc(object = object[[assay]]), y = list())
  if (!inherits(x = misc.slot, what = 'list')) {
    warning("Misc slot already occupied")
  } else{
    misc.slot[[save.as]] <- new.data
    object[[assay]]@misc <- misc.slot
  }
  return(object)
}

globalVariables(names = c('group', 'readcount'), package = 'Signac')
#' Average Counts
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
#' @return Returns a dataframe
#' @examples
#' AverageCounts(atac_small)
AverageCounts <- function(
  object,
  assay = NULL,
  group.by = NULL,
  verbose = TRUE
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
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

#' Cells per group
#'
#' Count the number of cells in each group
#'
#' @param object A Seurat object
#' @param group.by A grouping variable. Default is the active identities
#' @importFrom Seurat Idents
#' @export
#' @return Returns a vector
#' @examples
#' CellsPerGroup(atac_small)
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

#' Closest Feature
#'
#' Find the closest feature to a given set of genomic regions
#'
#' @param regions A set of genomic regions to query
#' @param annotation Annotation information. Can be a GRanges object or an EnsDb object.
#' If an EnsDb object is provided, protein-coding genes will be extracted from the
#' object and only the closest protein coding genes are reported. If a GRanges
#' object is provided, no filtering is performed and the closest genomic range
#' is reported.
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors subjectHits mcols
#' @importFrom GenomicFeatures genes
#' @importFrom GenomeInfoDb seqlevelsStyle "seqlevelsStyle<-"
#' @importFrom methods is
#' @return Returns a dataframe with the name of each region, the closest feature in the annotation,
#' and the distance to the feature.
#' @export
#' @examples
#' \donttest{
#' ClosestFeature(
#'   regions = head(rownames(atac_small)),
#'   annotation = StringToGRanges(head(rownames(atac_small)), sep = c(':', '-')),
#'   sep = c(":", "-")
#' )
#' }
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
  df$closest_region <- GRangesToString(grange = feature_hits, ...)
  df$query_region <- GRangesToString(grange = regions, ...)
  df$distance <- mcols(x = nearest_feature)$distance
  return(df)
}

#' Find interesecting regions between two objects
#'
#' Intersects the regions stored in the rownames of two objects and
#' returns a vector containing the names of rows that interesect
#' for each object. The order of the row names return corresponds
#' to the intersecting regions, ie the nth feature of the first vector
#' will intersect the nth feature in the second vector. A distance
#' parameter can be given, in which case features within the given
#' distance will be called as intersecting.
#'
#' @param object.1 The first Seurat object
#' @param object.2 The second Seurat object
#' @param assay.1 Name of the assay to use in the first object. If NULL, use
#' the default assay
#' @param assay.2 Name of the assay to use in the second object. If NULL, use
#' the default assay
#' @param sep.1 Genomic coordinate separators to use for the first object
#' @param sep.2 Genomic coordinate separators to use for the second object
#' @param distance Maximum distance between regions allowed for an intersection to
#' be recorded. Default is 0.
#' @param verbose Display messages
#'
#' @importFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors subjectHits queryHits mcols
#' @importFrom Seurat DefaultAssay
#' @export
#' @return Returns a list of two character vectors containing the row names
#' in each object that overlap each other.
#' @examples
#' GetIntersectingFeatures(
#'   object.1 = atac_small,
#'   object.2 = atac_small,
#'   assay.1 = 'peaks',
#'   assay.2 = 'bins',
#'   sep.1 = c(":", "-"),
#'   sep.2 = c("-", "-")
#' )
GetIntersectingFeatures <- function(
  object.1,
  object.2,
  assay.1 = NULL,
  assay.2 = NULL,
  distance = 0,
  sep.1 = c("-", "-"),
  sep.2 = c("-", "-"),
  verbose = TRUE
) {
  assay.1 <- SetIfNull(x = assay.1, y = DefaultAssay(object = object.1))
  assay.2 <- SetIfNull(x = assay.2, y = DefaultAssay(object = object.2))
  regions.1 <- StringToGRanges(regions = rownames(x = object.1[[assay.1]]), sep = sep.1)
  regions.2 <- StringToGRanges(regions = rownames(x = object.2[[assay.2]]), sep = sep.2)
  if (verbose) {
    message("Intersecting regions across objects")
  }
  region.intersections <- distanceToNearest(x = regions.1, subject = regions.2)
  keep.intersections <- mcols(x = region.intersections)$distance <= distance
  region.intersections <- region.intersections[keep.intersections, ]
  intersect.object1 <- regions.1[queryHits(x = region.intersections)]
  intersect.object2 <- regions.2[subjectHits(x = region.intersections)]
  regions.obj1 <- GRangesToString(grange = intersect.object1, sep = sep.1)
  regions.obj2 <- GRangesToString(grange = intersect.object2, sep = sep.2)
  return(list(regions.obj1, regions.obj2))
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
#' @export
#' @return Returns a Seurat object
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' SetFragments(object = atac_small, file = fpath)
SetFragments <- function(
  object,
  file,
  assay = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
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

#' String to GRanges
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

#' GRanges to String
#'
#' Convert GRanges object to a vector of strings
#'
#' @param grange A GRanges object
#' @param sep Vector of separators to use for genomic string. First element is used to separate chromosome
#' and coordinates, second separator is used to separate start and end coordinates.
#' @importFrom GenomicRanges seqnames start end
#' @examples
#' GRangesToString(grange = blacklist_hg19)
#' @return Returns a character vector
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

# Chunk GRanges
#
# Split a genomic ranges object into evenly sized chunks
#
# @param granges A GRanges object
# @param nchunk Number of chunks to split into
#
# @return Returns a list of GRanges objects
# @examples
# ChunkGRanges(blacklist_hg19, n = 10)
ChunkGRanges <- function(granges, nchunk) {
  chunksize <- as.integer(x = (length(granges) / nchunk))
  range.list <- sapply(X = seq_len(length.out = nchunk), FUN = function(x) {
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

#' Generate matrix of integration sites
#'
#' Generates a cell-by-position matrix of Tn5 integration sites
#' centered on a given region (usually a DNA sequence motif). This
#' matrix can be used for downstream footprinting analysis.
#'
#' @param object A Seurat object
#' @param region A GRanges object containing the region of interest
#' @param assay Name of the assay to use
#' @param cells Which cells to include in the matrix. If NULL (default), use all
#' cells in the object
#' @param tabix.file A TabixFile object. If NULL, the file specified in \code{fragment.path}
#' will be opened and closed after the function completes. If iterating over many regions, providing an
#' open TabixFile is much faster as it avoids opening and closing the connection each time.
#' @param verbose Display messages
#' @importFrom BiocGenerics width start end
#' @return Returns a sparse matrix
#' @export
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' atac_small <- SetFragments(atac_small, file = fpath)
#' CutMatrix(
#'  object = atac_small,
#'  region = StringToGRanges("chr1-10245-762629")
#' )
CutMatrix <- function(
  object,
  region,
  tabix.file = NULL,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  if (!inherits(x = region, what = 'GRanges')) {
    stop("Region is not a GRanges object.")
  }
  all.cells <- SetIfNull(x = cells, y = colnames(x = object))
  if (is.null(x = tabix.file)) {
    fragment.path <- GetFragments(object = object, assay = assay)
  }
  fragments <- GetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    tabix.file = tabix.file,
    verbose = verbose
  )
  # if there are no reads in the region, create an empty matrix of the correct dimension
  if (nrow(x = fragments) == 0) {
    cut.matrix <- sparseMatrix(
      i = NULL,
      j = NULL,
      dims = c(length(x = all.cells), width(x = region))
    )
  } else {
    cut.df <- data.frame(
      position = c(fragments$start, fragments$end) - start(x = region) + 1,
      cell = c(fragments$cell, fragments$cell),
      stringsAsFactors = FALSE
    )
    cut.df <- cut.df[cut.df$position > 0 & cut.df$position <= width(x = region), ]
    cell.vector <- seq_along(along.with = all.cells)
    names(x = cell.vector) <- all.cells
    cell.matrix.info <- cell.vector[cut.df$cell]
    cut.matrix <- sparseMatrix(
      i = cell.matrix.info,
      j = cut.df$position,
      x = 1,
      dims = c(length(x = all.cells), width(x = region))
    )
  }
  rownames(x = cut.matrix) <- all.cells
  colnames(x = cut.matrix) <- start(region):end(region)
  return(cut.matrix)
}

#' Extend
#'
#' Resize GenomicRanges upstream and or downstream.
#' From \url{https://support.bioconductor.org/p/78652/}
#'
#' @param x A range
#' @param upstream Length to extend upstream
#' @param downstream Length to extend downstream
#' @param from.midpoint Count bases from region midpoint,
#' rather than the 5' or 3' end for upstream and downstream
#' respectively.
#'
#' @importFrom GenomicRanges strand start end trim
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' @examples
#' Extend(x = blacklist_hg19, upstream = 100, downstream = 100)
Extend <- function(
  x,
  upstream = 0,
  downstream = 0,
  from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x)/2)
    new_start <- midpoints - ifelse(test = on_plus, yes = upstream, no = downstream)
    new_end <- midpoints + ifelse(test = on_plus, yes = downstream, no = upstream)
  } else {
    new_start <- start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
    new_end <- end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  }
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
#' @return Returns a list
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' GetCellsInRegion(tabix = fpath, region = "chr1-10245-762629")
GetCellsInRegion <- function(tabix, region, sep = c("-", "-"), cells = NULL) {
  if (!is(object = region, class2 = 'GRanges')) {
    region <- StringToGRanges(regions = region)
  }
  bin.reads <- scanTabix(file = tabix, param = region)
  reads <- sapply(X = bin.reads, FUN = ExtractCell, simplify = FALSE)
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
#' @param tabix.file A TabixFile object. If NULL, the file specified in \code{fragment.path}
#' will be opened and closed after the function completes. If iterating over many regions, providing an
#' open TabixFile is much faster as it avoids opening and closing the connection each time.
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
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' atac_small <- SetFragments(object = atac_small, file = fpath)
#' region <- StringToGRanges(regions = "chr1-10245-762629")
#' GetReadsInRegion(object = atac_small, region = region)
GetReadsInRegion <- function(
  object,
  region,
  assay = NULL,
  tabix.file = NULL,
  group.by = NULL,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (is.null(x = group.by)) {
    group.by <- Idents(object = object)
  } else {
    meta.data <- object[[]]
    group.by <- meta.data[[group.by]]
    names(x = group.by) <- rownames(x = meta.data)
  }
  if (verbose) {
    message('Extracting reads in requested region')
  }
  if (!is(object = region, class2 = 'GRanges')) {
    region <- StringToGRanges(regions = region, ...)
  }
  if (is.null(x = tabix.file)) {
    fragment.path <- GetFragments(object = object, assay = assay)
    tabix.file <- TabixFile(file = fragment.path)
    tbx <- open(con = tabix.file)
    close.file <- TRUE
  } else {
    close.file <- FALSE
  }
  reads <- scanTabix(file = tabix.file, param = region)
  if (close.file) {
    close(con = tabix.file)
  }
  reads <- TabixOutputToDataFrame(reads = reads)
  reads <- reads[reads$cell %in% names(group.by), ]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    return(reads)
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
#' @importFrom methods slot
#' @importFrom Seurat DefaultAssay
#' @return Returns the path to a fragments file stored in the Assay if present
#' @export
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' atac_small <- SetFragments(object = atac_small, file = fpath)
#' GetFragments(object = atac_small)
GetFragments <- function(
  object,
  assay = NULL
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
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
#' @return Returns a numeric vector
#' @examples
#' \donttest{
#' CountsInRegion(
#'   object = atac_small,
#'   assay = 'bins',
#'   regions = blacklist_hg19
#' )
#' }
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
#' @return Returns a string
#' @examples
#' ExtractCell(x = "chr1\t1\t10\tatcg\t1")
ExtractCell <- function(x) {
  if (length(x = x) == 0) {
    return(NULL)
  } else {
    tmp <- strsplit(x = x, split = "\t")
    return(unlist(x = tmp)[5*(seq_along(along.with = tmp))-1])
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
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData
#'
#' @export
#' @return Returns a numeric vector
#' @examples
#' FractionCountsInRegion(
#'   object = atac_small,
#'   assay = 'bins',
#'   regions = blacklist_hg19
#' )
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

# Get vector of cell names and associated identity
# @param object A Seurat object
# @param group.by Identity class to group cells by
# @param idents which identities to include
# @return Returns a named vector
#' @importFrom Seurat Idents
GetGroups <- function(
  object,
  group.by,
  idents
) {
  if (is.null(x = group.by)) {
    obj.groups <- Idents(object = object)
  } else {
    obj.md <- object[[group.by]]
    obj.groups <- obj.md[, 1]
    names(obj.groups) <- rownames(x = obj.md)
  }
  if (!is.null(idents)) {
    obj.groups <- obj.groups[obj.groups %in% idents]
  }
  return(obj.groups)
}

#' Intersect genomic coordinates with matrix rows
#'
#' Remove or retain matrix rows that intersect given genomic regions
#'
#' @param matrix A matrix with genomic regions in the rows
#' @param regions A set of genomic regions to intersect with regions in the matrix.
#' Either a vector of strings encoding the genomic coordinates, or a GRanges object.
#' @param invert Discard rows intersecting the genomic regions supplied, rather than retain. Default FALSE.
#' @param sep A length-2 character vector containing the separators to be used for
#' extracting genomic coordinates from a string. The first element will be used to separate the
#' chromosome name from coordinates, and the second element used to separate start and end
#' coordinates.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[IRanges]{findOverlaps}}
#'
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#'
#' @export
#' @return Returns a sparse matrix
#' @examples
#' counts <- matrix(data = rep(0, 12), ncol = 2)
#' rownames(counts) <- c("chr1-565107-565550","chr1-569174-569639",
#' "chr1-713460-714823","chr1-752422-753038",
#' "chr1-762106-763359","chr1-779589-780271")
#' IntersectMatrix(matrix = counts, regions = blacklist_hg19)
IntersectMatrix <- function(
  matrix,
  regions,
  invert = FALSE,
  sep = c("-", "-"),
  verbose = TRUE,
  ...
) {
  if (is(object = regions, class2 = 'character')) {
    regions <- StringToGRanges(regions = regions, sep = sep)
  }
  rowranges <- StringToGRanges(regions = rownames(x = matrix), sep = sep)
  if (verbose) {
    message("Intersecting genomic regions")
  }
  region.overlaps <- findOverlaps(query = rowranges, subject = regions, ...)
  keep.rows <- queryHits(x = region.overlaps)
  if (invert) {
    all.rows <- seq_len(length.out = nrow(matrix))
    keep.rows <- setdiff(x = all.rows, y = keep.rows)
  }
  if (verbose) {
    message("Subsetting matrix")
  }
  matrix <- matrix[keep.rows, ]
  return(matrix)
}

#' Match DNA sequence characteristics
#'
#' Return a vector if genomic regions that match the distribution of a set of query regions
#' for any given set of characteristics, specified in the input \code{meta.feature} dataframe.
#'
#' @param meta.feature A dataframe containing DNA sequence information
#' @param regions Set of query regions. Must be present in rownames.
#' @param n Number of regions to select, with characteristics matching the query
#' @param features.match Which features of the query to match when selecting a set of
#' regions. A vector of column names present in the feature metadata can be supplied to
#' match multiple characteristics at once. Default is GC content.
#' @param verbose Display messages
#' @param ... Arguments passed to other functions
#' @return Returns a character vector
#'
#' @importFrom stats density approx
#' @export
#' @examples
#' metafeatures <- Seurat::GetAssayData(object = atac_small[['peaks']], slot = 'meta.features')
#' MatchRegionStats(
#'   meta.feature = metafeatures,
#'   regions = head(rownames(metafeatures), 10),
#'   features.match = "percentile",
#'   n = 10
#' )
MatchRegionStats <- function(
  meta.feature,
  regions,
  features.match = c('GC.percent'),
  n = 10000,
  verbose = TRUE,
  ...
) {
  if (length(x = features.match) == 0) {
    stop("Must supply at least one sequence characteristic to match")
  }
  mf.query <- meta.feature[regions, ]
  choosefrom <- setdiff(x = rownames(x = meta.feature), y = rownames(x = mf.query))
  if (length(x = choosefrom) < n) {
    n <- length(x = choosefrom)
    warning("Requested more features than present in supplied data. Returning ", n, " features")
  }
  features.choose <- meta.feature[choosefrom, ]
  feature.weights <- rep(0, nrow(features.choose))
  for (i in features.match) {
    if (verbose) {
      message("Matching ", i, " distribution")
    }
    density.estimate <- density(x = mf.query[[i]], kernel = "gaussian", bw = 1)
    weights <- approx(
      x = density.estimate$x,
      y = density.estimate$y,
      xout = features.choose[[i]],
      yright = 0.0001,
      yleft = 0.0001
    )$y
    feature.weights <- feature.weights + weights
  }
  feature.select <- sample(
    x = rownames(x = features.choose),
    size = n,
    prob = feature.weights
  )
  return(feature.select)
}

#' Region-aware object merging
#'
#' This will find intersecting regions in both objects and rename the
#' overlapping features with the region coordinates of the first object
#' (by default; this can be changed with the regions.use parameter).
#'
#' This allows a merged object to be constructed with common feature names.
#'
#' @param object.1 The first Seurat object
#' @param object.2 The second Seurat object
#' @param assay.1 Name of the assay to use in the first object. If NULL, use
#' the default assay
#' @param assay.2 Name of the assay to use in the second object. If NULL, use
#' the default assay
#' @param sep.1 Genomic coordinate separators to use for the first object
#' @param sep.2 Genomic coordinate separators to use for the second object
#' @param regions.use Which regions to use when naming regions in the merged object.
#' Options are:
#' \itemize{
#'  \item{1}: Use the region coordinates from the first object
#'  \item{2}: Use the region coordinates from the second object
#' }
#' @param distance Maximum distance between regions allowed for an intersection to
#' be recorded. Default is 0.
#' @param new.assay.name Name for the merged assay. Default is 'peaks'
#' @param project Project name for the new object
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[Seurat]{CreateAssayObject}}
#'
#' @importFrom Seurat DefaultAssay CreateAssayObject GetAssayData
#' @importFrom utils packageVersion
#'
#' @export
#' @return Returns a \code{\link[Seurat]{Seurat}} object
#' @examples
#' MergeWithRegions(
#'   object.1 = atac_small,
#'   object.2 = atac_small,
#'   assay.1 = 'peaks',
#'   assay.2 = 'bins',
#'   sep.1 = c(":","-"),
#'   sep.2 = c("-","-")
#' )
MergeWithRegions <- function(
  object.1,
  object.2,
  assay.1 = NULL,
  assay.2 = NULL,
  sep.1 = c("-", "-"),
  sep.2 = c("-", "-"),
  regions.use = 1,
  distance = 0,
  new.assay.name = 'peaks',
  project = 'SeuratProject',
  verbose = TRUE,
  ...
) {
  assay.1 <- SetIfNull(x = assay.1, y = DefaultAssay(object = object.1))
  assay.2 <- SetIfNull(x = assay.2, y = DefaultAssay(object = object.2))
  intersecting.regions <- GetIntersectingFeatures(
    object.1 = object.1,
    object.2 = object.2,
    assay.1 = assay.1,
    assay.2 = assay.2,
    sep.1 = sep.1,
    sep.2 = sep.2,
    distance = distance,
    verbose = verbose
  )
  regions.obj1 <- intersecting.regions[[1]]
  regions.obj2 <- intersecting.regions[[2]]
  # TODO add option to keep non-overlapping regions
  if (regions.use == 1) {
    region.names <- regions.obj1
  } else if (regions.use == 2) {
    region.names <-regions.obj2
  } else {
    # TODO add option to rename regions as coordinate merge
    # TODO add option to rename regions as coordinate intersect
    stop("Choose either 1 or 2 for regions.use")
  }
  combined.meta.data <- data.frame(row.names = c(colnames(object.1, colnames(object.2))))
  new.idents <- c()
  for (object in c(object.1, object.2)) {
    old.meta.data <- object[[]]
    if (any(!colnames(x = old.meta.data) %in% colnames(x = combined.meta.data))) {
      cols.to.add <- colnames(x = old.meta.data)[!colnames(x = old.meta.data) %in% colnames(x = combined.meta.data)]
      combined.meta.data[, cols.to.add] <- NA
    }
    i <- sapply(X = old.meta.data, FUN = is.factor)
    old.meta.data[i] <- lapply(X = old.meta.data[i], FUN = as.vector)
    combined.meta.data[rownames(x = old.meta.data), colnames(x = old.meta.data)] <- old.meta.data
    new.idents <- c(new.idents, as.vector(Idents(object = object)))
  }
  names(x = new.idents) <- rownames(x = combined.meta.data)
  new.idents <- factor(x = new.idents)
  if (verbose) {
    message("Constructing merged object")
  }
  counts.1 <- GetAssayData(object = object.1, assay = assay.1, slot = 'counts')[regions.obj1, ]
  counts.2 <- GetAssayData(object = object.2, assay = assay.2, slot = 'counts')[regions.obj2, ]
  rownames(counts.1) <- region.names
  rownames(counts.2) <- region.names
  allcounts <- cbind(counts.1, counts.2)
  assays <- list()
  new.assay <- CreateAssayObject(counts = allcounts, ...)
  assays[[new.assay.name]] <- new.assay
  merged.object <- new(
    Class = 'Seurat',
    assays = assays,
    meta.data = combined.meta.data,
    active.assay = new.assay.name,
    active.ident = new.idents,
    project.name = project,
    version = packageVersion(pkg = 'Seurat')
  )
  return(merged.object)
}

# Generate cut matrix for many regions
#
# Run CutMatrix on multiple regions and add them together.
# Assumes regions are pre-aligned.
#
# @param object A Seurat object
# @param regions A set of GRanges
# @param assay Name of the assay to use
# @param cells Vector of cells to include
# @param verbose Display messages
#' @importFrom Rsamtools TabixFile
MultiRegionCutMatrix <- function(
  object,
  regions,
  assay = NULL,
  cells = NULL,
  verbose = FALSE
) {
  fragment.path <- GetFragments(object = object, assay = assay)
  tabix.file <- TabixFile(file = fragment.path)
  open(con = tabix.file)
  cm.list <- lapply(
    X = seq_along(along.with = regions),
    FUN = function(x) {
      CutMatrix(
        object = object,
        assay = assay,
        tabix.file = tabix.file,
        region = regions[x, ],
        verbose = verbose
      )
    }
  )
  cm <- Reduce(f = `+`, x = cm.list)
  close(con = tabix.file)
  return(cm)
}

# Create cut site pileup matrix
#
# For a set of aligned genomic ranges, find the total number of
# integration sites per cell per base.
#
# @param object A Seurat object
# @param regions A GRanges object
# @param upstream Number of bases to extend upstream
# @param downstream Number of bases to extend downstream
# @param assay Name of the assay to use
# @param cells Which cells to include. If NULL, use all cells
# @param verbose Display messages
#' @importFrom BiocGenerics strand
CreateRegionPileupMatrix <- function(
  object,
  regions,
  upstream = 1000,
  downstream = 1000,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  # extend upstream and downstream from midpoint
  regions <- Extend(
    x = regions,
    upstream = upstream,
    downstream = downstream,
    from.midpoint = TRUE
  )
  # split into strands
  on_plus <- strand(x = regions) == "+" | strand(x = regions) == "*"
  plus.strand <- regions[on_plus, ]
  minus.strand <- regions[!on_plus, ]

  # get cut matrices for each strand
  if (verbose) {
    message("Finding + strand cut sites")
  }
  cut.matrix.plus <- MultiRegionCutMatrix(
    regions = plus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  if (verbose) {
    message("Finding - strand cut sites")
  }
  cut.matrix.minus <- MultiRegionCutMatrix(
    regions = minus.strand,
    object = object,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )

  # reverse minus strand and add together
  full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
  colnames(full.matrix) <- -upstream:downstream
  return(full.matrix)
}

# Apply function to integration sites per base per group
#
# Perform colSums on a cut matrix with cells in the rows
# and position in the columns, for each group of cells
# separately.
#
# @param mat A cut matrix. See \code{\link{CutMatrix}}
# @param groups A vector of group identities, with the name
# of each element in the vector set to the cell name.
# @param fun Function to apply to each group of cells.
# For example, colSums or colMeans.
# @param group.scale.factors Scaling factor for each group. Should
# be computed using the number of cells in the group and the average number of counts
# in the group.
# @param normalize Perform sequencing depth and cell count normalization (default is TRUE)
# @param scale.factor Scaling factor to use. If NULL (default), will use the median normalization
# factor for all the groups.
ApplyMatrixByGroup <- function(
  mat,
  groups,
  fun,
  normalize = TRUE,
  group.scale.factors = NULL,
  scale.factor = NULL
) {
  if (normalize) {
    if (is.null(x = group.scale.factors) | is.null(x = scale.factor)) {
      stop("If normalizing counts, supply group scale factors")
    }
  }
  results <- list()
  all.groups <- unique(x = groups)
  for (i in seq_along(along.with = all.groups)) {
    pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    results[[i]] <- data.frame(
      group = all.groups[[i]],
      count = totals,
      position = as.numeric(colnames(x = mat)),
      stringsAsFactors = FALSE
    )
  }
  coverages <- as.data.frame(x = do.call(what = rbind, args = results), stringsAsFactors = FALSE)
  if (normalize) {
    scale.factor <- SetIfNull(x = scale.factor, y = median(x = group.scale.factors))
    coverages$norm.value <- coverages$count / group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
}

# TabixOutputToDataFrame
#
# Create a single dataframe from list of character vectors
#
# @param reads List of character vectors (the output of \code{\link{scanTabix}})
# @param record.ident Add a column recording which region the reads overlapped with (default TRUE)
#' @importFrom data.table rbindlist
#' @importFrom utils read.table
# @return Returns a data.frame
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  # TODO rewrite this without rbindlist
  df.list <- lapply(X = seq_along(along.with = reads), FUN = function(x) {
    if (length(x = reads[[x]]) == 0) {
      return(NULL)
    }
    df <- read.table(
      file = textConnection(object = reads[[x]]),
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      comment.char = ""
    )
    colnames(x = df) <- c('chr', 'start', 'end', 'cell', 'count')
    if (record.ident) {
      df$ident <- x
    }
    return(df)
  })
  return(rbindlist(l = df.list))
}

# Convert PFMMatrix to
# @param x A PFMatrix
PFMatrixToList <- function(x) {
  if (!requireNamespace('TFBSTools', quietly = TRUE)) {
    stop("Please install TFBSTools. https://www.bioconductor.org/packages/TFBSTools/")
  }
  position.matrix <- TFBSTools::Matrix(x = x)
  name.use <- TFBSTools::name(x = x)
  return(list("matrix" = position.matrix, "name" = name.use))
}

#' Unify genomic ranges
#'
#' Create a unified set of non-overlapping genomic ranges
#' from multiple Seurat objects containing single-cell
#' chromatin data.
#'
#' @param object.list A list of Seurat objects
#' @param mode Function to use when combining genomic ranges. Can be "reduce" (default)
#' or "disjoin". See \code{\link[GenomicRanges]{reduce}} and \code{\link[GenomicRanges]{disjoin}}
#' for more information on these functions.
#' @param sep Separators to use to extract genomic ranges from object row names. To specify different
#' separators for different objects, pass a list of length equal to the length of \code{object.list}.
#'
#' @importFrom GenomicRanges reduce disjoin
#' @export
#' @return Returns a GRanges object
#' @examples
#' UnifyPeaks(object.list = list(atac_small, atac_small))
UnifyPeaks <- function(object.list, mode = 'reduce', sep = c(":", "-")) {
  if (inherits(x = sep, what = "list")) {
    if (length(x = sep) != length(x = object.list)) {
      stop("Must specify separators for each object in the input list")
    }
  } else {
    sep <- rep(x = list(sep), length(x = object.list))
  }
  peak.ranges <- list()
  for (i in seq_along(along.with = object.list)) {
    peak.ranges[[i]] <- StringToGRanges(regions = rownames(object.list[[i]]), sep = sep[[i]])
  }
  peak.ranges <- Reduce(f = c, x = peak.ranges)
  if (mode == 'reduce') {
    return(reduce(x = peak.ranges))
  } else if (mode == 'disjoin') {
    return(disjoin(x = peak.ranges))
  } else {
    stop("Unknown mode requested")
  }
}

#' Subset matrix rows and columns
#'
#' Subset the rows and columns of a matrix by removing
#' rows and columns with less than the specified number of
#' non-zero elements.
#'
#' @param mat A matrix
#' @param min.rows Minimum number of non-zero elements for
#' the row to be retained
#' @param min.cols Minimum number of non-zero elements for
#' the column to be retained
#' @return Returns a matrix
#' @export
#' @importFrom Matrix colSums rowSums
#' @examples
#' SubsetMatrix(mat = volcano)
SubsetMatrix <- function(mat, min.rows = 1, min.cols = 1) {
  rowcount <- rowSums(mat > 0)
  colcount <- colSums(mat > 0)
  keeprows <- rowcount > min.rows
  keepcols <- colcount > min.cols
  return(mat[keeprows, keepcols])
}

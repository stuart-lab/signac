#' @include generics.R
#' @importFrom utils globalVariables
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

#' @importFrom Seurat DefaultAssay Misc
AddToMisc <- function(
  object,
  new.data,
  save.as,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  misc.slot <- Misc(object = object[[assay]]) %||% list()
  if (!inherits(x = misc.slot, what = 'list')) {
    warning("Misc slot already occupied")
  } else{
    misc.slot[[save.as]] <- new.data
    object[[assay]]@misc <- misc.slot
  }
  return(object)
}

globalVariables(names = c('group', 'readcount'), package = 'Signac')
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
#' @examples
#' AverageCounts(atac_small)
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
#' @examples
#' \dontrun{
#' library(EnsDb.Hsapiens.v75)
#' ClosestFeature(
#'   regions = head(rownames(atac_small)),
#'   annotation = EnsDb.Hsapiens.v75,
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
  df$region <- GRangesToString(grange = feature_hits)
  df$distance <- mcols(x = nearest_feature)$distance
  return(df)
}

# Calculate nCount and nFeature
#
# From Seurat
#
# @param object An Assay object
#
# @return A named list with nCount and nFeature
#
#' @importFrom Matrix colSums
#
CalcN <- function(object) {
  if (IsMatrixEmpty(x = GetAssayData(object = object, slot = "counts"))) {
    return(NULL)
  }
  return(list(
    nCount = colSums(x = object, slot = 'counts'),
    nFeature = colSums(x = GetAssayData(object = object, slot = 'counts') > 0)
  ))
}

# Extract delimiter information from a string.
#
# From Seurat
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

#' Compute Tn5 insertion bias
#'
#' Counts the Tn5 insertion frequency for each DNA hexamer.
#'
#' @param object A Seurat object
#' @param genome A BSgenome object
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param region Region to use when assessing bias. Default is human chromosome 1.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom Seurat SetAssayData
#' @export
#'
#' @return Returns a Seurat object
InsertionBias <- function(
  object,
  genome,
  assay = NULL,
  region = 'chr1-1-249250621',
  verbose = TRUE,
  ...
) {
  reads <- GetReadsInRegion(
    object = object,
    region = region,
    assay = assay,
    verbose = verbose,
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
  sequences <- getSeq(x = genome, insertions)
  insertion_hex_freq <- as.matrix(x = table(as.vector(x = sequences)))
  genome_freq <- oligonucleotideFrequency(
    x = getSeq(x = genome, names = 'chr1'),
    width = 6
  )
  insertion_hex_freq <- insertion_hex_freq[names(x = genome_freq), ]
  bias <- insertion_hex_freq / genome_freq
  object <- SetAssayData(object = object, assay = assay, slot = 'bias', new.data = bias)
  return(object)
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
#'
GetIntersectingFeatures <- function(
  object.1,
  object.2,
  assay.1 = NULL,
  assay.2 = NULL,
  distance = 0,
  verbose = TRUE
) {
  regions.1 <- GetAssayData(object = object.1, assay = assay.1, slot = 'ranges')
  regions.2 <- GetAssayData(object = object.2, assay = assay.2, slot = 'ranges')
  if (verbose) {
    message("Intersecting regions across objects")
  }
  region.intersections <- distanceToNearest(x = regions.1, subject = regions.2)
  keep.intersections <- mcols(x = region.intersections)$distance <= distance
  region.intersections <- region.intersections[keep.intersections, ]
  intersect.object1 <- queryHits(x = region.intersections)
  intersect.object2 <- subjectHits(x = region.intersections)
  return(list(intersect.object1, intersect.object2))
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

#' ChunkGRanges
#'
#' Split a genomic ranges object into evenly sized chunks
#'
#' @param granges A GRanges object
#' @param nchunk Number of chunks to split into
#'
#' @return Returns a list of GRanges objects
#' @export
#' @examples
#' ChunkGRanges(blacklist_hg19, n = 10)
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
#' \dontrun{
#' CutMatrix(
#'  object = atac_small,
#'  region = StringToGRanges("chr15-102404831-102407364")
#' )
#' }
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
  all.cells <- cells %||% colnames(x = object)
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
#' @param from.midpoint Extend from the midpoint of the region, rather than the
#' start and end coordinates (default FALSE).
#'
#' @importFrom GenomicRanges strand start end trim
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
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
  assay <- assay %||% DefaultAssay(object = object)
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
    fragment.path <- GetAssayData(object = object, assay = assay, slot = 'fragments')
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

#' CountsInRegion
#'
#' Count reads per cell overlapping a given set of regions
#'
#' @param object A Seurat object
#' @param assay Name of a chromatin assay in the object to use
#' @param regions A GRanges object
#' @param ... Additional arguments passed to \code{\link[IRanges]{findOverlaps}}
#'
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData
#'
#' @export
#' @examples
#' CountsInRegion(
#'   object = atac_small,
#'   assay = 'bins',
#'   regions = blacklist_hg19
#' )
CountsInRegion <- function(
  object,
  assay,
  regions,
  ...
) {
  if (!is(object = object[[assay]], class2 = 'ChromatinAssay')) {
    stop("Must supply a ChromatinAssay")
  }
  obj.granges <- GetAssayData(object = object, assay = assay, slot = 'ranges')
  overlaps <- findOverlaps(query = obj.granges, subject = regions, ...)
  hit.regions <- queryHits(x = overlaps)
  # TODO define a function to return matrix rows corresponding to granges
  # Might not even need it if granges are same order as rows
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
#' @examples
#' ExtractCell(x = "chr1\t1\t10\tatcg\t1")
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
#' @param ... Additional arguments passed to \code{\link{CountsInRegion}}
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData
#'
#' @export
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
  ...
) {
  reads.in.region <- CountsInRegion(
    object = object,
    regions = regions,
    assay = assay,
    ...
  )
  total.reads <- colSums(x = GetAssayData(object = object, assay = assay, slot = 'counts'))
  return(reads.in.region / total.reads)
}

#' Get vector of cell names and associated identity
#' @param object A Seurat object
#' @param group.by Identity class to group cells by
#' @param idents which identities to include
#' @return Returns a named vector
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


# Check if a matrix is empty
#
# From Seurat
#
# Takes a matrix and asks if it's empty (either 0x0 or 1x1 with a value of NA)
#
# @param x A matrix
#
# @return Whether or not \code{x} is empty
#
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
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
#' @examples
#' \dontrun{
#' library(Seurat)
#' counts <- GetAssayData(object = atac_small, assay = 'bins', slot = 'counts')
#' IntersectMatrix(matrix = counts, regions = blacklist_hg19)
#' }
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
#' @param regions Set of query regions. Must be present in rown
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
#' @param regions.use Which regions to use when naming regions in the merged object.
#' Options are:
#' \itemize{
#'  \item{1}: Use the region coordinates from the first object
#'  \item{2}: Use the region coordinates from the second object
#' }
#' @param distance Maximum distance between regions allowed for an intersection to
#' be recorded. Default is 0.
#' @param new.assay.name Name for the merged assay. Default is 'peaks'
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{CreateChromatinAssayObject}}
#'
#' @importFrom Seurat DefaultAssay CreateAssayObject GetAssayData Project
#' @importFrom utils packageVersion
#' @importFrom GenomicRanges granges
#'
#' @export
#' @return Returns a Seurat object
#'
MergeWithRegions <- function(
  object.1,
  object.2,
  assay.1 = NULL,
  assay.2 = NULL,
  regions.use = 1,
  distance = 0,
  new.assay.name = 'ATAC',
  verbose = TRUE,
  ...
) {
  # TODO write as a method merge.ChromatinAssay
  assay.1 <- assay.1 %||% DefaultAssay(object = object.1)
  assay.2 <- assay.2 %||% DefaultAssay(object = object.2)
  intersecting.regions <- GetIntersectingFeatures(
    object.1 = object.1,
    object.2 = object.2,
    assay.1 = assay.1,
    assay.2 = assay.2,
    distance = distance,
    verbose = verbose
  )
  regions.obj1 <- intersecting.regions[[1]]
  regions.obj2 <- intersecting.regions[[2]]
  # TODO add option to keep non-overlapping regions
  if (regions.use == 1) {
    region.names <- rownames(x = object.1)[regions.obj1]
    project <- Project(object = object.1)
    regions <- granges(object.1[[assay.1]])[regions.obj1]
  } else if (regions.use == 2) {
    region.names <-rownames(x = object.2)[regions.obj2]
    project <- Project(object = object.2)
    regions <- granges(object.2[[assay.2]])[regions.obj2]
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
  new.assay <- CreateChromatinAssayObject(counts = allcounts, ranges = regions, ...)
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

####################
### Not exported ###
####################

#' Generate cut matrix for many regions
#'
#' Run CutMatrix on multiple regions and add them together.
#' Assumes regions are pre-aligned.
#'
#' @param object A Seurat object
#' @param regions A set of GRanges
#' @param assay Name of the assay to use
#' @param cells Vector of cells to include
#' @param verbose Display messages
#' @importFrom Rsamtools TabixFile
#' @importFrom Seurat GetAssayData
MultiRegionCutMatrix <- function(
  object,
  regions,
  assay = NULL,
  cells = NULL,
  verbose = FALSE
) {
  fragment.path <- GetAssayData(object = object, assay = assay, slot = 'fragments')
  tabix.file <- TabixFile(file = fragment.path)
  open(con = tabix.file)
  cm.list <- lapply(
    X = 1:length(x = regions),
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

#' Create cut site pileup matrix
#'
#' For a set of aligned genomic ranges, find the total number of
#' integration sites per cell per base.
#'
#' @param object A Seurat object
#' @param regions A GRanges object
#' @param upstream Number of bases to extend regions upstream
#' @param downstream Number of bases to extend regions downstream
#' @param assay Name of the assay to use
#' @param cells Which cells to include. If NULL, use all cells
#' @param verbose Display messages
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

#' Apply function to integration sites per base per group
#'
#' Perform colSums on a cut matrix with cells in the rows
#' and position in the columns, for each group of cells
#' separately.
#'
#' @param mat A cut matrix. See \code{\link{CutMatrix}}
#' @param groups A vector of group identities, with the name
#' of each element in the vector set to the cell name.
#' @param fun Function to apply to each group of cells.
#' For example, colSums or colMeans.
#' @param group.scale.factors Scaling factor for each group. Should
#' be computed using the number of cells in the group and the average number of counts
#' in the group.
#' @param normalize Perform sequencing depth and cell count normalization (default is TRUE)
#' @param scale.factor Scaling factor to use. If NULL (default), will use the median normalization
#' factor for all the groups.
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
    scale.factor <- scale.factor %||% median(x = group.scale.factors)
    coverages$norm.value <- coverages$count / group.scale.factors[coverages$group] * scale.factor
  } else {
    coverages$norm.value <- coverages$count
  }
  return(coverages)
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
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  # TODO rewrite this without rbindlist
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

# Merge multiple rows within a matrix that intersect
# with a single region in another matrix
# First finds the intersecting ranges, then finds rows
# in each matrix that intersect the same region in the
# other matrix then merges the groups of matrix rows
# that intersect the same region of the other matrix.
# The merged row is re-named with the row name of the
# first row that was merged. Rows are returned to their
# original order (missing the merged rows).
# Most of the work done by MergeInternalRows function.
#
# @param mat.a First matrix
# @param ranges.a Ranges associated with the rows of the first matrix
# @param mat.b Second matrix. If NULL, use the first matrix (self-comparison).
# @param ranges.b Ranges associated with the rows of the second matrix.
# @param verbose Display messages
# @return Returns a list of two sparse matrices, A and B,
# with rows that intersect multiple regions in the other matrix merged,
# and two granges objects. If no mat.b is passed, just returns the first
# two elements (matrix and granges).
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicRanges findOverlaps
MergeIntersectingRows <- function(
  mat.a,
  ranges.a,
  mat.b = NULL,
  ranges.b = NULL,
  verbose = TRUE
) {
  if (verbose) {
    message("Finding overlapping ranges")
  }
  if (is.null(x = mat.b)) {
    self = TRUE
  } else {
    self = FALSE
  }
  mat.b <- mat.b %||% mat.a
  ranges.b <- ranges.b %||% ranges.a
  overlaps <- findOverlaps(query = ranges.b, subject = ranges.a)
  b.hits <- queryHits(x = overlaps)
  a.hits <- subjectHits(x = overlaps)
  queryhits.multi.a <-  queryHits(x = overlaps[which(!isUnique(x = a.hits))])
  subjecthits.multi.a <- subjectHits(x = overlaps[which(!isUnique(x = a.hits))])
  queryhits.multi.b <- queryHits(x = overlaps[which(!isUnique(x = b.hits))])
  subjecthits.multi.b <- subjectHits(x = overlaps[which(!isUnique(x = b.hits))])
  a.hits <- ResolveBridge(query = subjecthits.multi.a, subject = queryhits.multi.a)
  b.hits <- ResolveBridge(query = queryhits.multi.b, subject = subjecthits.multi.b)
  multihit.a <- rle(x = b.hits$query)
  multihit.b <- rle(x = a.hits$query)
  if (verbose) {
    message("Merging multiple rows of matrix B that intersect single row in matrix A")
  }
  if (length(x = multihit.b$lengths) != 0) {
    b.mod <- MergeInternalRows(
      mat = mat.b,
      multihit = multihit.b,
      queryhits = a.hits$subject,
      rowranges = ranges.b,
      verbose = verbose
    )
  } else {
    b.mod <- list(mat.b, ranges.b)
  }
  if (self) {
    return(b.mod)
  } else {
    if (verbose) {
      message("Merging multiple rows of matrix A that intersect single row in matrix B")
    }
    if (length(x = multihit.a$lengths) != 0) {
      a.mod <- MergeInternalRows(
        mat = mat.a,
        multihit = multihit.a,
        queryhits = b.hits$subject,
        rowranges = ranges.a,
        verbose = verbose
      )
    } else {
      a.mod <- list(mat.a, ranges.a)
    }
    return(c(a.mod, b.mod))
  }
}

# Resolve cases where a row overlaps two
# ranges that each have multiple overlapping
# rows. This would cause the middle row,
# that overlaps multiple, to become duplicated
# @param query output of queryHits(findOverlaps)
# @param subject output of subjectHits(findOverlaps)
# @return Returns a named list with the modified
# query and subject
ResolveBridge <- function(query, subject) {
  rle.hits <- rle(subject)
  nonunique.hits <- which(rle.hits$lengths > 1)
  runlengths <- rle.hits$lengths[nonunique.hits]
  for (i in seq_along(along.with = nonunique.hits)) {
    q <- nonunique.hits[[i]]
    rl <- runlengths[[i]]
    subject[q+rl-1] <- subject[q+rl]
    query[q+rl-1] <- query[q-1]
  }
  return(list('query' = query, 'subject' = subject))
}

# Set intersecting matrix rows to a common name
#
# ranges.a should correspond to the rows of mat.a
# ranges.b should correspond to the rows of mat.b
# Finds intersecting genomic ranges and renames rows
# mat.b with the name of row in mat.a that it intersects
#
# @param mat.a The first matrix
# @param mat.b The second matrix
# @param ranges.a Genomic ranges corresponding to the rows of mat.a
# @param ranges.b Genomic ranges corresponding to the rows of mat.b
# @param verbose Display messages
# @return Returns mat.b with altered row names
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicRanges findOverlaps
RenameIntersectingRows <- function(
  mat.a,
  mat.b,
  ranges.a,
  ranges.b,
  verbose = TRUE
) {
  overlaps <- findOverlaps(query = ranges.a, subject = ranges.b)
  a.hits <- queryHits(x = overlaps)
  b.hits <- subjectHits(x = overlaps)
  rownames(mat.b)[b.hits] <- rownames(mat.a)[a.hits]
  return(mat.b)
}

# Condense overlapping GenomicRanges
#
# Take the union of two sets of genomic ranges,
# then split back into the ranges that were present
# in each original set.
#
# @param ranges.a First set of GenomicRanges
# @param ranges.b Second set of GenomicRanges
# @return Returns a list of two GenomicRanges objects
#' @importFrom GenomicRanges union findOverlaps
#' @importFrom S4Vectors queryHits
CondenseOverlappingGRanges <- function(
  ranges.a,
  ranges.b
) {
  grange.union <- union(x = granges(x = ranges.a), y = granges(x = ranges.b))
  a.overlap.union <- findOverlaps(
    query = grange.union,
    subject = ranges.a,
    ignore.strand = TRUE,
    select = 'first' # select = 'first' is essential to have match matrix rows
  )
  a.overlap.union <- a.overlap.union[!is.na(x = a.overlap.union)]
  b.overlap.union <- findOverlaps(
    query = grange.union,
    subject = ranges.b,
    ignore.strand = TRUE,
    select = 'first'
  )
  b.overlap.union <- b.overlap.union[!is.na(x = b.overlap.union)]
  a.subset <- ranges.a[a.overlap.union]
  b.subset <- ranges.b[b.overlap.union]
  # TODO still doesn't match results from matrix row merge
  return(list(a.subset, b.subset))
}

# Merge rows of matrix based on intersection with a second matrix.

# TODO this can be slow if there are a lot of rows to merge
#      because we repeatedly extract and sum rows of a sparse matrix.
#      could try to make parallel, or think about other ways of doing it.

# @param mat A sparse matrix
# @param multihit a run-length encoded list describing the number of rows to be merged in each merge step
# @param queryhits Row indices to be merged
# @param rowranges GenomicRanges associated with matrix rows
# @param verbose Display progress
# @return Returns a list containing a sparse matrix,
# with target rows merged in-place (order same but condensed rows),
# and a set of GenomicRanges
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix colSums
MergeInternalRows <- function(
  mat,
  multihit,
  queryhits,
  rowranges,
  verbose = TRUE
) {
  mat.rows <- rownames(x = mat)
  i <- 1
  newmat <- list()
  todelete <- c()
  todelete.range <- c()
  if (verbose) {
    pb <- txtProgressBar(
      min = 1,
      max = length(x = multihit$lengths),
      style = 3,
      file = stderr()
    )
  }
  while (i < length(x = multihit$lengths)) {
    rowrun <- multihit$lengths[[i]]
    # find the matrix rows to merge
    rowindex <- mat.rows[queryhits[i:(i+rowrun-1)]]
    if (rowrun > 1) {
      rangeindex <- queryhits[(i+1):(i+rowrun-1)]
    }
    # merge rows and add to list, name will become the name of the row
    if (length(x = rowindex) > 1) {
      newmat[[rowindex[[1]]]] <- colSums(mat[rowindex, ])
    } else {
      newmat <- NULL
    }
    # record which rows need to be removed
    todelete <- c(todelete, rowindex)
    todelete.range <- c(todelete.range, rangeindex)
    i <- i + rowrun
    if (verbose) setTxtProgressBar(pb = pb, value = i)
  }
  # contruct matrix
  if (verbose) {
    message("\nBinding matrix rows")
  }
  merged.mat <- Reduce(f = rbind, x = newmat)
  rownames(merged.mat) <- names(newmat)
  merged.mat <- as(object = merged.mat, Class = 'dgCMatrix')
  # remove rows from A that were merged
  tokeep <- setdiff(mat.rows, todelete)
  mat.mod <- mat[tokeep, ]
  rowranges[todelete.range] <- NULL
  # add new merged rows to A
  mat.mod <- rbind(mat.mod, merged.mat)
  # put back in original order
  rows.order <- mat.rows[mat.rows %in% rownames(x = mat.mod)]
  mat.mod <- mat.mod[rows.order, ]
  return(list(mat.mod, rowranges))
}

# Check if fragment file exists and is indexed
# @param f A fragment file path
# @return Returns TRUE is present and indexed, otherwise FALSE
ValidFragments <- function(f) {
  index.file <- paste0(f, ".tbi")
  if (all(file.exists(f, index.file))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

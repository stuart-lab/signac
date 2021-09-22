#' @include generics.R
#' @importFrom utils globalVariables
NULL


#' Add chromatin module
#'
#' Compute chromVAR deviations for groups of peaks. The goal of this function is
#' similar to that of \code{\link[Seurat]{AddModuleScore}} except that it is
#' designed for single-cell chromatin data. The chromVAR deviations for each
#' group of peaks will be added to the object metadata.
#'
#' @param object A Seurat object
#' @param features A named list of features to include in each module. The name
#' of each element in the list will be used to name the modules computed, which
#' will be stored in the object metadata.
#' @param genome A BSgenome object
#' @param assay Name of assay to use. If NULL, use the default assay.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{RunChromVAR}
#'
#' @return Returns a Seurat object
#'
#' @importFrom fastmatch fmatch
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat DefaultAssay GetAssayData AddMetaData
#'
#' @export
#' @concept utilities
AddChromatinModule <- function(
  object,
  features,
  genome,
  assay = NULL,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }

  # first find index of each feature
  feat.idx <- sapply(X = features, FUN = fmatch, rownames(x = object[[assay]]))
  j <- sapply(X = seq_along(along.with = features), FUN = function(x) {
    rep(x = x, length(x = features[[x]]))
  })

  # construct sparse matrix with features
  mat <- sparseMatrix(
    i = unlist(x = feat.idx, use.names = FALSE),
    j = unlist(x = j, use.names = FALSE),
    x = 1,
    dims = c(nrow(x = object[[assay]]), length(x = features))
  )
  rownames(x = mat) <- rownames(x = object[[assay]])
  colnames(x = mat) <- names(x = features)

  # run chromVAR
  cv <- RunChromVAR(
    object = object[[assay]],
    motif.matrix = mat,
    genome = genome,
    verbose = verbose,
    ...
  )

  # add module scores to metadata
  chromvar.data <- GetAssayData(object = cv, slot = "data")
  object <- AddMetaData(
    object = object,
    metadata = as.data.frame(x = t(x = chromvar.data))
  )
  return(object)
}

globalVariables(names = c("group", "readcount"), package = "Signac")
#' Average Counts
#'
#' Compute the mean counts per group of cells for a given assay
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. Default is the active assay
#' @param group.by Grouping variable to use. Default is the active identities
#' @param verbose Display messages
#'
#' @importFrom Seurat DefaultAssay Idents GetAssayData
#' @importFrom dplyr group_by summarize
#' @export
#' @concept utilities
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
  # pull nCount_ column
  col.use <- paste0("nCount_", assay)
  total.df <- object[[col.use]]
  colnames(x = total.df) <- "readcount"
  total.df$cell <- rownames(x = total.df)
  total.df$group <- group.by[total.df$cell]
  total.df <- group_by(total.df, group)
  if (verbose) {
    message("Computing average counts per group")
  }
  group.means <- summarize(.data = total.df, mn = mean(x = readcount))
  results <- group.means$mn
  names(x = results) <- group.means$group
  return(results)
}

#' Accessible peaks
#'
#' Find accessible peaks in a set of cells
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param idents A set of identity classes to find accessible peaks for
#' @param cells A vector of cells to find accessible peaks for
#' @param min.cells Minimum number of cells with the peak accessible (>0 counts)
#' for the peak to be called accessible
#' @export
#' @concept utilities
#' @importFrom Seurat WhichCells DefaultAssay GetAssayData
#' @importFrom Matrix rowSums
#' @return Returns a vector of peak names
AccessiblePeaks <- function(
  object,
  assay = NULL,
  idents = NULL,
  cells = NULL,
  min.cells = 10
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  cells <- SetIfNull(x = cells, y = WhichCells(object, idents = idents))
  open.peaks <- GetAssayData(
    object = object,
    assay = assay,
    slot = "counts"
  )[, cells]
  peaks <- names(x = which(x = rowSums(x = open.peaks > 0) > min.cells))
  return(peaks)
}

#' Cells per group
#'
#' Count the number of cells in each group
#'
#' @param object A Seurat object
#' @param group.by A grouping variable. Default is the active identities
#' @importFrom Seurat Idents
#' @export
#' @concept utilities
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
  cells.per.group <- table(cellgroups, useNA = "always")
  lut <- as.vector(x = cells.per.group)
  names(x = lut) <- names(x = cells.per.group)
  return(lut)
}

#' Closest Feature
#'
#' Find the closest feature to a given set of genomic regions
#'
#' @param object A Seurat object
#' @param regions A set of genomic regions to query
#' @param annotation A GRanges object containing annotation information. If
#' NULL, use the annotations stored in the object.
#' @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#'
#' @importMethodsFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors subjectHits mcols
#' @importFrom methods is
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomeInfoDb dropSeqlevels
#' @return Returns a dataframe with the name of each region, the closest feature
#' in the annotation, and the distance to the feature.
#' @export
#' @concept utilities
#' @examples
#' \donttest{
#' ClosestFeature(
#'   object = atac_small,
#'   regions = head(granges(atac_small))
#' )
#' }
ClosestFeature <- function(
  object,
  regions,
  annotation = NULL,
  ...
) {
  if (!is(object = regions, class2 = 'GRanges')) {
    regions <- StringToGRanges(regions = regions, ...)
  }
  if (inherits(x = object, what = "Seurat")) {
    # running on Seurat object, extract the assay
    assay <- DefaultAssay(object = object)
    object <- object[[assay]]
  }
  if (!inherits(x = object, what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  if (length(x = regions) == 0) {
    stop("No query regions supplied")
  }
  annotation <- SetIfNull(x = annotation, y = Annotation(object = object))
  missing_seqlevels <- setdiff(
    x = seqlevels(x = regions), y = seqlevels(x = annotation)
  )
  if (length(x = missing_seqlevels) > 0) {
    warning(
      "The following seqlevels present in query regions are not present\n ",
      "in the supplied gene annotations and will be removed: ",
      paste(missing_seqlevels, collapse = ", ")
    )
    regions <- dropSeqlevels(
      x = regions,
      value = missing_seqlevels,
      pruning.mode = "coarse"
    )
    if (length(x = regions) == 0) {
      stop("None of the supplied regions were found in the supplied annotation")
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

#' Create gene activity matrix
#'
#' Compute counts per cell in gene body and promoter region.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use. If NULL, use the default assay
#' @param features Genes to include. If NULL, use all protein-coding genes in
#' the annotations stored in the object
#' @param extend.upstream Number of bases to extend upstream of the TSS
#' @param extend.downstream Number of bases to extend downstream of the TTS
#' @param biotypes Gene biotypes to include. If NULL, use all biotypes in the
#' gene annotation.
#' @param max.width Maximum allowed gene width for a gene to be quantified.
#' Setting this parameter can avoid quantifying extremely long transcripts that
#' can add a relatively long amount of time. If NULL, do not filter genes based
#' on width.
#' @param verbose Display messages
#' @param ... Additional options passed to \code{\link{FeatureMatrix}}
#'
#' @return Returns a sparse matrix
#'
#' @concept utilities
#' @export
#' @importFrom Seurat DefaultAssay
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' fragments <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   validate.fragments = FALSE
#' )
#' Fragments(atac_small) <- fragments
#' GeneActivity(atac_small)
GeneActivity <- function(
  object,
  assay = NULL,
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = "protein_coding",
  max.width = 500000,
  verbose = TRUE,
  ...
) {
  if (!is.null(x = features)) {
    if (length(x = features) == 0) {
      stop("Empty list of features provided")
    }
  }
  # collapse to longest protein coding transcript
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  annotation <- Annotation(object = object[[assay]])
  if (length(x = annotation) == 0) {
    stop("No gene annotations present in object")
  }
  if (verbose) {
    message("Extracting gene coordinates")
  }
  transcripts <- CollapseToLongestTranscript(ranges = annotation)
  if (!is.null(x = biotypes)) {
    transcripts <- transcripts[transcripts$gene_biotype %in% biotypes]
    if (length(x = transcripts) == 0) {
      stop("No genes remaining after filtering for requested biotypes")
    }
  }

  # filter genes if provided
  if (!is.null(x = features)) {
    transcripts <- transcripts[transcripts$gene_name %in% features]
    if (length(x = transcripts) == 0) {
      stop("None of the requested genes were found in the gene annotation")
    }
  }
  if (!is.null(x = max.width)) {
    transcript.keep <- which(x = width(x = transcripts) < max.width)
    transcripts <- transcripts[transcript.keep]
    if (length(x = transcripts) == 0) {
      stop("No genes remaining after filtering for max.width")
    }
  }

  # extend to include promoters
  transcripts <- Extend(
    x = transcripts,
    upstream = extend.upstream,
    downstream = extend.downstream
  )

  # quantify
  frags <- Fragments(object = object[[assay]])
  if (length(x = frags) == 0) {
    stop("No fragment information found for requested assay")
  }
  cells <- colnames(x = object[[assay]])
  counts <- FeatureMatrix(
    fragments = frags,
    features = transcripts,
    cells = cells,
    verbose = verbose,
    ...
  )

  # set row names
  gene.key <- transcripts$gene_name
  names(x = gene.key) <- GRangesToString(grange = transcripts)
  rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
  counts <- counts[rownames(x = counts) != "", ]

  return(counts)
}

#' Extract genomic ranges from EnsDb object
#'
#' Pulls the transcript information for all chromosomes from an EnsDb object.
#' This wraps \code{\link[biovizBase]{crunch}} and applies the extractor
#' function to all chromosomes present in the EnsDb object.
#'
#' @param ensdb An EnsDb object
#' @param standard.chromosomes Keep only standard chromosomes
#' @param biotypes Biotypes to keep
#' @param verbose Display messages
#'
#' @importFrom GenomeInfoDb keepStandardChromosomes seqinfo
#' @concept utilities
#' @export
GetGRangesFromEnsDb <- function(
  ensdb,
  standard.chromosomes = TRUE,
  biotypes = c("protein_coding", "lincRNA", "rRNA", "processed_transcript"),
  verbose = TRUE
) {
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = ensdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  }

  # extract genes from each chromosome
  if (verbose) {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      biovizBase::crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype"))
    })
  } else {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      suppressMessages(expr = biovizBase::crunch(
        obj = ensdb,
        which = whole.genome[x],
        columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")))
    })
  }

  # combine
  tx <- do.call(what = c, args = tx)
  tx <- tx[tx$gene_biotype %in% biotypes]
  return(tx)
}

#' Find transcriptional start sites
#'
#' Get the TSS positions from a set of genomic ranges containing gene positions.
#' Ranges can contain exons, introns, UTRs, etc, rather than the whole
#' transcript. Only protein coding gene biotypes are included in output.
#'
#' @param ranges A GRanges object containing gene annotations.
#' @param biotypes Gene biotypes to include. If NULL, use all biotypes in the
#' supplied gene annotation.
#' @importFrom GenomicRanges resize
#' @importFrom S4Vectors mcols
#' @export
#' @concept utilities
GetTSSPositions <- function(ranges, biotypes = "protein_coding") {
  if (!("gene_biotype" %in% colnames(x = mcols(x = ranges)))) {
    stop("Gene annotation does not contain gene_biotype information")
  }
  if (!is.null(x = biotypes)){
    ranges <- ranges[ranges$gene_biotype == "protein_coding"]
  }
  gene.ranges <- CollapseToLongestTranscript(ranges = ranges)
  # shrink to TSS position
  tss <- resize(gene.ranges, width = 1, fix = 'start')
  return(tss)
}

#' Find intersecting regions between two objects
#'
#' Intersects the regions stored in the rownames of two objects and
#' returns a vector containing the names of rows that intersect
#' for each object. The order of the row names return corresponds
#' to the intersecting regions, i.e. the nth feature of the first vector
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
#' @param distance Maximum distance between regions allowed for an intersection
#' to be recorded. Default is 0.
#' @param verbose Display messages
#'
#' @importMethodsFrom GenomicRanges distanceToNearest
#' @importFrom S4Vectors subjectHits queryHits mcols
#' @importFrom Seurat DefaultAssay
#' @export
#' @concept utilities
#' @return Returns a list of two character vectors containing the row names
#' in each object that overlap each other.
#' @examples
#' GetIntersectingFeatures(
#'   object.1 = atac_small,
#'   object.2 = atac_small,
#'   assay.1 = 'peaks',
#'   assay.2 = 'bins'
#' )
GetIntersectingFeatures <- function(
  object.1,
  object.2,
  assay.1 = NULL,
  assay.2 = NULL,
  distance = 0,
  verbose = TRUE
) {
  assay.1 <- SetIfNull(x = assay.1, y = DefaultAssay(object = object.1))
  assay.2 <- SetIfNull(x = assay.2, y = DefaultAssay(object = object.2))
  if (!inherits(x = object.1[[assay.1]], what = "ChromatinAssay")) {
    stop("Requested assay in object 1 is not a ChromatinAssay.")
  }
  if (!inherits(x = object.2[[assay.2]], what = "ChromatinAssay")) {
    stop("Requested assay in object 2 is not a ChromatinAssay")
  }
  regions.1 <- GetAssayData(object = object.1, assay = assay.1, slot = "ranges")
  regions.2 <- GetAssayData(object = object.2, assay = assay.2, slot = "ranges")
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

#' String to GRanges
#'
#' Convert a genomic coordinate string to a GRanges object
#'
#' @param regions Vector of genomic region strings
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param ... Additional arguments passed to
#' \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}
#' @return Returns a GRanges object
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate
#' @examples
#' regions <- c('chr1-1-10', 'chr2-12-3121')
#' StringToGRanges(regions = regions)
#' @export
#' @concept utilities
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}

#' GRanges to String
#'
#' Convert GRanges object to a vector of strings
#'
#' @param grange A GRanges object
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @importMethodsFrom GenomicRanges start end seqnames
#' @examples
#' GRangesToString(grange = blacklist_hg19)
#' @return Returns a character vector
#' @export
#' @concept utilities
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
#' @importFrom GenomicRanges trim
#' @importFrom BiocGenerics start strand end width
#' @importMethodsFrom GenomicRanges strand start end width
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
#' @concept utilities
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
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  ranges(x = x) <- IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}

#' Get cells in a region
#'
#' Extract cell names containing reads mapped within a given genomic region
#'
#' @param tabix Tabix object
#' @param region A string giving the region to extract from the fragments file
#' @param cells Vector of cells to include in output. If NULL, include all cells
#'
#' @importFrom Rsamtools scanTabix
#' @importFrom methods is
#' @importFrom fastmatch fmatch
#' @export
#' @concept utilities
#' @return Returns a list
#' @examples
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' GetCellsInRegion(tabix = fpath, region = "chr1-10245-762629")
GetCellsInRegion <- function(tabix, region, cells = NULL) {
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  reads <- scanTabix(file = tabix, param = region)
  reads <- lapply(X = reads, FUN = ExtractCell)
  if (!is.null(x = cells)) {
    reads <- sapply(X = reads, FUN = function(x) {
      x <- x[fmatch(x = x, table = cells, nomatch = 0L) > 0L]
      if (length(x = x) == 0) {
        return(NULL)
      } else {
        return(x)
      }
    })
  }
  return(reads)
}

#' Counts in region
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
#' @concept utilities
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
  ...
) {
  if (!is(object = object[[assay]], class2 = "ChromatinAssay")) {
    stop("Must supply a ChromatinAssay")
  }
  obj.granges <- GetAssayData(object = object, assay = assay, slot = "ranges")
  overlaps <- findOverlaps(query = obj.granges, subject = regions, ...)
  hit.regions <- queryHits(x = overlaps)
  data.matrix <- GetAssayData(
    object = object, assay = assay, slot = "counts"
  )[hit.regions, , drop = FALSE]
  return(colSums(data.matrix))
}

#' Fraction of counts in a genomic region
#'
#' Find the fraction of counts per cell that overlap a given set of genomic
#' ranges
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param regions A GRanges object containing a set of genomic regions
#' @param ... Additional arguments passed to \code{\link{CountsInRegion}}
#' @importFrom Matrix colSums
#' @importFrom Seurat GetAssayData DefaultAssay
#'
#' @export
#' @concept utilities
#' @return Returns a numeric vector
#' @examples
#' FractionCountsInRegion(
#'   object = atac_small,
#'   assay = 'bins',
#'   regions = blacklist_hg19
#' )
FractionCountsInRegion <- function(
  object,
  regions,
  assay = NULL,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  reads.in.region <- CountsInRegion(
    object = object,
    regions = regions,
    assay = assay,
    ...
  )
  total.reads <- colSums(x = GetAssayData(
    object = object, assay = assay, slot = "counts"
  ))
  return(reads.in.region / total.reads)
}

#' Intersect genomic coordinates with matrix rows
#'
#' Remove or retain matrix rows that intersect given genomic regions
#'
#' @param matrix A matrix with genomic regions in the rows
#' @param regions A set of genomic regions to intersect with regions in the
#' matrix. Either a vector of strings encoding the genomic coordinates, or a
#' GRanges object.
#' @param invert Discard rows intersecting the genomic regions supplied, rather
#' than retain.
#' @param sep A length-2 character vector containing the separators to be used
#' for extracting genomic coordinates from a string. The first element will be
#' used to separate the chromosome name from coordinates, and the second element
#' used to separate start and end coordinates.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link[IRanges]{findOverlaps}}
#'
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#'
#' @export
#' @concept utilities
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
  if (is(object = regions, class2 = "character")) {
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

#' Get gene coordinates
#'
#' Extract the coordinates of the longest transcript for a gene stored in the
#' annotations within an object.
#'
#' @param object A Seurat object
#' @param gene Name of a gene to extract
#' @param assay Name of assay to use
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @export
#' @concept utilities
#' @examples
#' LookupGeneCoords(atac_small, gene = "MIR1302-10")
LookupGeneCoords <- function(object, gene, assay = NULL) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  annotations <- Annotation(object = object[[assay]])
  isgene <- annotations$gene_name == gene
  isgene <- !is.na(x = isgene) & isgene
  annot.sub <- annotations[isgene]
  if (length(x = annot.sub) == 0) {
    return(NULL)
  } else {
    gr <- GRanges(seqnames = as.character(x = seqnames(x = annot.sub))[[1]],
                  ranges = IRanges(start = min(start(x = annot.sub)),
                                   end = max(end(x = annot.sub))))
    return(gr)
  }
}

#' Match DNA sequence characteristics
#'
#' Return a vector if genomic regions that match the distribution of a set of
#' query regions for any given set of characteristics, specified in the input
#' \code{meta.feature} dataframe.
#'
#' For each requested feature to match, a density
#' distribution is estimated using the \code{\link[stats]{density}} function,
#' and a set of weights for each feature in the dataset estimated based on the
#' density distribution. If multiple features are to be matched (for example,
#' GC content and overall accessibility), a joint density distribution is then
#' computed by multiplying the individual feature weights. A set of features
#' with characteristics matching the query regions is then selected using the
#' \code{\link[base]{sample}} function, with the probability of randomly
#' selecting each feature equal to the joint density distribution weight.
#'
#' @param meta.feature A dataframe containing DNA sequence information for
#' features to choose from
#' @param query.feature A dataframe containing DNA sequence information for
#' features to match.
#' @param n Number of regions to select, with characteristics matching the query
#' @param features.match Which features of the query to match when selecting a
#' set of regions. A vector of column names present in the feature metadata can
#' be supplied to match multiple characteristics at once. Default is GC content.
#' @param verbose Display messages
#' @param ... Arguments passed to other functions
#' @return Returns a character vector
#'
#' @importFrom stats density approx
#' @export
#' @concept utilities
#' @concept motifs
#' @examples
#' metafeatures <- Seurat::GetAssayData(
#'   object = atac_small[['peaks']], slot = 'meta.features'
#' )
#' query.feature <- metafeatures[1:10, ]
#' features.choose <- metafeatures[11:nrow(metafeatures), ]
#' MatchRegionStats(
#'   meta.feature = features.choose,
#'   query.feature = query.feature,
#'   features.match = "percentile",
#'   n = 10
#' )
MatchRegionStats <- function(
  meta.feature,
  query.feature,
  features.match = c("GC.percent"),
  n = 10000,
  verbose = TRUE,
  ...
) {
  if (!inherits(x = meta.feature, what = 'data.frame')) {
    stop("meta.feature should be a data.frame")
  }
  if (!inherits(x = query.feature, what = "data.frame")) {
    stop("query.feature should be a data.frame")
  }
  if (length(x = features.match) == 0) {
    stop("Must supply at least one sequence characteristic to match")
  }
  if (nrow(x = meta.feature) < n) {
    n <- nrow(x = meta.feature)
    warning("Requested more features than present in supplied data.
            Returning ", n, " features")
  }
  # features.choose <- meta.feature[choosefrom, ]
  for (i in seq_along(along.with = features.match)) {
    featmatch <- features.match[[i]]
    if (!(featmatch %in% colnames(x = query.feature))) {
      if (featmatch == "GC.percent") {
        stop("GC.percent not present in meta.features.",
             " Run RegionStats to compute GC.percent for each feature.")
      } else {
        stop(featmatch, " not present in meta.features")
      }
    }
    if (verbose) {
      message("Matching ", featmatch, " distribution")
    }
    density.estimate <- density(
      x = query.feature[[featmatch]], kernel = "gaussian", bw = 1
    )
    weights <- approx(
      x = density.estimate$x,
      y = density.estimate$y,
      xout = meta.feature[[featmatch]],
      yright = 0.0001,
      yleft = 0.0001
    )$y
    if (i > 1) {
      feature.weights <- feature.weights * weights
    } else {
      feature.weights <- weights
    }
  }
  feature.select <- sample.int(
    n = nrow(x = meta.feature),
    size = n,
    prob = feature.weights
  )
  feature.select <- rownames(x = meta.feature)[feature.select]
  return(feature.select)
}

#' Unify genomic ranges
#'
#' Create a unified set of non-overlapping genomic ranges
#' from multiple Seurat objects containing single-cell
#' chromatin data.
#'
#' @param object.list A list of Seurat objects or ChromatinAssay objects
#' @param mode Function to use when combining genomic ranges. Can be "reduce"
#' (default) or "disjoin".
#' See \code{\link[GenomicRanges]{reduce}}
#' and \code{\link[GenomicRanges]{disjoin}}
#' for more information on these functions.
#' @importFrom GenomicRanges reduce disjoin
#' @export
#' @concept utilities
#' @concept preprocessing
#' @return Returns a GRanges object
#' @examples
#' UnifyPeaks(object.list = list(atac_small, atac_small))
UnifyPeaks <- function(object.list, mode = "reduce") {
  peak.ranges <- list()
  for (i in seq_along(along.with = object.list)) {
    peak.ranges[[i]] <- granges(x = object.list[[i]])
  }
  peak.ranges <- Reduce(f = c, x = peak.ranges)
  if (mode == "reduce") {
    return(reduce(x = peak.ranges))
  } else if (mode == "disjoin") {
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
#' @param max.row.val Maximum allowed value in a row for the
#' row to be retained. If NULL, don't set any limit.
#' @param max.col.val Maximum allowed value in a column for
#' the column to be retained. If NULL, don't set any limit.
#'
#' @return Returns a matrix
#' @export
#' @concept utilities
#' @importFrom Matrix colSums rowSums
#' @examples
#' SubsetMatrix(mat = volcano)
SubsetMatrix <- function(
  mat,
  min.rows = 1,
  min.cols = 1,
  max.row.val = 10,
  max.col.val = NULL
) {
  rowcount <- rowSums(mat > 0)
  colcount <- colSums(mat > 0)
  keeprows <- rowcount > min.rows
  keepcols <- colcount > min.cols
  if (!is.null(x = max.row.val)) {
    rowmax <- apply(X = mat, MARGIN = 1, FUN = max)
    keeprows <- keeprows & (rowmax < max.row.val)
  }
  if (!is.null(x = max.col.val)) {
    colmax <- apply(X = mat, MARGIN = 2, FUN = max)
    keepcols <- keepcols & (colmax < max.col.val)
  }
  return(mat[keeprows, keepcols])
}

#### Not exported ####

#' @importFrom IRanges isDisjoint
NonOverlapping <- function(x, all.features) {
  # x is list of assays
  diff.features <- names(x = all.features[all.features < length(x = x)])
  if (length(x = diff.features) == 0) {
    return(TRUE)
  } else {
    diff.ranges <- StringToGRanges(regions = diff.features)
    return(isDisjoint(x = diff.ranges))
  }
}

#' @importFrom Matrix sparseMatrix
AddMissingCells <- function(x, cells) {
  # add columns with zeros for cells not in matrix
  missing.cells <- setdiff(x = cells, y = colnames(x = x))
  if (!(length(x = missing.cells) == 0)) {
    null.mat <- sparseMatrix(
      i = c(),
      j = c(),
      dims = c(nrow(x = x), length(x = missing.cells))
    )
    rownames(x = null.mat) <- rownames(x = x)
    colnames(x = null.mat) <- missing.cells
    x <- cbind(x, null.mat)
  }
  x <- x[, cells, drop = FALSE]
  return(x)
}

#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom Matrix Diagonal tcrossprod rowSums
AverageCountMatrix <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL
) {
  assay = SetIfNull(x = assay, y = DefaultAssay(object = object))
  countmatrix <- GetAssayData(object = object[[assay]], slot = "counts")
  ident.matrix <- BinaryIdentMatrix(
    object = object,
    group.by = group.by,
    idents = idents
  )
  collapsed.counts <- tcrossprod(x = countmatrix, y = ident.matrix)
  avg.counts <- tcrossprod(
    x = collapsed.counts,
    y = Diagonal(x = 1 / rowSums(x = ident.matrix))
  )
  return(as.matrix(x = avg.counts))
}

# Create binary cell x class matrix of group membership
#' @importFrom Matrix sparseMatrix
BinaryIdentMatrix <- function(object, group.by = NULL, idents = NULL) {
  group.idents <- GetGroups(object = object, group.by = group.by, idents = idents)
  cell.idx <- seq_along(along.with = names(x = group.idents))
  unique.groups <- as.character(x = unique(x = group.idents))
  ident.idx <- seq_along(along.with = unique.groups)
  names(x = ident.idx) <- unique.groups
  ident.matrix <- sparseMatrix(
    i = ident.idx[as.character(x = group.idents)],
    j = cell.idx,
    x = 1
  )
  colnames(x = ident.matrix) <- names(x = group.idents)
  rownames(x = ident.matrix) <- unique.groups
  ident.matrix <- as(object = ident.matrix, Class = "dgCMatrix")
  return(ident.matrix)
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
    nCount = colSums(x = object, slot = "counts"),
    nFeature = colSums(x = GetAssayData(object = object, slot = "counts") > 0)
  ))
}

globalVariables(
  names = c("start", "end", "seqnames", "strand", "gene_biotype"),
  package = "Signac"
)
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# Chunk GRanges
#
# Split a genomic ranges object into evenly sized chunks
# @param granges A GRanges object
# @param nchunk Number of chunks to split into
#
# @return Returns a list of GRanges objects
# @examples
# ChunkGRanges(blacklist_hg19, n = 10)
ChunkGRanges <- function(granges, nchunk) {
  if (length(x = granges) < nchunk) {
    nchunk <- length(x = granges)
  }
  chunksize <- as.integer(x = (length(granges) / nchunk))
  range.list <- sapply(X = seq_len(length.out = nchunk), FUN = function(x) {
    chunkupper <- (x * chunksize)
    if (x == 1) {
      chunklower <- 1
    } else {
      chunklower <- ((x - 1) * chunksize) + 1
    }
    if (x == nchunk) {
      chunkupper <- length(x = granges)
    }
    return(granges[chunklower:chunkupper])
  })
  return(range.list)
}

# Extract cell
#
# Extract cell barcode from list of tab delimited character
# vectors (output of \code{\link{scanTabix}})
#
# @param x List of character vectors
# @return Returns a string
#' @importFrom stringi stri_split_fixed
ExtractCell <- function(x) {
  if (length(x = x) == 0) {
    return(NULL)
  } else {
    x <- stri_split_fixed(str = x, pattern = "\t")
    n <- length(x = x)
    x <- unlist(x = x)
    return(unlist(x = x)[5 * (1:n) - 1])
  }
}

# Run groupCommand for the first n lines, convert the cell barcodes in the file
# to the cell names that appear in the fragment object, and subset the output to
# cells present in the fragment object
#
# Every cell in the fragment file will be present in the output dataframe. If
# the cell information is not set, every cell barcode that appears in the first
# n lines will be present.
#
# @param fragments A Fragment object
# @param n Number of lines to read from the beginning of the fragment file
# @param verbose Display messages
#
# @return Returns a data.frame
ExtractFragments <- function(fragments, n = NULL, verbose = TRUE) {
  fpath <- GetFragmentData(object = fragments, slot = "path")
  if (isRemote(x = fpath)) {
    stop("Remote fragment files not supported")
  }
  fpath <- normalizePath(path = fpath, mustWork = TRUE)
  cells <- GetFragmentData(object = fragments, slot = "cells")
  if (!is.null(x = cells)) {
    cells.use <- as.character(x = cells)
  } else {
    cells.use <- NULL
  }
  verbose <- as.logical(x = verbose)
  n <- SetIfNull(x = n, y = 0)
  n <- as.numeric(x = n)
  n <- round(x = n, digits = 0)
  counts <- groupCommand(
    fragments = fpath,
    some_whitelist_cells = cells.use,
    max_lines = n,
    verbose = verbose
  )
  # convert cell names
  if (!is.null(x = cells)) {
    # every cell will be present in the output, even if 0 counts
    converter <- names(x = cells)
    names(x = converter) <- cells
    counts$CB <- converter[counts$CB]
  }
  return(counts)
}

# convert region argument to genomic coordinates
# region can be a string, name of a gene, or GRanges object
FindRegion <- function(
  object,
  region,
  sep = c("-", "-"),
  assay = NULL,
  extend.upstream = 0,
  extend.downstream = 0
) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  )
  )
  return(region)
}

# GetReadsInRegion
#
# Extract reads for each cell within a given genomic region or set of regions
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A genomic region, specified as a string in the format
# 'chr:start-end'. Can be a vector of regions.
# @param tabix.file A TabixFile object.
# @param cells Cells to include. Default is all cells present in the object.
# @param verbose Display messages
# @param ... Additional arguments passed to \code{\link{StringToGRanges}}
#
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom Seurat Idents
#' @importFrom fastmatch fmatch
#
# @return Returns a data frame
GetReadsInRegion <- function(
  cellmap,
  region,
  tabix.file,
  cells = NULL,
  verbose = TRUE,
  ...
) {
  file.to.object <- names(x = cellmap)
  names(x = file.to.object) <- cellmap

  if (verbose) {
    message("Extracting reads in requested region")
  }
  if (!is(object = region, class2 = "GRanges")) {
    region <- StringToGRanges(regions = region, ...)
  }
  # remove regions that aren't in the fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region),
    y = seqnamesTabix(file = tabix.file)
  )
  region <- keepSeqlevels(
    x = region,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  reads <- scanTabix(file = tabix.file, param = region)
  reads <- TabixOutputToDataFrame(reads = reads)
  reads <- reads[
    fmatch(x = reads$cell, table = cellmap, nomatch = 0L) > 0,
  ]
  # convert cell names to match names in object
  reads$cell <- file.to.object[reads$cell]
  if (!is.null(x = cells)) {
    reads <- reads[reads$cell %in% cells, ]
  }
  if (nrow(reads) == 0) {
    return(reads)
  }
  reads$length <- reads$end - reads$start
  return(reads)
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

# Check if path is remote
# @param x path/s to check
isRemote <- function(x) {
  return(grepl(pattern = "^http|^ftp", x = x))
}

# row merge list of matrices
# @param mat.list list of sparse matrices
# @param new.rownames rownames to assign merged matrix
#' @importFrom Seurat RowMergeSparseMatrices
MergeMatrixParts <- function(mat.list, new.rownames) {
  # RowMergeSparseMatrices only exported in Seurat release Dec-2019 (3.1.2)
  merged.all <- mat.list[[1]]
  for (i in 2:length(x = mat.list)) {
    merged.all <- RowMergeSparseMatrices(
      mat1 = merged.all,
      mat2 = mat.list[[i]]
    )
  }
  # reorder rows to match genomic ranges
  merged.all <- merged.all[new.rownames, ]
  return(merged.all)
}

# Run GetReadsInRegion for a list of Fragment objects
# concatenate the output dataframes and return
# @param object A Seurat or ChromatinAssay object
# @param region Genomic region to extract fragments for
# @param fragment.list A list of Fragment objects. If NULL, pull them from the
# object
# @param assay Name of assay to use if supplying a Seurat object
#' @importFrom Seurat DefaultAssay
#' @importFrom Rsamtools TabixFile
#' @importFrom GenomeInfoDb keepSeqlevels
MultiGetReadsInRegion <- function(
  object,
  region,
  fragment.list = NULL,
  assay = NULL,
  ...
) {
  if (inherits(x = object, what = "Seurat")) {
    # pull the assay
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    object <- object[[assay]]
  }
  fragment.list <- SetIfNull(
    x = fragment.list,
    y = Fragments(object = object)
  )
  if (length(x = fragment.list) == 0) {
    # no fragments set
    stop("No fragment files found")
  }
  res <- data.frame()
  for (i in seq_along(along.with = fragment.list)) {
    tbx.path <- GetFragmentData(object = fragment.list[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragment.list[[i]], slot = "cells")
    tabix.file <- TabixFile(file = tbx.path)
    open(con = tabix.file)
    reads <- GetReadsInRegion(
      cellmap = cellmap,
      region = region,
      tabix.file = tabix.file,
      ...
    )
    res <- rbind(res, reads)
    close(con = tabix.file)
  }
  return(res)
}

# Generate matrix of integration sites
#
# Generates a cell-by-position matrix of Tn5 integration sites
# centered on a given region (usually a DNA sequence motif). This
# matrix can be used for downstream footprinting analysis.
#
# @param cellmap A mapping of cell names in the fragment file to cell names in
# the Seurat object. Should be a named vector where each element is a cell name
# that appears in the fragment file and the name of each element is the
# name of the cell in the Seurat object.
# @param region A set of GRanges containing the regions of interest
# @param cells Which cells to include in the matrix. If NULL, use all cells in
# the cellmap
# @param tabix.file A \code{\link[Rsamtools]{TabixFile}} object.
# @param verbose Display messages
#' @importFrom Matrix sparseMatrix
#' @importFrom Rsamtools TabixFile
#' @importMethodsFrom GenomicRanges width start end
# @return Returns a sparse matrix
SingleFileCutMatrix <- function(
  cellmap,
  region,
  cells = NULL,
  tabix.file,
  verbose = TRUE
) {
  # if multiple regions supplied, must be the same width
  cells <- SetIfNull(x = cells, y = names(x = cellmap))
  if (length(x = region) == 0) {
    return(NULL)
  }
  fragments <- GetReadsInRegion(
    region = region,
    cellmap = cellmap,
    cells = cells,
    tabix.file = tabix.file,
    verbose = verbose
  )
  start.lookup <- start(x = region)
  names(start.lookup) <- seq_along(region)
  # if there are no reads in the region
  # create an empty matrix of the correct dimension
  if (nrow(x = fragments) == 0) {
    cut.matrix <- sparseMatrix(
      i = NULL,
      j = NULL,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  } else {
    fragstarts <- start.lookup[fragments$ident] + 1
    cut.df <- data.frame(
      position = c(fragments$start, fragments$end) - fragstarts,
      cell = c(fragments$cell, fragments$cell),
      stringsAsFactors = FALSE
    )
    cut.df <- cut.df[
      (cut.df$position > 0) & (cut.df$position <= width(x = region)[[1]]),
      ]
    cell.vector <- seq_along(along.with = cells)
    names(x = cell.vector) <- cells
    cell.matrix.info <- cell.vector[cut.df$cell]
    cut.matrix <- sparseMatrix(
      i = cell.matrix.info,
      j = cut.df$position,
      x = 1,
      dims = c(length(x = cells), width(x = region)[[1]])
    )
  }
  rownames(x = cut.matrix) <- cells
  colnames(x = cut.matrix) <- seq_len(width(x = region)[[1]])
  return(cut.matrix)
}

# Generate matrix of integration sites
#
# Generates a cell-by-position matrix of Tn5 integration sites.
#
# @param object A Seurat object
# @param region A GRanges object containing the region of interest
# @param assay A name of assay to use. Must be a \code{\link{ChromatinAssay}}
# containing a list of \code{\link{Fragment}} objects.
# @param cells Which cells to include in the matrix. If NULL (default), use all
# cells in the object
# @param group.by Name of grouping variable to use
# @param verbose Display messages
# @return Returns a sparse matrix
#' @importFrom Seurat DefaultAssay
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom GenomeInfoDb keepSeqlevels
CutMatrix <- function(
  object,
  region,
  group.by = NULL,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  # run SingleFileCutMatrix for each fragment file and combine results
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  fragments <- Fragments(object = object[[assay]])
  res <- list()
  for (i in seq_along(along.with = fragments)) {
    fragment.path <- GetFragmentData(object = fragments[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
    tabix.file <- TabixFile(file = fragment.path)
    open(con = tabix.file)
    # remove regions that aren't in the fragment file
    seqnames.in.both <- intersect(
      x = seqnames(x = region),
      y = seqnamesTabix(file = tabix.file)
    )
    region <- keepSeqlevels(
      x = region,
      value = seqnames.in.both,
      pruning.mode = "coarse"
    )
    if (length(x = region) != 0) {
      cm <- SingleFileCutMatrix(
        region = region,
        cellmap = cellmap,
        tabix.file = tabix.file,
        cells = cells,
        verbose = FALSE
      )
      res[[i]] <- cm
    }
    close(con = tabix.file)
  }
  res <- Reduce(f = `+`, x = res)
  return(res)
}

# Generate cut matrix for many regions
#
# Run CutMatrix on multiple regions and add them together.
# Assumes regions are pre-aligned.
#
# @param object A Seurat object
# @param regions A set of GRanges
# @param group.by Name of grouping variable to use
# @param fragments A list of Fragment objects
# @param assay Name of the assay to use
# @param cells Vector of cells to include
# @param verbose Display messages
#' @importFrom Rsamtools TabixFile seqnamesTabix
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomeInfoDb keepSeqlevels
MultiRegionCutMatrix <- function(
  object,
  regions,
  group.by = NULL,
  fragments = NULL,
  assay = NULL,
  cells = NULL,
  verbose = FALSE
) {
  if (inherits(x = object, what = "Seurat")) {
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    object <- object[[assay]]
  }
  fragments <- SetIfNull(x = fragments, y = Fragments(object = object))
  res <- list()
  if (length(x = fragments) == 0) {
    stop("No fragment files present in assay")
  }
  for (i in seq_along(along.with = fragments)) {
    frag.path <- GetFragmentData(object = fragments[[i]], slot = "path")
    cellmap <- GetFragmentData(object = fragments[[i]], slot = "cells")
    if (is.null(x = cellmap)) {
      cellmap <- colnames(x = object)
      names(x = cellmap) <- cellmap
    }
    tabix.file <- TabixFile(file = frag.path)
    open(con = tabix.file)
    # remove regions that aren't in the fragment file
    common.seqlevels <- intersect(
      x = seqlevels(x = regions),
      y = seqnamesTabix(file = tabix.file)
    )
    regions <- keepSeqlevels(
      x = regions,
      value = common.seqlevels,
      pruning.mode = "coarse"
    )
    cm <- SingleFileCutMatrix(
      cellmap = cellmap,
      tabix.file = tabix.file,
      region = regions,
      verbose = verbose
    )
    close(con = tabix.file)
    res[[i]] <- cm
  }
  # each matrix contains data for different cells at same positions
  # bind all matrices together
  res <- do.call(what = rbind, args = res)
  return(res)
}

# Create cut site pileup matrix
#
# For a set of aligned genomic ranges, find the total number of
# integration sites per cell per base.
#
# @param object A Seurat object
# @param regions A GRanges object
# @param assay Name of the assay to use
# @param cells Which cells to include. If NULL, use all cells
# @param verbose Display messages
#' @importMethodsFrom GenomicRanges strand
CreateRegionPileupMatrix <- function(
  object,
  regions,
  assay = NULL,
  cells = NULL,
  verbose = TRUE
) {
  if (length(x = regions) == 0) {
    stop("No regions supplied")
  }
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
  if (is.null(x = cut.matrix.plus)) {
    full.matrix <- cut.matrix.minus[, rev(x = colnames(x = cut.matrix.minus))]
  } else if (is.null(x = cut.matrix.minus)) {
    full.matrix <- cut.matrix.plus
  } else {
    full.matrix <- cut.matrix.plus + cut.matrix.minus[, rev(
      x = colnames(x = cut.matrix.minus)
    )]
  }
  # rename so 0 is center
  region.width <- width(x = regions)[[1]]
  midpoint <- round(x = (region.width / 2))
  colnames(full.matrix) <- seq_len(length.out = region.width) - midpoint
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
# be computed using the number of cells in the group and the average number of
# counts in the group.
# @param normalize Perform sequencing depth and cell count normalization
# @param scale.factor Scaling factor to use. If NULL (default), will use the
# median normalization factor for all the groups.
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
  all.groups <- as.character(x = unique(x = groups))
  if (any(is.na(x = groups))) {
    all.groups <- c(all.groups, NA)
  }
  ngroup <- length(x = all.groups)
  npos <- ncol(x = mat)

  group <- unlist(
    x = lapply(X = all.groups, FUN = function(x) rep(x, npos))
  )
  position <- rep(x = as.numeric(x = colnames(x = mat)), ngroup)
  count <- vector(mode = "numeric", length = npos * ngroup)

  for (i in seq_along(along.with = all.groups)) {
    grp <- all.groups[[i]]
    if (is.na(x = grp)) {
      pos.cells <- names(x = groups)[is.na(x = groups)]
    } else {
      pos.cells <- names(x = groups)[groups == all.groups[[i]]]
    }
    if (length(x = pos.cells) > 1) {
      totals <- fun(x = mat[pos.cells, ])
    } else {
      totals <- mat[pos.cells, ]
    }
    count[((i - 1) * npos + 1):((i * npos))] <- totals
  }

  # construct dataframe
  coverages <- data.frame(
    "group" = group, "position" = position, "count" = count,
    stringsAsFactors = FALSE
  )

  if (normalize) {
    scale.factor <- SetIfNull(
      x = scale.factor, y = median(x = group.scale.factors)
    )
    coverages$norm.value <- coverages$count /
      group.scale.factors[coverages$group] * scale.factor
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
# @param record.ident Add a column recording which region the reads overlapped
# with
#' @importFrom stringi stri_split_fixed
#' @importFrom S4Vectors elementNROWS
TabixOutputToDataFrame <- function(reads, record.ident = TRUE) {
  if (record.ident) {
    nrep <- elementNROWS(x = reads)
  }
  reads <- unlist(x = reads, use.names = FALSE)
  if (length(x = reads) == 0) {
    df <- data.frame(
      "chr" = "",
      "start" = "",
      "end" = "",
      "cell" = "",
      "count" = ""
    )
    df <- df[-1, ]
    return(df)
  }
  reads <- stri_split_fixed(str = reads, pattern = "\t")
  n <- length(x = reads[[1]])
  unlisted <- unlist(x = reads)
  e1 <- unlisted[n * (seq_along(along.with = reads)) - (n - 1)]
  e2 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 2)])
  e3 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 3)])
  e4 <- unlisted[n * (seq_along(along.with = reads)) - (n - 4)]
  e5 <- as.numeric(x = unlisted[n * (seq_along(along.with = reads)) - (n - 5)])
  df <- data.frame(
    "chr" = e1,
    "start" = e2,
    "end" = e3,
    "cell" = e4,
    "count" = e5,
    stringsAsFactors = FALSE,
    check.rows = FALSE,
    check.names = FALSE
  )
  if (record.ident) {
    df$ident <- rep(x = seq_along(along.with = nrep), nrep)
  }
  return(df)
}

# Extract delimiter information from a string.
#
# From Seurat
#
# Parses a string (usually a cell name) and extracts fields based on a delimiter
#
# @param string String to parse.
# @param field Integer(s) indicating which field(s) to extract. Can be a vector
# multiple numbers.
# @param delim Delimiter to use, set to underscore by default.
#
# @return A new string, that parses out the requested fields, and (if multiple),
# rejoins them with the same delimiter
#' @importFrom stringi stri_split_fixed
#
ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(
    x = unlist(x = stri_split_fixed(
      str = as.character(x = field), pattern = ",")
    )
  )
  if (length(x = fields) == 1) {
    return(stri_split_fixed(str = string, pattern = delim)[[1]][field])
  }
  return(paste(
    stri_split_fixed(str = string, pattern = delim)[[1]][fields],
    collapse = delim))
}

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

# Find matrix indices corresponding to overlapping genomic ranges
# @param assay.list A list of ChromatinAssay objects
# @param all.ranges Combined genomic ranges for all objects. This should be the
# set of ranges that \code{reduce} was run on to get \code{reduced.ranges}
# @param reduced.ranges A set of reduced genomic ranges containing the rev.map
# information
GetRowsToMerge <- function(assay.list, all.ranges, reduced.ranges) {
  revmap <- as.vector(x = reduced.ranges$revmap)

  # get indices of ranges that changed
  revmap.lengths <- sapply(X = revmap, FUN = length)
  changed.ranges <- which(x = revmap.lengths > 1)
  grange.string <- GRangesToString(grange = reduced.ranges[changed.ranges])

  # preallocate
  offsets <- list()
  results <- list()
  matrix.indices <- vector(
    mode = "numeric",
    length = length(x = changed.ranges) * 2
  )
  granges <- vector(
    mode = "character",
    length = length(x = changed.ranges) * 2
  )
  for (i in seq_along(along.with = assay.list)) {
    indices <- which(x = all.ranges$dataset == i)
    offsets[[i]] <- min(indices) - 1
    offsets[[i]][[2]] <- max(indices) + 1
    results[['matrix']][[i]] <- matrix.indices
    results[['grange']][[i]] <- granges
  }

  # find sets of ranges for each dataset
  counter <- vector(mode = "numeric", length = length(x = assay.list))
  for (x in seq_along(along.with = changed.ranges)) {
    idx <- changed.ranges[[x]]
    all.assay <- revmap[[idx]]
    for (i in seq_along(along.with = assay.list)) {
      this.assay <- all.assay[
        (all.assay > offsets[[i]][1]) & (all.assay < offsets[[i]][2])
      ]
      mat.idx <- this.assay - offsets[[i]][1]
      mat.idx <- mat.idx[mat.idx < offsets[[i]][2] & mat.idx > 0]
      for (y in seq_along(along.with = mat.idx)) {
        counter[i] <- counter[i] + 1
        results[['matrix']][[i]][[counter[i]]] <- mat.idx[[y]]
        results[['grange']][[i]][[counter[i]]] <- grange.string[[x]]
      }
    }
  }
  # remove trailing extra values in each vector
  for (i in seq_along(along.with = assay.list)) {
    results$matrix[[i]] <- results$matrix[[i]][1:counter[i]]
    results$grange[[i]] <- results$grange[[i]][1:counter[i]]
  }
  return(results)
}

# Merge rows of count matrices with overlapping genomic ranges
# @param mergeinfo The output of GetRowsToMerge: a list of matrix indices
#  and matrix rownames to be merged for each assay
# @param assay.list List of assays
# @param verbose Display messages
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom Matrix rowSums
#' @importMethodsFrom Matrix t
MergeOverlappingRows <- function(
  mergeinfo,
  assay.list,
  slot = "counts",
  verbose = TRUE
) {
  merge.counts <- list()
  for (i in seq_along(along.with = assay.list)) {
    # get count matrix
    counts <- GetAssayData(object = assay.list[[i]], slot = slot)

    if (nrow(x = counts) == 0) {
      # no counts, only data
      # skip row merge and return empty counts matrices
      merge.counts <- lapply(
        X = seq_along(along.with = assay.list),
        FUN = matrix,
        nrow = 0,
        ncol = 0
      )
      return(merge.counts)
    }

    # transpose for faster access since matrix is column major
    counts <- t(x = counts)

    # get rows to merge
    mrows <- mergeinfo$matrix[[i]]
    new.rownames <- mergeinfo$grange[[i]]
    nrep <- rle(x = new.rownames)

    # allocate
    todelete <- c()
    newmat <- vector(
      mode = "list",
      length = length(new.rownames)
    )
    newmat.names <- vector(
      mode = "character",
      length = length(x = new.rownames)
    )
    x <- 1  # row index for matrix
    y <- 1  # counter for list index
    if (verbose & length(x = nrep$lengths) > 1) {
      pb <- txtProgressBar(
        min = 1,
        max = length(x = nrep$lengths),
        style = 3,
        file = stderr()
      )
    }
    to.rename.idx <- vector(
      mode = "numeric", length = length(x = nrep$lengths)
    )
    to.rename.names <- vector(
      mode = "character", length = length(x = nrep$lengths)
    )
    idx.counter <- 0
    for (j in seq_along(along.with = nrep$lengths)) {
      rowrun <- nrep$lengths[[j]]
      new.feature.name <- nrep$values[[j]]
      index.range <- x:(x + rowrun - 1)
      matrix.index <- mrows[index.range]
      if (rowrun < 2) {
        idx.counter <- idx.counter + 1
        # no merge needed, just rename row in-place
        # store row indices and names to do the change in one step at the end
        to.rename.idx[idx.counter] <- matrix.index
        to.rename.names[idx.counter] <- new.feature.name
      } else {
        # merge multiple rows and add to list
        newmat[[y]] <- rowSums(x = counts[, matrix.index])
        # mark merged row for deletion
        todelete <- c(todelete, matrix.index)
        # add row names
        newmat.names[y] <- new.feature.name
        y <- y + 1
      }
      if (verbose & length(x = nrep$lengths) > 1) {
        setTxtProgressBar(pb = pb, value = j)
      }
      x <- x + rowrun
    }
    # remove extra elements in vectors
    to.rename.idx <- to.rename.idx[1:idx.counter]
    to.rename.names <- to.rename.names[1:idx.counter]
    newmat <- newmat[1:(y - 1)]
    newmat.names <- newmat.names[1:(y - 1)]

    # transpose back
    counts <- t(x = counts)

    # rename matrix rows that weren't merged
    rownames(counts)[to.rename.idx] <- to.rename.names

    if (y == 1) {
      # no rows were merged, can return counts
      merge.counts[[i]] <- counts
    } else if (y == 2) {
      # only one element
      tomerge <- matrix(data = newmat[[1]], nrow = 1)
      colnames(x = tomerge) <- names(x = newmat[[1]])
      rownames(x = tomerge) <- newmat.names
      tomerge <- tomerge[, colnames(x = counts)]
      counts <- rbind(counts, tomerge)
      merge.counts[[i]] <- counts
    } else {
      # construct sparse matrix
      if (verbose) {
        message("\nBinding matrix rows")
      }
      merged.mat <- Reduce(f = rbind, x = newmat)
      rownames(merged.mat) <- newmat.names
      merged.mat <- as(object = merged.mat, Class = "dgCMatrix")

      # remove rows from count matrix that were merged
      mat.rows <- seq_len(length.out = nrow(x = counts))
      tokeep <- setdiff(mat.rows, todelete)
      counts <- counts[tokeep, ]

      # add new merged rows to counts
      counts <- rbind(counts, merged.mat)
      merge.counts[[i]] <- counts
    }
  }
  return(merge.counts)
}

#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors elementNROWS
PartialMatrix <- function(tabix, regions, sep = c("-", "-"), cells = NULL) {
  # construct sparse matrix for one set of regions
  # names of the cells vector can be ignored here, conversion is handled in
  # the parent functions
  open(con = tabix)
  cells.in.regions <- GetCellsInRegion(
    tabix = tabix,
    region = regions,
    cells = cells
  )
  close(con = tabix)
  gc(verbose = FALSE)
  nrep <- elementNROWS(x = cells.in.regions)
  if (all(nrep == 0) & !is.null(x = cells)) {
    # no fragments
    # zero for all requested cells
    featmat <- sparseMatrix(
      dims = c(length(x = regions), length(x = cells)),
      i = NULL,
      j = NULL
    )
    rownames(x = featmat) <- GRangesToString(grange = regions)
    colnames(x = featmat) <- cells
    featmat <- as(object = featmat, Class = "dgCMatrix")
    return(featmat)
  } else if (all(nrep == 0)) {
    # no fragments, no cells requested
    # create empty matrix
    featmat <- sparseMatrix(
      dims = c(length(x = regions), 0),
      i = NULL,
      j = NULL
    )
    rownames(x = featmat) <- GRangesToString(grange = regions)
    featmat <- as(object = featmat, Class = "dgCMatrix")
    return(featmat)
  } else {
    # fragments detected
    if (is.null(x = cells)) {
      all.cells <- unique(x = unlist(x = cells.in.regions))
      cell.lookup <- seq_along(along.with = all.cells)
      names(x = cell.lookup) <- all.cells
    } else {
      cell.lookup <- seq_along(along.with = cells)
      names(cell.lookup) <- cells
    }
    # convert cell name to integer
    cells.in.regions <- unlist(x = cells.in.regions)
    cells.in.regions <- unname(obj = cell.lookup[cells.in.regions])
    all.features <- GRangesToString(grange = regions, sep = sep)
    feature.vec <- rep(x = seq_along(along.with = all.features), nrep)
    featmat <- sparseMatrix(
      i = feature.vec,
      j = cells.in.regions,
      x = rep(x = 1, length(x = cells.in.regions))
    )
    featmat <- as(Class = "dgCMatrix", object = featmat)
    rownames(x = featmat) <- all.features[1:max(feature.vec)]
    colnames(x = featmat) <- names(x = cell.lookup)[1:max(cells.in.regions)]
    # add zero columns for missing cells
    if (!is.null(x = cells)) {
      featmat <- AddMissingCells(x = featmat, cells = cells)
    }
    # add zero rows for missing features
    missing.features <- all.features[!(all.features %in% rownames(x = featmat))]
    if (length(x = missing.features) > 0) {
      null.mat <- sparseMatrix(
        i = c(),
        j = c(),
        dims = c(length(x = missing.features), ncol(x = featmat))
      )
      rownames(x = null.mat) <- missing.features
      null.mat <- as(object = null.mat, Class = "dgCMatrix")
      featmat <- rbind(featmat, null.mat)
    }
    return(featmat)
  }
}

# Convert PFMMatrix to
# @param x A PFMatrix
PFMatrixToList <- function(x) {
  if (!requireNamespace("TFBSTools", quietly = TRUE)) {
    stop("Please install TFBSTools.
         https://www.bioconductor.org/packages/TFBSTools/")
  }
  position.matrix <- TFBSTools::Matrix(x = x)
  name.use <- TFBSTools::name(x = x)
  return(list("matrix" = position.matrix, "name" = name.use))
}

#' @importFrom Matrix rowMeans rowSums
SparseRowVar <- function(x) {
  return(rowSums(x = (x - rowMeans(x = x)) ^ 2) / (dim(x = x)[2] - 1))
}

#' @importMethodsFrom Matrix t
SparseColVar <- function(x) {
  return(SparseRowVar(x = t(x = x)))
}

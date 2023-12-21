#' @include generics.R
#'
NULL

#' Compute TSS enrichment score per cell
#'
#' Compute the transcription start site (TSS) enrichment score for each cell,
#' as defined by ENCODE:
#' \url{https://www.encodeproject.org/data-standards/terms/}.
#'
#' The computed score will be added to the object metadata as "TSS.enrichment".
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param tss.positions A GRanges object containing the TSS positions. If NULL,
#' use the genomic annotations stored in the assay.
#' @param n Number of TSS positions to use. This will select the first _n_
#' TSSs from the set. If NULL, use all TSSs (slower).
#' @param cells A vector of cells to include. If NULL (default), use all cells
#' in the object
#' @param fast Just compute the TSS enrichment score, without storing the
#' base-resolution matrix of integration counts at each site. This reduces the
#' memory required to store the object but does not allow plotting the
#' accessibility profile at the TSS.
#' @param process_n Number of regions to process at a time if using \code{fast}
#' option.
#' @param verbose Display messages
#' @param region_extension Distance extended upstream and downstream from TSS
#' in which to calculate enrichment and background.
#'
#' @importFrom Matrix rowMeans
#' @importFrom methods slot
#' @importFrom stats ecdf
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges start width strand
#' @importFrom SeuratObject DefaultAssay
#'
#' @return Returns a \code{\link[SeuratObject]{Seurat}} object
#' @export
#' @concept qc
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
#' Fragments(atac_small) <- CreateFragmentObject(
#'   path = fpath,
#'   cells = colnames(atac_small),
#'   tolerance = 0.5
#' )
#' TSSEnrichment(object = atac_small)
#' }
TSSEnrichment <- function(
  object,
  tss.positions = NULL,
  n = NULL,
  fast = TRUE,
  assay = NULL,
  cells = NULL,
  process_n = 2000,
  verbose = TRUE,
  region_extension = 1000
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  # first check that fragments are present
  frags <- Fragments(object = object[[assay]])
  if (length(x = frags) == 0) {
    stop("No fragment files present in assay")
  }
  if (is.null(x = tss.positions)) {
    if (verbose) {
      message("Extracting TSS positions")
    }
    # work out TSS positions from gene annotations
    annotations <- Annotation(object = object[[assay]])
    if (is.null(x = annotations)) {
      stop("No gene annotations present in assay")
    }
    tss.positions <- GetTSSPositions(ranges = annotations)
  }
  if (!is.null(x = n)) {
    if (n > length(x = tss.positions)) {
      n <- length(x = tss.positions)
    }
    tss.positions <- tss.positions[1:n, ]
  }
  
  # exclude chrM
  sn <- seqnames(x = tss.positions)
  tss.positions <- tss.positions[!as.character(sn) %in% c("chrM", "Mt", "MT")]
  
  if (fast) {
    # just compute the TSS enrichment score without storing the full matrix
    object <- TSSFast(
      object = object,
      assay = assay,
      tss.positions = tss.positions,
      process_n = process_n,
      verbose = verbose,
      region_extension = region_extension
    )
    return(object)
  }
  
  tss.positions <- Extend(
    x = tss.positions,
    upstream = region_extension,
    downstream = region_extension,
    from.midpoint = TRUE
  )
  cutmatrix <- CreateRegionPileupMatrix(
    object = object,
    regions = tss.positions,
    assay = assay,
    cells = cells,
    verbose = verbose
  )
  
  # compute mean read counts in 100 bp at each flank for each cell
  # (200 bp total averaged)
  if (verbose) {
    message("Computing mean insertion frequency in flanking regions")
  }
  total_region_length <- (2 * region_extension) + 1
  right_flank <- seq.int(from = (total_region_length - 99), to = total_region_length)
  flanking.mean <- rowMeans(x = cutmatrix[, c(1:100, right_flank)])
  
  # if the flanking mean is 0 for any cells, the enrichment score will be zero.
  # instead replace with the mean from the whole population
  flanking.mean[is.na(x = flanking.mean)] <- 0
  flanking.mean[flanking.mean == 0] <- mean(flanking.mean, na.rm = TRUE)
  
  # compute fold change at each position relative to flanking mean
  # (flanks should start at 1)
  if (verbose) {
    message("Normalizing TSS score")
  }
  
  norm.matrix <- cutmatrix / flanking.mean
  
  # Take signal value at center of distribution after normalization as
  # TSS enrichment score, average the 1001 bases at the center
  center_region <- seq.int(from = (region_extension - 500), to = (region_extension + 500))
  object$TSS.enrichment <- rowMeans(x = norm.matrix[, center_region], na.rm = TRUE)
  e.dist <- ecdf(x = object$TSS.enrichment)
  object$TSS.percentile <- round(
    x = e.dist(object$TSS.enrichment),
    digits = 2
  )
  
  # store expected as one additional row in the matrix
  expected.insertions <- rep(1, ncol(x = cutmatrix))
  expected.insertions <- t(x = as.matrix(x = expected.insertions))
  rownames(x = expected.insertions) <- "expected"
  
  # encode motif position as additional row in matrix
  motif.vec <- t(x = matrix(
    data = c(
      rep(x = 0, region_extension),
      1,
      rep(x = 0, region_extension)
    )
  )
  )
  rownames(x = motif.vec) <- "motif"
  
  # append
  norm.matrix <- rbind(norm.matrix, expected.insertions)
  norm.matrix <- rbind(norm.matrix, motif.vec)
  
  # store the normalized TSS matrix
  object <- suppressWarnings(SetAssayData(
    object = object,
    assay = assay,
    slot = "positionEnrichment",
    new.data = norm.matrix,
    key = "TSS"
  ))
  return(object)
}

#' @importFrom Rsamtools TabixFile
#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @importFrom stats ecdf
#' @importFrom Matrix rowSums
#' @importFrom SeuratObject DefaultAssay
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pblapply
TSSFast <- function(
  object,
  tss.positions,
  assay = NULL,
  process_n = 2000,
  verbose = TRUE,
  region_extension = 1000
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  
  # extract fragments
  frags <- Fragments(object = object[[assay]])
  if (length(x = frags) == 0) {
    stop("Fragments file not set for assay ", assay)
  }
  
  # get regions
  upstream.flank <- Extend(
    x = tss.positions,
    upstream = region_extension,
    downstream = -1 * (region_extension - 99),
    from.midpoint = TRUE
  )
  downstream.flank <- Extend(
    x = tss.positions,
    upstream = -1 * (region_extension - 99),
    downstream = region_extension,
    from.midpoint = TRUE
  )
  centers <- Extend(
    x = tss.positions,
    upstream = 500,
    downstream = 500,
    from.midpoint = TRUE
  )
  
  # chunk ranges
  process_n <- SetIfNull(x = process_n, y = length(x = centers))
  nchunk <- ceiling(x = length(x = upstream.flank) / process_n)
  upstream.flank <- ChunkGRanges(
    granges = upstream.flank,
    nchunk = nchunk
  )
  downstream.flank <- ChunkGRanges(
    granges = downstream.flank,
    nchunk = nchunk
  )
  centers <- ChunkGRanges(
    granges = centers,
    nchunk = nchunk
  )
  
  # iterate over fragment files and parts of region
  if (verbose) {
    message("Extracting fragments at TSSs")
  }
  
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  
  center.counts <- vector(mode = "numeric", length = ncol(x = object))
  flank.counts <- vector(mode = "numeric", length = ncol(x = object))
  for (i in seq_along(along.with = frags)) {
    # open fragment file
    tbx.path <- GetFragmentData(object = frags[[i]], slot = "path")
    cellmap <- GetFragmentData(object = frags[[i]], slot = "cells")
    if (is.null(x = cellmap)) {
      cellmap <- colnames(x = object)
      names(x = cellmap) <- cellmap
    } else {
      cellmap <- cellmap[intersect(names(x = cellmap), colnames(x = object))]
    }
    tbx <- TabixFile(
      file = tbx.path,
      index = GetIndexFile(fragment = tbx.path, verbose = FALSE)
    )
    # iterate over chunked ranges
    res <- mylapply(
      X = seq_along(along.with = centers),
      FUN = function(x) {
        extract_tss_counts(
          cellnames = colnames(x = object),
          region.centers = centers[[x]],
          upstream = upstream.flank[[x]],
          downstream = downstream.flank[[x]],
          tabix.file = tbx,
          cell.name.map = cellmap
        )
      }
    )
    
    # sum results from each chunk of granges
    cc <- lapply(X = res, FUN = `[[`, 1)
    fc <- lapply(X = res, FUN = `[[`, 2)
    center.counts <- center.counts + Reduce(f = `+`, x = cc)
    flank.counts <- flank.counts + Reduce(f = `+`, x = fc)
  }
  
  if (verbose) {
    message("\nComputing TSS enrichment score")
  }
  
  # take mean accessibility per base
  flank.mean <- flank.counts / 200
  flank.mean[flank.counts == 0] <- mean(x = flank.mean, na.rm = TRUE)
  
  center.norm <- center.counts / flank.mean
  
  # replace NA with 0
  center.norm[is.na(center.norm)] <- 0
  
  # compute TSS enrichment score and add to object
  object$TSS.enrichment <- center.norm / 1001
  e.dist <- ecdf(x = object$TSS.enrichment)
  object$TSS.percentile <- round(
    x = e.dist(object$TSS.enrichment),
    digits = 2
  )
  return(object)
}

## Not exported ##

#' @importFrom GenomeInfoDb seqlevels keepSeqlevels
#' @importFrom Rsamtools seqnamesTabix
#' @importFrom Matrix rowSums
extract_tss_counts <- function(
  cellnames,
  region.centers,
  tabix.file,
  upstream,
  downstream,
  cell.name.map
) {
  tabix.file <- open(con = tabix.file)
  # initialize vectors
  fc <- vector(mode = "numeric", length = length(x = cellnames))
  names(x = fc) <- cellnames
  cc <- vector(mode = "numeric", length = length(x = cellnames))
  names(x = cc) <- cellnames
  
  # remove seqlevels not present in fragment file
  common.seqlevels <- intersect(
    x = seqlevels(x = region.centers),
    y = seqnamesTabix(file = tabix.file)
  )
  if (length(x = common.seqlevels) == 0) {
    close(con = tabix.file)
    return(list(cc, fc))
  }
  uflanks.use <- keepSeqlevels(
    x = upstream,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  dflanks.use <- keepSeqlevels(
    x = downstream,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  centers.use <- keepSeqlevels(
    x = region.centers,
    value = common.seqlevels,
    pruning.mode = "coarse"
  )
  
  # count integration events
  cuts.center <- SingleFileCutMatrix(
    cellmap = cell.name.map,
    tabix.file = tabix.file,
    region = centers.use,
    verbose = FALSE
  )
  counts.center <- rowSums(x = cuts.center)
  cuts.flank <- SingleFileCutMatrix(
    cellmap = cell.name.map,
    tabix.file = tabix.file,
    region = c(uflanks.use, dflanks.use),
    verbose = FALSE
  )
  counts.flank <- rowSums(x = cuts.flank)
  cc[names(x = counts.center)] <- cc + as.vector(x = counts.center)
  fc[names(x = counts.flank)] <- fc + as.vector(x = counts.flank)
  close(con = tabix.file)
  return(list(cc, fc))
}

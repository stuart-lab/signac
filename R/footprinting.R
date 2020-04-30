#' @include generics.R
#'
NULL


#' @param regions A set of genomic ranges containing the motif instances
#' @param genome A \code{\link[BSgenome]{BSgenome}} object
#' @param motif.name Name of a motif stored in the assay to footprint. If not
#' supplied, must supply a set of regions.
#' @param group.by Grouping variable for the cells
#' @param idents Which identities to include
#' @param upstream Number of bases to extend upstream
#' @param downstream Number of bases to extend downstream
#' @param verbose Display messages
#' @param ... Arguments passed to other methods
#' @importFrom Biostrings getSeq
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
#' @export
#' @rdname Footprint
#' @method Footprint ChromatinAssay
Footprint.ChromatinAssay <- function(
  object,
  genome,
  motif.name = NULL,
  regions = NULL,
  group.by = NULL,
  idents = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
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
  }
  regions.use <- Extend(
    x = regions,
    upstream = upstream + 3,
    downstream = downstream + 3
  )
  dna.sequence <- getSeq(x = genome, regions.use)
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
  if (verbose) {
    message("Computing base composition at motif sites")
  }
  dna.string <- as.character(dna.sequence)
  row.index <- c()
  total.bases <- upstream + downstream + 1
  # TODO add progress bar here, add future. think about other ways to do this.
  for (i in 1:length(x = dna.string)) {
    for (j in 1:total.bases) {
      row.index <- c(row.index, substring(text = dna.string[[i]], first = j, last = j + 5))
    }
  }
  unique.hexamer <- unique(x = row.index)
  hexamer.row.index <- match(x = row.index, table = unique.hexamer)
  hexamer.col.index <- rep(1:total.bases, length(x = dna.string))
  hexamer.matrix <- sparseMatrix(
    i = hexamer.row.index,
    j = hexamer.col.index,
    x = 1
  )
  rownames(hexamer.matrix) <- unique.hexamer
  colnames(hexamer.matrix) <- 1:total.bases
  if (verbose) {
    message("Computing expected Tn5 insertions per base")
  }
  hexamer.matrix <- hexamer.matrix[names(x = bias), ]
  expected.insertions <- crossprod(x = hexamer.matrix, y = as.matrix(x = bias))
  if (verbose) {
    message("Computing observed Tn5 insertions per base")
  }
  insertion.matrix <- CreateRegionPileupMatrix(
    object = object,
    regions = regions,
    upstream = upstream,
    downstream = downstream
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  group.counts <- ApplyMatrixByGroup(
    mat = insertion.matrix,
    groups = obj.groups,
    fun = colSums,
    normalize = FALSE
  )
  if (verbose) {
    message("Computing observed/expected Tn5 insertions per base")
  }
  flanks <- c(1:50, (total.bases - 50):total.bases)
  norm.factor.expected <- mean(expected.insertions[flanks, ])
  norm.expected <- expected.insertions / norm.factor.expected
  unique.groups <- unique(x = group.counts$group)
  norm.counts <- data.frame()
  for (i in seq_along(along.with = unique.groups)) {
    group.use <- group.counts[group.counts$group == unique.groups[[i]], 'count']
    norm.factor <- mean(x = group.use[flanks])
    normalized.group.counts <- group.use / norm.factor
    obs.expect <- normalized.group.counts / as.vector(x = norm.expected)
    norm.counts <- rbind(norm.counts, data.frame(
      postion = -upstream:downstream + 1,
      observed.over.expected = obs.expect,
      group = as.character(x = unique.groups[[i]])
    ))
    # norm.counts[[as.character(unique.groups[[i]])]] <- obs.expect
  }
  # norm.counts[['expected']] <- as.vector(x = norm.expected)
  return(norm.counts)
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


#' @include generics.R
#'
NULL

#' Transcription factor footprinting analysis
#'
#' Compute the normalized observed/expected Tn5 insertion frequency
#' for each position surrounding a set of motif instances.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param regions A set of GRanges containing the motif instances
#' @param genome A BSgenome object
#' @param group.by Grouping variable for the cells
#' @param idents Which identities to include
#' @param upstream Number of bases to extend upstream
#' @param downstream Number of bases to extend downstream
#' @param verbose Display messages
#' @importFrom Seurat Misc DefaultAssay
#' @importFrom Biostrings getSeq
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix colSums
Footprint <- function(
  object,
  regions,
  genome,
  group.by = NULL,
  idents = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  regions.use <- Extend(x = regions, upstream = upstream + 3, downstream = downstream + 3)
  dna.sequence <- getSeq(x = genome, regions.use)
  misc.slot <- Misc(object = object[[assay]])
  if (!("Tn5.bias" %in% names(x = misc.slot))) {
    object <- InsertionBias(
      object = object,
      genome = genome
    )
    misc.slot <- Misc(object = object[[assay]])
  }
  bias <- misc.slot$Tn5.bias
  if (verbose) {
    message("Computing base composition at motif sites")
  }
  dna.string <- as.character(dna.sequence)
  row.index <- c()
  total.bases <- upstream + downstream + 1
  for (i in 1:length(x = dna.string)) {
    for (j in 1:total.bases) {
      row.index <- c(row.index, substring(text = dna.string[[i]], first = j, last = j + 5))
    }
  }
  unique.hexamer <- unique(x = row.index)
  hexamer.row.index <- match(x = row.index, table = unique.hexamer)
  hexamer.col.index <- rep(1:total.bases, length(mef.peaks))
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
  flanks <- c(1:50, (total.bases-50):total.bases)
  norm.factor.expected <- mean(expected.insertions[flanks,])
  norm.expected <- expected.insertions / norm.factor.expected
  unique.groups <- unique(x = group.counts$group)
  norm.counts <- data.frame()
  for (i in seq_along(along.with = unique.groups)) {
    group.use <- group.counts[group.counts$group == unique.groups[[i]], 'count']
    norm.factor <- mean(x = group.use[flanks])
    normalized.group.counts <- group.use / norm.factor
    obs.expect <- normalized.group.counts / as.vector(x = norm.expected)
    data.frame <- rbind(norm.counts, data.frame(
      postion = -upstream:downstream+1,
      observed.over.expected = obs.expect,
      group = as.character(x = unique.groups[[i]])
    ))
    # norm.counts[[as.character(unique.groups[[i]])]] <- obs.expect
  }
  # norm.counts[['expected']] <- as.vector(x = norm.expected)
  return(norm.counts)
}

#' @include generics.R
#'
NULL

#' @param regions A [GenomicRanges::GRanges()] object containing the
#' set of genomic ranges to quantify
#' @param assay Name of assay to use. If NULL, use the default assay
#' @param group.by Grouping variable to use when aggregating data across cells.
#' If NULL, use the active cell identities
#' @param idents Cell identities to include. If NULL, include all identities
#' @param upstream Number of bases to extend regions upstream
#' @param downstream Number of bases to extend regions downstream
#' @param verbose Display messages
#' @concept heatmap
#' @method RegionMatrix Seurat
#' @export
#' @rdname RegionMatrix
#' @importFrom SeuratObject DefaultAssay
RegionMatrix.Seurat <- function(
  object,
  regions,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  upstream = 3000,
  downstream = 3000,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (!inherits(x = object[[assay]], what = "ChromatinAssay5")) {
    stop("Requested assay is not a ChromatinAssay5 or GRangesAssay")
  }
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  mat <- RegionMatrix(
    object = object[[assay]],
    regions = regions,
    group.by = obj.groups,
    upstream = upstream,
    downstream = downstream,
    verbose = verbose,
    ...
  )
  return(mat)
}

#' @method RegionMatrix ChromatinAssay5
#' @export
#' @rdname RegionMatrix
#' @concept heatmap
RegionMatrix.ChromatinAssay5 <- function(
  object,
  regions,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  upstream = 3000,
  downstream = 3000,
  verbose = TRUE,
  ...
) {
  # check that group.by is expected format
  if (!all(names(x = group.by) %in% colnames(x = object))) {
    stop("Incorrect cell grouping information supplied")
  }
  # get relevant fragment files
  allfrag <- Fragments(object = object)
  if (length(x = allfrag) == 0) {
    stop("No fragment files present in assay")
  }
  frag.keep <- sapply(X = allfrag, FUN = function(x) {
    any(names(x = group.by) %in% Cells(x = x))
  })
  if (sum(frag.keep) == 0) {
    stop("No fragments files for requested cells")
  }
  allfrag <- allfrag[frag.keep]
  matlist <- RegionMatrix(
    object = allfrag,
    regions = regions,
    group.by = group.by,
    upstream = upstream,
    downstream = downstream,
    verbose = verbose,
    ...
  )

  return(matlist)
}

#' @method RegionMatrix default
#' @export
#' @concept heatmap
#' @rdname RegionMatrix
#' @importFrom GenomicRanges resize strand
#' @importFrom Rsamtools TabixFile scanTabix seqnamesTabix
#' @importFrom Seqinfo seqnames
#' @importFrom fastmatch fmatch
RegionMatrix.default <- function(
  object,
  regions,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  upstream = 3000,
  downstream = 3000,
  verbose = TRUE,
  ...
) {
  # object is a list of fragment objects
  if (!inherits(x = object, what = "list")) {
    object <- list(object)
  }
  all.valid <- sapply(X = object, FUN = inherits, what = "Fragment2")
  if (!all(all.valid)) {
    stop("Must supply a list of Fragment objects")
  }

  # get center of region and extend
  regions <- resize(x = regions, width = 1, fix = "center")
  regions <- suppressWarnings(expr = Extend(
    x = regions,
    upstream = upstream,
    downstream = downstream,
    from.midpoint = TRUE
  ))

  total_bases <- upstream + downstream

  # separate matrix for each group of cells
  matlist <- list()
  unique.groups <- as.character(x = unique(x = group.by))

  for (i in seq_along(along.with = object)) {
    tmplist <- list()

    # open tabix connection
    fragfile <- GetFragmentData(object = object[[i]], slot = "file.path")
    fragindex <- GetFragmentData(object = object[[i]], slot = "file.index")
    cellnames <- GetFragmentData(object = object[[i]], slot = "cells")
    tabix.file <- TabixFile(file = fragfile, index = fragindex)
    open(con = tabix.file)

    # initialize empty matrix for each group of cells
    # each row is region
    # each column is a base in region
    # (TODO: implement window sum to reduce size of matrix)
    ncol.mat <- upstream + downstream + 1
    for (j in unique.groups) {
      # TODO make this dgCMatrix instead
      tmplist[[j]] <- matrix(
        data = 0,
        nrow = length(x = regions),
        ncol = ncol.mat
      )
    }

    # remove regions on sequences that aren't in this fragment file, otherwise
    # scanTabix aborts. orig.idx maps back to the input regions so that each
    # matrix row still corresponds to the region supplied by the caller
    in.file <- as.character(x = seqnames(x = regions)) %in%
      seqnamesTabix(file = tabix.file)
    regions.use <- regions[in.file]
    orig.idx <- which(x = in.file)
    n.dropped <- sum(!in.file)
    if (n.dropped > 0) {
      warning(
        n.dropped,
        ifelse(test = n.dropped == 1, yes = " region is", no = " regions are"),
        " on seqnames not present in the fragment file ", fragfile,
        ". These will be counted as zero.",
        call. = FALSE
      )
    }

    if (length(x = regions.use) > 0) {
      frags <- scanTabix(file = tabix.file, param = regions.use)
      res <- TabixOutputToDataFrame(
        reads = frags, record.ident = TRUE
      )

      # assign counts to cell groups
      for (j in unique(x = res$ident)) { # for each region
        res_region <- res[res$ident == j, ]
        # subtract start from fragment position
        on_plus <- as.logical(strand(x = regions.use[j]) == "+" |
                                strand(x = regions.use[j]) == "*")
        res_region$start <- res_region$start - start(x = regions.use[j])
        res_region$end <- res_region$end - start(x = regions.use[j])

        # positions outside the window are ignored by tabulate() below, so
        # they are not filtered here. Dropping whole fragments that had one
        # end outside the window discarded the in-window insertion too
        for (cell in unique.groups) {
          cells.keep <- names(x = group.by[group.by == cell])
          subfrag <- res_region[
            fmatch(
              x = res_region$cell,
              table = cellnames[cells.keep],
              nomatch = 0L
            ) > 0, ,
            drop = FALSE
          ]
          if (nrow(x = subfrag) > 0) {
            if (on_plus) {
              startpos <- subfrag$start
              endpos <- subfrag$end
            } else {
              # reverse for minus strand
              startpos <- total_bases - subfrag$start
              endpos <- total_bases - subfrag$end
            }
            # accumulate insertions per position. tabulate() ignores indices
            # <= 0 and > nbins, clipping positions outside the window
            tmplist[[cell]][orig.idx[j], ] <- tmplist[[cell]][orig.idx[j], ] +
              tabulate(bin = startpos, nbins = ncol.mat) +
              tabulate(bin = endpos, nbins = ncol.mat)
          }
        }
      }
    }

    if (i == 1) {
      # one fragment file
      matlist <- tmplist
    } else {
      # sum across fragment files
      for (cell in unique.groups) {
        matlist[[cell]] <- matlist[[cell]] + tmplist[[cell]]
      }
    }
  }
  # get normalization factors
  cells.per.group <- table(group.by, useNA = "always")
  lut <- as.vector(x = cells.per.group)
  names(x = lut) <- names(x = cells.per.group)
  
  
  # store upstream and downstream parameters
  params <- list(
    "upstream" = upstream,
    "downstream" = downstream,
    "cells" = lut
  )
  results <- list("matrix" = matlist, "parameters" = params)
  return(results)
}

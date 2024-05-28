#' @include generics.R
#'
NULL

#' @param regions A \code{\link[GenomicRanges]{GRanges}} object containing the
#' set of genomic ranges to quantify
#' @param key Name to store resulting matrices under
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
  key,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  upstream = 3000,
  downstream = 3000,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("Requested assay is not a ChromatinAssay")
  }
  if (missing(x = key)) {
    stop("No key supplied")
  }
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  object[[assay]] <- RegionMatrix(
    object = object[[assay]],
    regions = regions,
    group.by = obj.groups,
    key = key,
    upstream = upstream,
    downstream = downstream,
    verbose = verbose,
    ...
  )
  return(object)
}

#' @method RegionMatrix ChromatinAssay
#' @export
#' @importFrom SeuratObject GetAssayData
#' @rdname RegionMatrix
#' @concept heatmap
RegionMatrix.ChromatinAssay <- function(
  object,
  regions,
  key,
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
  matlist$function.parameters <- params
  
  # assigning a list using SetAssayData will overwrite the whole slot
  # temporary solution
  if (key %in% names(GetAssayData(object, layer = "positionEnrichment"))) {
    warning("Requested name is already present, overwriting existing data")
  }
  object@positionEnrichment[[key]] <- matlist
  return(object)
}

#' @method RegionMatrix default
#' @export
#' @concept heatmap
#' @rdname RegionMatrix
#' @importFrom GenomicRanges resize strand
#' @importFrom fastmatch fmatch
RegionMatrix.default <- function(
  object,
  regions,
  key,
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
  all.valid <- sapply(X = object, FUN = inherits, what = "Fragment")
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
  unique.groups <- unique(x = group.by)
  
  on_plus <- strand(x = regions) == "+" | strand(x = regions) == "*"
  plus.strand <- regions[on_plus, ]
  minus.strand <- regions[!on_plus, ]
  for (i in seq_along(along.with = object)) {
    tmplist <- list()
    
    # open tabix connection
    fragfile <- GetFragmentData(object = object[[i]], slot = "path")
    cellnames <- GetFragmentData(object = object[[i]], slot = "cells")
    tabix.file <- TabixFile(
      file = fragfile,
      index = GetIndexFile(fragment = fragfile, verbose = FALSE)
    )
    open(con = tabix.file)
    
    # initialize empty matrix for each group of cells
    # each row is region
    # each column is a base in region (TODO: implement window sum to reduce size of matrix)
    ncol.mat <- upstream + downstream + 1
    for (j in unique.groups) {
      # TODO make this dgCMatrix instead
      tmplist[[j]] <- matrix(
        data = 0,
        nrow = length(x = regions),
        ncol = ncol.mat
      )
    }
    
    if (length(x = regions) > 0) {
      frags <- scanTabix(file = tabix.file, param = regions)
      res <- TabixOutputToDataFrame(
        reads = frags, record.ident = TRUE
      )
      
      # assign counts to cell groups
      for (j in unique(x = res$ident)) { # for each region
        res_region <- res[res$ident == j, ]
        # subtract start from fragment position
        on_plus <- as.logical(strand(x = regions[j]) == "+" |
                                strand(x = regions[j]) == "*")
        res_region$start <- res_region$start - start(x = regions[j])
        res_region$end <- res_region$end - start(x = regions[j])
        
        # remove out of bounds positions
        res_region <- res_region[
          res_region$start > 0 & res_region$start < ncol.mat, , drop = FALSE
        ]
        res_region <- res_region[
          res_region$end > 0 & res_region$end < ncol.mat, , drop = FALSE
        ]
        for (cell in unique.groups) {
          cells.keep <- names(x = group.by[group.by == cell])
          subfrag <- res_region[
            fmatch(
              x = res_region$cell,
              table = cellnames[cells.keep],
              nomatch = 0L
            ) > 0, ,
            drop = FALSE]
          if (nrow(x = subfrag) > 0) {
            if (on_plus) {
              startpos <- subfrag$start
              endpos <- subfrag$end
            } else {
              # reverse for minus strand
              startpos <- total_bases - subfrag$start
              endpos <- total_bases - subfrag$end
            }
            tmplist[[cell]][j, startpos] <- tmplist[[cell]][j, startpos] + 1
            tmplist[[cell]][j, endpos] <- tmplist[[cell]][j, endpos] + 1
          }
        }
      }
    }
    
    if (i == 1) {
      # one fragment file
      matlist <- tmplist
    } else {
      # sum acros fragment files
      for (cell in unique.groups) {
        matlist[[cell]] <- matlist[[cell]] + tmplist[[cell]]
      }
    }
  }
  return(matlist)
}
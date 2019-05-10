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

# Resize GenomicRanges upstream and or downstream
# from https://support.bioconductor.org/p/78652/
#
Extend <- function(x, upstream = 0, downstream = 0) {
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
  new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
  new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}

#' Add genomic annotation information to a Seurat object
#'
#' @param object A Seurat object
#' @param annotation A genomic feature annotation
#' @param key Key used to store the annotation information, eg "gene", "TE", "promoter", etc.
#'
#' @export
#'
SetAnnotations <- function(
  object,
  annotation,
  key
) {
  # TODO check what the input format should be, eg something from bioconductor. Will depend on plotting requirements
  current.tools <- slot(object = object, name = 'tools')
  if ('annotations' %in% names(x = current.tools)) {
    current.tools$annotations[[key]] <- annotation
  } else {
    current.tools$annotations <- list(key = annotation)
  }
  slot(object = object, name = 'tools') <- current.tools
  return(object)
}

#' Set the fragments file path for creating plots
#'
#' Give path of indexed fragments file that goes with data in the object.
#' Checks for a valid path and an index file with the same name (.tbi) at the same path.
#' Stores the path under the tools slot for access by visualization functions.
#' One fragments file can be stored for each assay.
#'
#' @param object A Seurat object
#' @param file Path to indexed fragment file. See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}
#' @param assay Assay used to generate the fragments. If NULL, use the active assay.
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

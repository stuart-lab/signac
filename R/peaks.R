#' @include generics.R
#'
NULL

#' @param object A Seurat object, ChromatinAssay object, Fragment object, or the
#' path to fragment file/s.
#' @param assay Name of assay to use
#' @param macs2.path Path to MACS program. If NULL, try to find MACS
#' automatically.
#' @param outdir Path for output files
#' @param effective.genome.size Effective genome size parameter for MACS
#' (\code{-g}). Default is the human effective genome size (2.7e9).
#' @param extsize \code{extsize} parameter for MACS.
#' @param additional.args Additional arguments passed to MACS. This should be a
#' single character string
#' @param name Name for output MACS files. This will also be placed in the
#' \code{name} field in the GRanges output.
#' @param cleanup Remove MACS output files
#' @param verbose Display messages
#' @param ... Arguments passed to other methods
#'
#' @method CallPeaks Seurat
#' @rdname CallPeaks
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Project
#'
#' @export
CallPeaks.Seurat <- function(
  object,
  assay = NULL,
  macs2.path = NULL,
  outdir = tempdir(),
  effective.genome.size = 2.7e9,
  extsize = 200,
  additional.args = NULL,
  name = Project(object),
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  gr <- CallPeaks(
    object = object[[assay]],
    macs2.path = macs2.path,
    outdir = outdir,
    effective.genome.size = effective.genome.size,
    extsize = extsize,
    additional.args = additional.args,
    name = name,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  return(gr)
}

#' @method CallPeaks ChromatinAssay
#' @rdname CallPeaks
#' @export
CallPeaks.ChromatinAssay <- function(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  effective.genome.size = 2.7e9,
  extsize = 200,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  # get fragment files
  frags <- Fragments(object = object)
  # get all fragment file paths
  allfragpaths <- sapply(X = frags, FUN = GetFragmentData, slot = "path")
  allfragpaths <- Reduce(f = paste, x = allfragpaths)
  gr <- CallPeaks(
    object = allfragpaths,
    macs2.path = macs2.path,
    outdir = outdir,
    effective.genome.size = effective.genome.size,
    extsize = extsize,
    additional.args = additional.args,
    name = name,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  return(gr)
}

#' @method CallPeaks Fragment
#' @rdname CallPeaks
#' @export
CallPeaks.Fragment <- function(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  effective.genome.size = 2.7e9,
  extsize = 200,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  fragpath <- GetFragmentData(object = object, slot = "path")
  gr <- CallPeaks(
    object = fragpath,
    macs2.path = macs2.path,
    outdir = outdir,
    effective.genome.size = effective.genome.size,
    extsize = extsize,
    additional.args = additional.args,
    name = name,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  return(gr)
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @method CallPeaks default
#' @rdname CallPeaks
#' @export
CallPeaks.default <- function(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  effective.genome.size = 2.7e9,
  extsize = 200,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  # find macs2
  macs2.path <- SetIfNull(
    x = macs2.path,
    y = unname(obj = Sys.which(names = "macs2"))
  )
  if (nchar(x = macs2.path) == 0) {
    stop("MACS2 not found. Please install MACS:",
         "https://macs3-project.github.io/MACS/")
  }

  # if list of paths given, collapse to a single space-separated string
  if (length(x = object) > 1) {
    object <- Reduce(f = paste, x = object)
  }

  cmd <- paste0(
    macs2.path,
    " callpeak -t ",
    object,
    " -g ",
    as.character(x = effective.genome.size),
    " -f BED --nomodel --extsize ",
    as.character(x = extsize),
    " -n ",
    as.character(x = name),
    " --outdir ",
    outdir,
    " ",
    additional.args
  )

  # call macs2
  system(
    command = cmd,
    wait = TRUE,
    ignore.stderr = !verbose,
    ignore.stdout = !verbose
  )

  # read in narrowpeak file and create granges
  df <- read.table(
    file = paste0(outdir, "/", name, "_peaks.narrowPeak"),
    col.names = c("chr", "start", "end", "name",
                  "score", "strand", "fold_change",
                  "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                  "relative_summit_position")
  )
  gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  if (cleanup) {
    files.to.remove <- paste0(
      name,
      c("_peaks.narrowPeak", "_peaks.xls", "_summits.bed")
    )
    files.to.remove <- paste0(outdir, "/", files.to.remove)
    for (i in files.to.remove) {
      if (file.exists(i)) {
        file.remove(i)
      }
    }
  }
  return(gr)
}

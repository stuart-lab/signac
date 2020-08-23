
#' Call peaks
#'
#' Call peaks using MACS. Fragment files linked to the specified assay will be
#' used to call peaks. If multiple fragment files are present, all will be used
#' in a single MACS invocation. Returns the \code{.narrowPeak} MACS output as a
#' \code{GRanges} object.
#'
#' See \url{https://macs3-project.github.io/MACS/} for MACS documentation.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param macs2.path Path to MACS program. If NULL, try to find MACS
#' automatically.
#' @param outpath Path for output files. If NULL, use temp directory and remove
#' files after peak calling.
#' @param effective.genome.size Effective genome size parameter for MACS
#' (\code{-g}). Default is the human effective genome size (2.7e9).
#' @param extsize \code{extsize} parameter for MACS.
#' @param additional.args Additional arguments passed to MACS. This should be a
#' single character string
#' @param name Name for output MACS files
#' @param cleanup Remove MACS output files
#' @param verbose Display messages
#'
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom Seurat Project
#'
CallPeaks <- function(
  object,
  assay = NULL,
  macs2.path = NULL,
  outpath = tempdir(),
  effective.genome.size = 2.7e9,
  extsize = 200,
  additional.args = NULL,
  name = Project(object),
  cleanup = TRUE,
  verbose = TRUE
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
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  # get fragment files
  frags <- Fragments(object = object[[assay]])
  # get all fragment file paths
  allfragpaths <- sapply(X = frags, FUN = GetFragmentData, slot = "path")
  allfragpaths <- Reduce(f = paste, x = allfragpaths)

  cmd <- paste0(
    macs2.path,
    " callpeak -t ",
    allfragpaths,
    " -g ",
    as.character(x = effective.genome.size),
    " -f BED --nomodel --extsize ",
    as.character(x = extsize),
    " -n ",
    as.character(x = name),
    " --outdir ",
    outpath,
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
    file = paste0(outpath, "/", name, "_peaks.narrowPeak"),
    col.names = c("chr", "start", "end", "name",
                  "score", "strand", "fold_change",
                  "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                  "relative_summit_position")
  )
  gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  if (cleanup) {
    files.to.remove <- paste0(
      name,
      c("_peaks.narrowPeak", "_peaks.xls", "_summits.bed"
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

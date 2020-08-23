
#' Call peaks
#'
#' Call peaks using MACS2
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param macs2.path Path to MACS2 program. If NULL, try to find MACS2
#' automatically.
#' @param outpath Path for output files. If NULL, use temp directory and remove
#' files after peak calling.
#' @param effective.genome.size Effective genome size parameter for MACS2
#' (\code{-g}). Default is the human effective genome size (2.7e9).
#' @param extsize \code{extsize} parameter for MACS2.
#'
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#'
#' @importFrom Seurat DefaultAssay
#'
CallPeaks <- function(
  object,
  assay = NULL,
  macs2.path = NULL,
  outpath = tempdir(),
  effective.genome.size = 2.7e9,
  extsize = 200
) {
  # find macs2
  macs2.path <- SetIfNull(
    x = macs2.path,
    y = unname(obj = Sys.which(names = "macs2"))
  )
  if (nchar(x = macs2.path) == 0) {
    stop("MACS2 not found. Please install MACS2:",
         "https://macs3-project.github.io/MACS/")
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  # get fragment files
  # TODO just running on first frag file to test. need to iterate over all and merge peaks
  frags <- Fragments(object = object[[assay]])
  fragpath <- GetFragmentData(object = frags[[1]], slot = "path") # TODO update to iterate over all files

  cmd <- paste0(
    macs2.path,
    " callpeak -t ",
    fragpath,
    " -g ",
    as.character(x = effective.genome.size),
    " -f BED --nomodel --extsize ",
    as.character(x = extsize),
    " -n test "
  )

  # call macs2
  system(
    command = cmd,
    wait = TRUE
  )
  # return()
}

#' @include generics.R
#'
NULL

#' @param object A Seurat object, ChromatinAssay object, Fragment object, or the
#' path to fragment file/s.
#' @param assay Name of assay to use
#' @param group.by Grouping variable to use. If set, peaks will be called
#' independently on each group of cells and then combined. Note that to call
#' peaks using subsets of cells we first split the fragment file/s used, so
#' using a grouping variable will require extra time to split the files and
#' perform multiple MACS peak calls, and will store additional files on-disk
#' that may be large. Note that we store split fragment files in the temp
#' directory (\code{\link[base]{tempdir}}) by default, and if the program is
#' interrupted before completing these temporary files will not be removed. If
#' NULL, peaks are called using all cells together (pseudobulk).
#' @param idents List of identities to include if grouping cells (only valid if
#' also setting the \code{group.by} parameter). If NULL, peaks will be called
#' for all cell identities.
#' @param macs2.path Path to MACS program. If NULL, try to find MACS
#' automatically.
#' @param combine.peaks Controls whether peak calls from different groups of
#' cells are combined using \code{GenomicRanges::reduce} when calling peaks for
#' different groups of cells (\code{group.by} parameter). If FALSE, a list of
#' \code{GRanges} object will be returned. Note that metadata fields such as the
#' p-value, q-value, and fold-change information for each peak will be lost if
#' combining peaks.
#' @param broad Call broad peaks (\code{--broad} parameter for MACS)
#' @param format File format to use. Should be either "BED" or "BEDPE" (see 
#' MACS documentation).
#' @param outdir Path for output files
#' @param fragment.tempdir Path to write temporary fragment files. Only used if
#' \code{group.by} is not NULL.
#' @param effective.genome.size Effective genome size parameter for MACS
#' (\code{-g}). Default is the human effective genome size (2.7e9).
#' @param extsize \code{extsize} parameter for MACS. Only relevant if 
#' format="BED"
#' @param shift \code{shift} parameter for MACS. Only relevant if format="BED"
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
#' @concept quantification
#'
#' @importFrom Seurat DefaultAssay
#' @importFrom Seurat Project
#' @importFrom GenomicRanges reduce
#'
#' @export
CallPeaks.Seurat <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  macs2.path = NULL,
  broad = FALSE,
  format = "BED",
  outdir = tempdir(),
  fragment.tempdir = tempdir(),
  combine.peaks = TRUE,
  effective.genome.size = 2.7e9,
  extsize = 200,
  shift = -extsize/2,
  additional.args = NULL,
  name = Project(object),
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  if (!dir.exists(paths = outdir)) {
    stop("Requested output directory does not exist")
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!is.null(x = group.by)) {
    # first check macs2 path before we spend time splitting the files
    macs2.path <- SetIfNull(
      x = macs2.path,
      y = unname(obj = Sys.which(names = "macs2"))
    )
    if (nchar(x = macs2.path) == 0) {
      stop("MACS2 not found. Please install MACS:",
           "https://macs3-project.github.io/MACS/")
    }
    if (fragment.tempdir != tempdir()) {
      if (!dir.exists(paths = fragment.tempdir)) {
        warning("Requested output directory does not exist, creating directory")
        dir.create(path = fragment.tempdir)
      }
    }
    # split fragment files
    SplitFragments(
      object = object,
      assay = assay,
      group.by = group.by,
      idents = idents,
      outdir = fragment.tempdir,
      verbose = verbose
    )

    # work out what all the file paths to use are
    groups <- GetGroups(
      object = object,
      group.by = group.by,
      idents = idents
    )
    groups <- gsub(pattern = " ", replacement = "_", x = groups)
    groups <- gsub(pattern = .Platform$file.sep, replacement = "_", x = groups)
    unique.groups <- unique(x = groups)

    # call peaks on each split fragment file separately
    grlist <- list()
    for (i in seq_along(along.with = unique.groups)) {
      fragpath <- paste0(
        fragment.tempdir,
        .Platform$file.sep,
        unique.groups[[i]],
        ".bed"
      )
      gr <- CallPeaks(
        object = fragpath,
        macs2.path = macs2.path,
        outdir = outdir,
        broad = broad,
        format = format,
        effective.genome.size = effective.genome.size,
        extsize = extsize,
        shift = shift,
        additional.args = additional.args,
        name = unique.groups[[i]],
        cleanup = cleanup,
        verbose = verbose
      )
      # remove split fragment file from temp dir
      file.remove(fragpath)
      # add ident
      if (length(x = gr) > 0) {
        gr$ident <- unique.groups[[i]]
        grlist[[i]] <- gr
      } else {
        message("No peaks found for ", unique.groups[[i]])
      }
    }
    if (combine.peaks) {
      # combine peaks and reduce, maintaining ident information
      gr.combined <- Reduce(f = c, x = grlist)
      gr <- reduce(x = gr.combined, with.revmap = TRUE)
      dset.vec <- vector(mode = "character", length = length(x = gr))
      ident.vec <- gr.combined$ident
      revmap <- gr$revmap
      for (i in seq_len(length.out = length(x = gr))) {
        datasets <- ident.vec[revmap[[i]]]
        dset.vec[[i]] <- paste(unique(x = datasets), collapse = ",")
      }
      gr$peak_called_in <- dset.vec
      gr$revmap <- NULL
    } else {
      gr <- grlist
    }
  } else {
    gr <- CallPeaks(
      object = object[[assay]],
      macs2.path = macs2.path,
      outdir = outdir,
      broad = broad,
      format = format,
      effective.genome.size = effective.genome.size,
      extsize = extsize,
      shift = shift,
      additional.args = additional.args,
      name = name,
      cleanup = cleanup,
      verbose = verbose,
      ...
    )
  }
  return(gr)
}

#' @method CallPeaks ChromatinAssay
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.ChromatinAssay <- function(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  broad = FALSE,
  format = "BED",
  effective.genome.size = 2.7e9,
  extsize = 200,
  shift = -extsize/2,
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
  gr <- CallPeaks(
    object = allfragpaths,
    macs2.path = macs2.path,
    outdir = outdir,
    broad = broad,
    format = format,
    effective.genome.size = effective.genome.size,
    extsize = extsize,
    shift = shift,
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
#' @concept quantification
#' @export
CallPeaks.Fragment <- function(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  broad = FALSE,
  format = "BED",
  effective.genome.size = 2.7e9,
  extsize = 200,
  shift = -extsize/2,
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
    broad = broad,
    format = format,
    effective.genome.size = effective.genome.size,
    extsize = extsize,
    shift = shift,
    additional.args = additional.args,
    name = name,
    cleanup = cleanup,
    verbose = verbose,
    ...
  )
  return(gr)
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils read.table
#' @method CallPeaks default
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.default <- function(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  broad = FALSE,
  format = "BED",
  effective.genome.size = 2.7e9,
  extsize = 200,
  shift = -extsize/2,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
  if (!dir.exists(paths = outdir)) {
    stop("Requested output directory does not exist")
  }
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
    object <- sapply(
      X = object, FUN = function(x) paste0("'", x, "'"), USE.NAMES = FALSE
    )
    object <- Reduce(f = paste, x = object)
  } else {
    object <- paste0("'", object, "'")
  }

  broadstring <- ifelse(test = broad, yes = " --broad ", no = "")
  nomod_str <- ifelse(
    test = format == "BED",
    yes = paste0(" --nomodel --extsize ",
    as.character(x = extsize),
    " --shift ",
    as.character(x = shift)
    ),
    no = ""
  )
  
  cmd <- paste0(
    macs2.path,
    " callpeak -t ",
    object,
    " -g ",
    as.character(x = effective.genome.size),
    broadstring,
    " -f ",
    format,
    nomod_str,
    " -n ",
    "'",
    as.character(x = name),
    "'",
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

  if (broad) {
    # read in broadpeak
    df <- read.table(
      file = paste0(outdir, .Platform$file.sep, name, "_peaks.broadPeak"),
      col.names = c("chr", "start", "end", "name",
                    "score", "strand", "fold_change",
                    "neg_log10pvalue_summit", "neg_log10qvalue_summit")
    )
    files.to.remove <- paste0(
      name,
      c("_peaks.broadPeak", "_peaks.xls", "_peaks.gappedPeak")
    )
  } else {
    # read in narrowpeak file
    df <- read.table(
      file = paste0(outdir, .Platform$file.sep, name, "_peaks.narrowPeak"),
      col.names = c("chr", "start", "end", "name",
                    "score", "strand", "fold_change",
                    "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                    "relative_summit_position")
    )
    files.to.remove <- paste0(
      name,
      c("_peaks.narrowPeak", "_peaks.xls", "_summits.bed")
    )
  }

  gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  if (cleanup) {
    files.to.remove <- paste0(outdir, .Platform$file.sep, files.to.remove)
    for (i in files.to.remove) {
      if (file.exists(i)) {
        file.remove(i)
      }
    }
  }
  return(gr)
}

#' @include generics.R
#'
NULL

# find macs3 & check if macs3 is executable
macs3_pathcheck <- function(macs3.path, mode) {
    macs3.path <- macs3.path %||% unname(obj = Sys.which(names = "macs3"))
    if (nchar(x = macs3.path) == 0) {
        stop(
            "MACS3 not found. Please install MACS3: ",
            "https://macs3-project.github.io/MACS/"
        )
    } else if (file.access(names = macs3.path, mode = 0) == -1) {
        stop(
            "MACS3 does not exist at specified path. Please install MACS3: ",
            "https://macs3-project.github.io/MACS/"
        )
    } else if (file.access(names = macs3.path, mode = 1) == -1) {
        stop("MACS3 exists but is not executable")
    }

    if (mode == 'hmmratac') {
        version <- system2("macs3", args = "--version", stdout = TRUE)
        version_n <- gsub(".*macs3 ", "", version)
        if (package_version(version_n) < "3.0.4") {
            stop("hmmratac mode requires macs3 >= 3.0.4. Please upgrade your macs3 installation or use `callpeak` mode.")
        }
    }
    
    return(macs3.path)
}

#' @importFrom GenomicRanges reduce
#' @importFrom IRanges extractList
#' @importFrom S4Vectors unstrsplit
CombinePeaks <- function(grlist) {
    # combine peaks and reduce, maintaining ident information
    if (is(object = grlist, class2 = "GRangesList")) {
      gr.combined <- unlist(x = grlist)
    } else {
      gr.combined <- do.call(c, grlist)
    }
    gr <- reduce(x = gr.combined, with.revmap = TRUE)
    ids_char <- as.character(x = gr.combined$ident)
    ids_list <- extractList(x = ids_char, i = gr$revmap)
    ids_list <- unique(x = ids_list)
    gr$peak_called_in <- unstrsplit(ids_list, sep = ",")
    gr$revmap <- NULL
    return(gr)
}

#' @param object A [SeuratObject::Seurat] object, [ChromatinAssay5-class]
#' object, [Fragment2-class] object, or the path to fragment file/s.
#' @param assay Name of assay to use.
#' @param group.by Grouping variable to use. If set, peaks will be called
#' independently on each group of cells. If multiple fragments are present in
#' the object,  peaks will be called separately for each fragment file and for
#' each group of cells. Note that to call peaks using a subset of cells, we
#' write the cell barcodes into `.txt` file/s that are stored in the temp
#' directory ([base::tempdir()]) by default and removed upon completion (if
#' `cleanup=TRUE`). If `NULL`, peaks are called using all cells.
#' @param idents List of identities to include if grouping cells (only valid if
#' also setting the `group.by` parameter). If `NULL`, peaks will be called
#' for all cell identities.
#' @param macs3.path Path to MACS3 program. If `NULL`, try to find MACS3
#' automatically.
#' @param mode MACS function to call, choose between `callpeak` or `hmmratac`.
#' Default is `callpeak`. If using `hmmratac` mode, runtime may be longer than
#' `callpeak` and additional arguments may be needed. See MACS HMMRATAC documentation 
#' for more details: https://deepwiki.com/macs3-project/MACS/3.2-atac-seq-analysis-with-hmmratac
#' @param combine.peaks Controls whether peak calls from different groups of
#' cells are combined using [GenomicRanges::reduce()] when calling peaks for
#' different groups of cells (`group.by` parameter). If `FALSE`, a list of
#' [GenomicRanges::GRanges] object will be returned. Note that metadata fields
#' such as the p-value, q-value, and fold-change information for each peak will
#' be lost if combining peaks.
#' @param broad Call broad peaks (`--broad` parameter for MACS).
#' @param outdir Path for output files.
#' @param cells Vector of cell barcodes to call peaks on. If `NULL`, use all
#' cells present in the fragment files, unless `group.by` or `idents` is set.
#' @param genome MACS3 built-in effective genome size. Default is `hs` (Human,
#' GRCh38). Genome sizes for `mm` (Mice, GRCm38), `ce` (C. elegans, WBcel235),
#' and `dm` (Drosophila M., dm6) are also available.
#' @param gsize Manually set effective genome size parameter. If specified,
#' overrides MACS3 built-in genome sizes.
#' @param additional.args Additional arguments passed to MACS. This should be a
#' single character string.
#' @param name Name for output MACS files. This will also be placed in the
#' `name` field in the GRanges output.
#' @param cleanup Remove MACS output files.
#' @param verbose Display messages.
#' @param ... Arguments passed to other methods.
#'
#' @method CallPeaks Seurat
#' @rdname CallPeaks
#'
#' @concept quantification
#'
#' @importFrom SeuratObject DefaultAssay Project
#' @importFrom GenomicRanges reduce
#'
#' @export
CallPeaks.Seurat <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  cells = NULL,
  macs3.path = NULL,
  mode = "callpeak",
  broad = FALSE,
  genome = "hs",
  gsize = NULL,
  outdir = tempdir(),
  combine.peaks = TRUE,
  additional.args = NULL,
  name = Project(object),
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
    if (!dir.exists(paths = outdir)) {
        stop("Requested output directory does not exist")
    }
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path, mode = mode)

    # check macs3 mode
    if (!mode %in% c("callpeak", "hmmratac")) {
        stop("Invalid macs3 command, choose between `callpeak` or `hmmratac`")
    }

    # check object assay
    assay <- assay %||% DefaultAssay(object = object)

    # check group.by & idents
    if (is.null(x = group.by) && !is.null(x = idents)) {
        stop("Set group.by parameter if calling peaks per ident")
    }

    if (!is.null(x = group.by)) {
        # get cell barcodes and the group they belong to
        groups <- GetGroups(object = object, group.by = group.by, idents = idents)

        # filter to list of cells
        if (!is.null(x = cells)) {
            groups <- groups[names(x = groups) %in% cells]
        }

        # iterate over cell groups and combine outputs
        unique_groups <- unique(x = groups)

        # check parallelization
        if (length(x = unique_groups) > 1) {
            # check n workers for parallelization
            if (future::nbrOfWorkers() > 1) {
                mylapply <- future_lapply
            } else {
                mylapply <- ifelse(
                    test = verbose,
                    yes = pbapply::pblapply,
                    no = lapply
                )
            }
        } else {
            mylapply <- lapply
        }

        pk.all <- mylapply(
            X = unique_groups,
            FUN = function(x) {
                gr <- CallPeaks(
                    object = object[[assay]],
                    cells = names(x = groups[groups == x]),
                    macs3.path = macs3.path,
                    combine.peaks = combine.peaks,
                    mode = mode,
                    outdir = outdir,
                    broad = broad,
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = x,
                    cleanup = cleanup,
                    verbose = FALSE
                )
                gr$ident <- x
                gr
            }
        )
        if (combine.peaks == TRUE) {
            peakcalls <- CombinePeaks(grlist = pk.all)
        } else {
            peakcalls <- pk.all
        }
    } else {
        peakcalls <- CallPeaks(
            object = object[[assay]],
            cells = cells,
            macs3.path = macs3.path,
            combine.peaks = combine.peaks,
            mode = mode,
            outdir = outdir,
            broad = broad,
            genome = genome,
            gsize = gsize,
            additional.args = additional.args,
            name = name,
            cleanup = cleanup,
            verbose = verbose
        )
    }
    return(peakcalls)
}

#' @method CallPeaks ChromatinAssay5
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.ChromatinAssay5 <- function(
  object,
  macs3.path = NULL,
  combine.peaks = TRUE,
  mode = "callpeak",
  outdir = tempdir(),
  broad = FALSE,
  cells = NULL,
  genome = "hs",
  gsize = NULL,
  additional.args = NULL,
  name = "macs3",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path, mode = mode)
    frags <- Fragments(object = object)

    if (is.null(x = cells)) {
        # cells not set, one combined MACS3 call with all fragment files
        allfragpaths <- lapply(X = frags, FUN = GetFragmentData, slot = "file.path")
        allfragpaths <- paste(allfragpaths, collapse = " ")
        peakcalls <- CallPeaks(
            object = allfragpaths,
            macs3.path = macs3.path,
            mode = mode,
            outdir = outdir,
            broad = broad,
            genome = genome,
            gsize = gsize,
            additional.args = additional.args,
            name = name,
            cleanup = cleanup,
            verbose = verbose
        )
    } else {
        # call out to the Fragment2 method on each fragment file
        pk.all <- list()
        for (i in seq_along(along.with = frags)) {
            gr <- CallPeaks(
                object = frags[[i]],
                cells = cells,
                macs3.path = macs3.path,
                mode = mode,
                outdir = outdir,
                broad = broad,
                genome = genome,
                gsize = gsize,
                additional.args = additional.args,
                name = paste0(name, as.character(x = i)),
                cleanup = cleanup,
                verbose = verbose
            )
            gr$ident <- i
            pk.all[[i]] <- gr
        }

        # combine output
        if (length(x = pk.all) > 1) {
            if (combine.peaks == TRUE) {
                peakcalls <- CombinePeaks(grlist = pk.all)
            } else {
                peakcalls <- pk.all
            }
        } else {
          peakcalls <- pk.all[[1]]
          peakcalls$ident <- NULL
        }
    }
    return(peakcalls)
}

#' @method CallPeaks Fragment2
#' @rdname CallPeaks
#' @concept quantification
#' @importFrom GenomicRanges GRanges
#' @export
CallPeaks.Fragment2 <- function(
  object,
  macs3.path = NULL,
  mode = "callpeak",
  outdir = tempdir(),
  broad = FALSE,
  cells = NULL,
  genome = "hs",
  gsize = NULL,
  additional.args = NULL,
  name = "macs3",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path, mode = mode)
    
    fpath <- GetFragmentData(object = object, slot = "file.path")

    # write cell barcodes file
    if (!is.null(x = cells)) {
        cell_barcodes <- GetFragmentData(object = object, slot = "cells")
        cell_barcodes <- unname(
            obj = cell_barcodes[names(x = cell_barcodes) %in% cells]
        )
        if (length(x = cell_barcodes) < 5) {
          warning("Insufficient cells present in fragment file: ", fpath)
          gr <- GRanges()
          return(gr)
        }
        barcodes <- paste0(
            outdir, .Platform$file.sep, paste0(name, "_barcodes.txt")
        )
        barcodes <- gsub(pattern = " ", replacement = "_", x = barcodes)
        writeLines(text = cell_barcodes, con = barcodes)
    } else {
      barcodes <- NULL
    }

    gr <- CallPeaks(
        object = fpath,
        macs3.path = macs3.path,
        mode = mode,
        outdir = outdir,
        broad = broad,
        barcodes = barcodes,
        genome = genome,
        gsize = gsize,
        additional.args = additional.args,
        name = name,
        cleanup = cleanup,
        verbose = verbose
    )

    return(gr)
}

#' @param barcodes Path to cell barcodes (`--barcodes` parameter for MACS).
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils read.table
#' @method CallPeaks default
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.default <- function(
  object,
  macs3.path = NULL,
  mode = "callpeak",
  outdir = tempdir(),
  broad = FALSE,
  barcodes = NULL,
  genome = "hs",
  gsize = NULL,
  additional.args = NULL,
  name = "macs3",
  cleanup = TRUE,
  verbose = TRUE,
  ...
) {
    if (!dir.exists(paths = outdir)) {
        stop("Requested output directory does not exist")
    }
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path, mode = mode)

    if (nchar(x = object) == 0 || object == " ") {
        stop("Empty path given for fragment file")
    }

    name <- gsub(pattern = " ", replacement = "_", x = name)
    name <- gsub(pattern = .Platform$file.sep, replacement = "_", x = name)

    # check macs3 mode
    if (!mode %in% c("callpeak", "hmmratac")) {
        stop("Invalid macs3 command, choose between `callpeak` or `hmmratac`")
    }

    # if list of paths given, collapse to a single space-separated string
    if (length(x = object) > 1) {
        object <- paste(" ", object, collapse = " ")
    }

    # buffer with whitespace
    object <- paste0(" ", object, " ")

    # check genome format
    if (!is.null(x = genome) && genome %in% c("hs", "mm", "ce", "dm")) {
        genome_string <- paste0(" -g ", genome, " ")
    } else if (is.null(x = genome) && !is.null(x = gsize) && !is.na(x = gsize)) {
        if (is.na(x = as.numeric(x = gsize))) {
            stop("Requested non-numeric gsize value")
        }
        genome_string <- paste0(" --gsize ", gsize, " ")
    } else {
        stop(
            "Invalid genome size, choose between `hs` (human, GRCh38)",
            "`mm` (mice, GRCm38), `ce` (c elegans, WBcel235), `dm` (drosophila m, dm6) for MACS3 built-in genome size,",
            " or manually set genome size with the gzise parameter"
        )
    }

    # broad setting
    broad_string <- ifelse(test = broad, yes = " --broad ", no = "")

    # add cell barcode
    barcode_string <- ifelse(test = !is.null(x = barcodes),
        yes = paste0(" --barcodes '", barcodes, "'"),
        no = ""
    )

    # object path
    if (mode == "callpeak") {
        object_string <- paste0(" -t ", object)
    } else if (mode == "hmmratac") {
        object_string <- paste0(" -i ", object)
        genome_string <- " "
        broad_string <- " "
        if (verbose) {
            message(paste0(
            "`hmmratac` mode selected. This will run slower than `callpeak` mode."
            )) 
        }
    } else {
        stop("Invalid macs3 mode")
    }

    # macs3 command
    cmd <- paste0(
        macs3.path, " ",
        mode,
        object_string,
        genome_string,
        broad_string,
        " -f FRAG ",
        barcode_string,
        " -n ", "'", as.character(x = name), "'",
        " --outdir ", outdir,
        " ",
        additional.args
    )
    
    # call macs3
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
            col.names = c(
                "chr", "start", "end", "name",
                "score", "strand", "fold_change",
                "neg_log10pvalue_summit", "neg_log10qvalue_summit"
            )
        )
        files.to.remove <- paste0(
            name,
            c("_peaks.broadPeak", "_peaks.xls", "_peaks.gappedPeak")
        )
    } else if (mode == "hmmratac") {
        # read in narrowpeak file
        df <- read.table(
            file = paste0(outdir, .Platform$file.sep, name, "_accessible_regions.narrowPeak"),
            col.names = c(
                "chr", "start", "end", "name",
                "score", "strand", "fold_change",
                "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                "relative_summit_position"
            )
        )
        files.to.remove <- paste0(
            name,
            c("_accessible_regions.narrowPeak", "_model.json", "_cutoff_analysis.tsv")
        )
    } else {
        # read in narrowpeak file
        df <- read.table(
            file = paste0(outdir, .Platform$file.sep, name, "_peaks.narrowPeak"),
            col.names = c(
                "chr", "start", "end", "name",
                "score", "strand", "fold_change",
                "neg_log10pvalue_summit", "neg_log10qvalue_summit",
                "relative_summit_position"
            )
        )
        files.to.remove <- paste0(
            name,
            c("_peaks.narrowPeak", "_peaks.xls", "_summits.bed")
        )
    }

    gr <- makeGRangesFromDataFrame(
        df = df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
    )

    if (cleanup) {
        # remove macs3 files
        files.to.remove <- paste0(outdir, .Platform$file.sep, files.to.remove)
        for (i in files.to.remove) {
            if (file.exists(i)) {
                file.remove(i)
            }
        }
        # remove barcode file
        if (!is.null(barcodes)) {
            if (file.exists(barcodes)) {
                file.remove(barcodes)
            }
        }
    }

    return(gr)
}

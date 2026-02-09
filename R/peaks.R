#' @include generics.R
#'
NULL

# find macs3 & check if macs3 is executable
macs3_pathcheck <- function(macs3.path) {
  macs3.path <- macs3.path %||% unname(obj = Sys.which(names = "macs3"))
  if (nchar(x = macs3.path) == 0) {
    stop(
      "MACS3 not found. Please install MACS3: ",
      "https://macs3-project.github.io/MACS/"
    )
  } else if(file.access(names = macs3.path, mode = 0) == -1) {
    stop(
      "MACS3 does not exist at specified path. Please install MACS3: ",
      "https://macs3-project.github.io/MACS/"
    )
  } else if(file.access(names = macs3.path, mode = 1) == -1) {
    stop("MACS3 exists but is not executable")
  }
  return(macs3.path)
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
#' Default is `callpeak`.
#' @param combine.peaks Controls whether peak calls from different groups of
#' cells are combined using [GenomicRanges::reduce()] when calling peaks for
#' different groups of cells (`group.by` parameter). If `FALSE`, a list of
#' [GenomicRanges::GRanges] object will be returned. Note that metadata fields
#' such as the p-value, q-value, and fold-change information for each peak will
#' be lost if combining peaks.
#' @param broad Call broad peaks (`--broad` parameter for MACS).
#' @param outdir Path for output files.
#' @param barcodes Path to cell barcodes (`--barcodes` parameter for MACS).
#' @param cells Vector of cell barcodes to call peaks on.
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
    barcodes = NULL,   
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
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path)
    
    # check macs3 mode
    if (!mode %in% c("callpeak", "hmmratac")) {
      stop("Invalid macs3 command, choose between `callpeak` or `hmmratac`")
    }
    
    # check object assay
    assay <- assay %||% DefaultAssay(object = object)
    
    # check group.by & idents
    if (is.null(group.by) && !is.null(idents)) {
        stop("Set `group.by` parameter if calling peaks per ident")
        }

    if (!is.null(group.by) && is.null(idents)) {
        idents <- unique(object[[group.by]])[, group.by]
    }

    # get number of fragments
    frags <- Fragments(object = object)
    allfragpaths <- as.list(sapply(X = frags, FUN = GetFragmentData, slot = "file.path"))

    # check parallelization
    if (length(allfragpaths) > 1 || !is.null(group.by)) {
        # check n workers for parallelization
        if (future::nbrOfWorkers() > 1) {
            mylapply <- future_lapply
        } else {
            mylapply <- ifelse(test = verbose, yes = pbapply::pblapply, no = lapply)
        }
    } else {
        mylapply <- lapply
    }

    # check callpeaks barcode input
    if (!is.null(group.by)) {
        # get cell barcodes per ident
        group_barcodes <- list()
        for (i in seq_along(idents)) {
            ident_barcodes <- names(GetGroups(object, group.by = group.by, idents = idents[[i]]))

            # if set cells, subset
            if (!is.null(cells)) {
                ident_barcodes <- ident_barcodes[ident_barcodes %in% cells]
            }

            group_barcodes[[i]] <- ident_barcodes
        }
    }

    # call peaks per fragment
    if (!is.null(group.by)) {
        idx <- expand.grid(
            frag = seq_along(frags),
            ident = seq_along(idents)
        )
        peakcalls <- mylapply(
            X = seq_len(nrow(idx)),
            FUN = function(k) {
                i <- idx$frag[k]
                j <- idx$ident[k]

                CallPeaks(
                    object = frags[[i]],
                    macs3.path = macs3.path,
                    mode = mode,
                    outdir = outdir,
                    broad = broad,
                    barcodes = barcodes,
                    cells = group_barcodes[[j]],
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = paste0(name, "_frag", i, "_ident", j),
                    cleanup = TRUE,
                    verbose = TRUE,
                    ...
                )
            }
        )
    } else if (is.null(group.by) && !is.null(cells)) {
        # separate cells per fragment
        frag_cells <- list()
        for (i in seq_along(frags)) {
            frag_cells[[i]] <- names(frags[[i]]@cells[cells])
        }

        # call peaks
        peakcalls <- mylapply(
            X = seq_along(frags),
            FUN = function(i) {
                CallPeaks(
                    object = frags[[i]],
                    macs3.path = macs3.path,
                    mode = mode,
                    outdir = outdir,
                    broad = broad,
                    barcodes = barcodes,
                    cells = frag_cells[[i]],
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = paste0(name,"_frag",i),
                    cleanup = TRUE,
                    verbose = TRUE,
                    ...)
            }
        )
    } else { 
        peakcalls <- mylapply(
            X = 1L,
            FUN = function(i) {
                CallPeaks(
                    object = object[[assay]],
                    macs3.path = macs3.path,
                    combine.peaks = combine.peaks,
                    mode = mode,
                    outdir = outdir,
                    broad = broad,
                    barcodes = barcodes,
                    cells = cells,
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = name,
                    cleanup = cleanup,
                    verbose = verbose,
                    ...
                )
            }
        )
    }

    # combine into 1 granges
    if (combine.peaks == TRUE && length(peakcalls) > 1) {
        peakcalls <- reduce(do.call(c, peakcalls))
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
    barcodes = NULL,
    cells = NULL,
    genome = "hs",
    gsize = NULL,
    additional.args = NULL,
    name = "macs3",
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path)
    # get fragment path(s)
    frags <- Fragments(object = object)
    allfragpaths <- as.list(sapply(X = frags, FUN = GetFragmentData, slot = "file.path"))

    # check number of fragments
    if (length(allfragpaths) > 1) {
        # check number of workers for parallelization
        if (future::nbrOfWorkers() > 1) {
            mylapply <- future_lapply
        } else {
            mylapply <- ifelse(test = verbose, yes = pbapply::pblapply, no = lapply)
        }
    } else {
        mylapply <- lapply
    }

    # write cell barcodes per fragment
    barcode_paths <- list()
    for (i in seq_along(frags)) {
        # get all cell barcodes in fragment
        bc <- frags[[i]]@cells

        if (!is.null(cells)) {
            bc <- frags[[i]]@cells[cells]
        }

        # write barcodes to file
        barcode_path <- paste0(outdir, .Platform$file.sep, paste0(name,"_frag",i,"_barcodes.txt"))
        writeLines(bc, con = barcode_path) 
        barcode_paths[[i]] <- barcode_path
    }

    # clean objects
    rm(object)
    rm(frags)
    rm(bc)
    gc()

    # call peaks
    peakcalls <- mylapply(
        X = seq_along(along.with = allfragpaths),
        FUN = function(i) {
            CallPeaks(
                object = allfragpaths[[i]],
                macs3.path = macs3.path,
                mode = mode,
                outdir = outdir,
                broad = broad,
                barcodes = barcode_paths[[i]],
                genome = genome,
                gsize = gsize,
                additional.args = additional.args,
                name = paste0(name,"_frag",i),
                cleanup = TRUE,
                verbose = TRUE)
        }
    )

    # combine into 1 granges
    if (combine.peaks == TRUE && length(allfragpaths) > 1) {
        peakcalls <- reduce(do.call(c, peakcalls))
    }

    return(peakcalls)
}

#' @method CallPeaks Fragment2
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.Fragment2 <- function( 
    object,
    macs3.path = NULL,
    mode = "callpeak",
    outdir = tempdir(),
    broad = FALSE,
    barcodes = NULL,
    cells = NULL,     
    genome = "hs",
    gsize = NULL,
    additional.args = NULL,
    name = "macs3",
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path)

    # get fragment file paths
    fragpath <- GetFragmentData(object)

    # write cell barcodes file
    if (!is.null(x = cells)) {
        cell_barcodes <- object@cells[cells]
    }
    barcodes <- paste0(outdir, .Platform$file.sep, paste0(name,"_barcodes.txt"))
    writeLines(cell_barcodes, con = barcodes)

    gr <- CallPeaks(
        object = fragpath,
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
        verbose = verbose,
        ...)

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
    macs3.path <- macs3_pathcheck(macs3.path = macs3.path)

    name <- gsub(pattern = " ", replacement = "_", x = name)
    name <- gsub(pattern = .Platform$file.sep, replacement = "_", x = name)

    # check macs3 mode
    if (!mode %in% c("callpeak", "hmmratac")) {
        stop("Invalid macs3 command, choose between `callpeak` or `hmmratac`")
    }

    # if list of paths given, collapse to a single space-separated string
    if (length(x = object) > 1) {
        object <- sapply(
        X = object, FUN = function(x) paste0(" ", x, " "), USE.NAMES = FALSE
        )
        object <- Reduce(f = paste, x = object)
    } else {
        object <- paste0(" ", object, " ")
        if (object == "  ") {
            stop("No fragment files supplied")
        }
    }

    # check genome format    
    if (!is.null(genome) && genome %in% c("hs", "mm", "ce", "dm")) {
        genome_string <- paste0(" -g ", genome, " ")
    } else if (is.null(genome) && !is.null(gsize) && is.numeric(gsize) && length(gsize) == 1 && !is.na(gsize)) {
        genome_string <- paste0(" --gsize ", gsize, " ")
    } else {
        stop("Invalid genome size, choose between `hs` (human, GRCh38)",
            "`mm` (mice, GRCm38), `ce` (c elegans, WBcel235), `dm` (drosophila m, dm6) for MACS3 built in genome size,",
            " or manually set genome size with the gzise parameter")
    }

    # broad setting
    broad_string <- ifelse(test = broad, yes = " --broad ", no = "")
    
    # add cell barcode
    barcode_string <- ifelse(test = !is.null(barcodes), 
                             yes = paste0(" --barcodes ", barcodes), 
                             no = "")
    
    # object path
    if (mode == "callpeak") {
        object_string <- paste0(" -t ",  object)
    } else if (mode == "hmmratac") {
        object_string <- paste0(" -i ", object)
        genome_string <- " "
    } else {
        stop("invalid macs3 mode")
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
        additional.args)

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

    gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
    
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

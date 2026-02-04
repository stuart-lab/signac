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
#' directory ([base::tempdir()]) by default, and if the program is
#' interrupted before completing these temporary files will not be removed. If
#' NULL, peaks are called using all cells together (pseudobulk).
#' @param idents List of identities to include if grouping cells (only valid if
#' also setting the `group.by` parameter). If NULL, peaks will be called
#' for all cell identities.
#' @param macs2.path Path to MACS program. If NULL, try to find MACS
#' automatically.
#' @param combine.peaks Controls whether peak calls from different groups of
#' cells are combined using `GenomicRanges::reduce` when calling peaks for
#' different groups of cells (`group.by` parameter). If FALSE, a list of
#' `GRanges` object will be returned. Note that metadata fields such as the
#' p-value, q-value, and fold-change information for each peak will be lost if
#' combining peaks.
#' @param broad Call broad peaks (`--broad` parameter for MACS)
#' @param format File format to use. Should be either "BED" or "BEDPE" (see
#' MACS documentation).
#' @param outdir Path for output files
#' @param fragment.tempdir Path to write temporary fragment files. Only used if
#' `group.by` is not NULL.
#' @param effective.genome.size Effective genome size parameter for MACS
#' (`-g`). Default is the human effective genome size (2.7e9).
#' @param extsize `extsize` parameter for MACS. Only relevant if
#' format="BED"
#' @param shift `shift` parameter for MACS. Only relevant if format="BED"
#' @param additional.args Additional arguments passed to MACS. This should be a
#' single character string
#' @param name Name for output MACS files. This will also be placed in the
#' `name` field in the GRanges output.
#' @param cleanup Remove MACS output files
#' @param verbose Display messages
#' @param ... Arguments passed to other methods
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
    group.by = NULL,   # call peaks separately per group.by variable (call multiple peaks)
    idents = NULL,     # call peaks by barcodes using idents (one peak call)
    macs3.path = NULL,
    mode = "callpeak", # add hmmratac mode once available
    broad = FALSE,
    genome = "hs",
    gsize = NULL,
    outdir = tempdir(),
    combine.peaks = TRUE,
    additional.args = NULL,
    name = Project(object),
    parallel = TRUE, 
    cleanup = TRUE,
    verbose = TRUE,
    ...
) {
    if (!dir.exists(paths = outdir)) {
        stop("Requested output directory does not exist")
    }
    master_outdir <- outdir
    
    # find macs3
    macs3.path <- macs3.path %||% unname(obj = Sys.which(names = "macs3"))
    if (nchar(x = macs3.path) == 0) {
        stop(
            "MACS3 not found. Please install MACS:",
            "https://macs3-project.github.io/MACS/"
        )
      }

    # get object fragment(s) path
    frags <- Fragments(object = object)
    allfragpaths <- as.list(sapply(X = frags, FUN = GetFragmentData, slot = "file.path"))
    
    assay <- assay %||% DefaultAssay(object = object)

    # call peaks in parallel if using group.by
    if (!is.null(x = group.by)) {
        unique_groups <- unique(object[[group.by]])[,1]

        # split barcodes per group
        barcode_path_list <- c()
        for (i in seq_along(unique_groups)) {
            group <- as.character(unique_groups[i])
            barcodes_vec <- WhichCells(pbmc, idents = group)
            barcodes_vec <- sub("^[^_]*_", "", barcodes_vec) # remove any cell id
            barcode_path <- paste0(master_outdir, .Platform$file.sep, paste0(group, "_barcodes.txt"))
            writeLines(barcodes_vec, con = barcode_path)
            barcode_path_list[[group]] <- barcode_path
        }

        # MACS3 parallelization
        if (parallel) {
            # detect maximum number of cores
            num_groups <- length(unique_groups)
            actual_cpus <- parallelly::availableCores()
            plan(future::multisession, workers = min(num_groups, actual_cpus))

            message(actual_cpus) ### DEBUG
            message(num_groups)  ### DEBUG
        
            # run CallPeaks in parallel
            gr_list <- future_lapply(unique_groups, function(i) {
                CallPeaks(
                    object = unlist(allfragpaths),  ### object fragment(s)
                    macs3.path = macs3.path,
                    mode = mode,
                    outdir = master_outdir,
                    broad = broad,
                    barcodes = barcode_path_list[[unique_groups[i]]],
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = paste0(unique_groups[i], "_", name),
                    cleanup = cleanup,
                    verbose = verbose
                )
            }, future.seed = TRUE,
                future.globals = c(
                    "CallPeaks", # for testing function only
                    "allfragpaths", "macs3.path", "mode", "master_outdir", "broad", 
                    "barcode_path_list", "unique_groups", "genome", "gsize", 
                    "additional.args", "name", "cleanup", "verbose"
            ))
            
            return(gr_list)
            
        } else {
            # If parallel = FALSE, run CallPeaks sequentially
            gr_list <- lapply(n_frags, function(i) {
                CallPeaks(
                    object = unlist(allfragpaths),
                    macs3.path = macs3.path,
                    mode = mode,
                    outdir = master_outdir,
                    broad = broad,
                    barcodes = barcode_path_list[[unique_groups[i]]],
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = paste0(unique_groups[i], "_", name),
                    cleanup = cleanup,
                    verbose = verbose
                )
            })
            return(gr_list)
        }
    }

    # call peaks by ident
    if (!is.null(x = idents)) {
        # get ident barcodes
        barcodes_vec <- WhichCells(pbmc, idents = idents)
        barcodes_vec <- sub("^[^_]*_", "", barcodes_vec) # remove any cell id
        barcode_path <- paste0(master_outdir, .Platform$file.sep, paste0(idents, "_barcodes.txt"))
        writeLines(barcodes_vec, con = barcode_path)
        
        # call peaks
        gr_list <- CallPeaks(
                        object = unlist(allfragpaths),  ### object fragment(s)
                        macs3.path = macs3.path,
                        mode = mode,
                        outdir = master_outdir,
                        broad = broad,
                        barcodes = barcode_path,
                        genome = genome,
                        gsize = gsize,
                        additional.args = additional.args,
                        name = paste0(idents, "_", name),
                        cleanup = cleanup,
                        verbose = verbose
        )

        return(gr_list)
    }

    # write obj cell barcodes
    barcodes_vec <- rownames(object[[assay]]@`cells`)
    barcodes_vec <- sub("^[^_]*_", "", barcodes_vec) # remove any cell id
    barcode_path <- paste0(master_outdir, .Platform$file.sep, paste0("barcodes.txt"))
    writeLines(barcodes_vec, con = barcode_path)

    # call peaks
    gr_list <- CallPeaks(
                    object = unlist(allfragpaths),  ### object fragment(s)
                    macs3.path = macs3.path,
                    mode = mode,
                    outdir = master_outdir,
                    broad = broad,
                    barcodes = barcode_path,
                    genome = genome,
                    gsize = gsize,
                    additional.args = additional.args,
                    name = name,
                    cleanup = cleanup,
                    verbose = verbose)
    
    return(gr_list)
}

#' @method CallPeaks ChromatinAssay
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.ChromatinAssay <- function(
    object,
    macs3.path = NULL,
    mode = "callpeak", # add callpeak/hmmratac command option
    outdir = tempdir(),
    broad = FALSE,
    barcodes = NULL,
    format = "FRAG",
    genome = "hs",
    gsize = NULL,
    additional.args = NULL,
    name = "macs3",
    parallel = TRUE,  ## add parallelization
    cleanup = TRUE,
    verbose = TRUE
) {
    # get fragment files
    frags <- Fragments(object = object)
    # get all fragment file paths
    allfragpaths <- as.list(sapply(X = frags, FUN = GetFragmentData, slot = "file.path"))

    # write cell barcodes
    n_frags <- seq_along(frags)
    barcode_paths <- c()
    
    for (i in seq_along(n_frags)) {
        cell_barcodes <- frags[[i]]@`cells`
        barcode_path <- paste0(outdir, .Platform$file.sep, paste0(i, "_barcodes.txt"))
        writeLines(cell_barcodes, con = barcode_path)
        barcode_paths[[i]] <- barcode_path
    }

    # clean objects
    rm(object)
    rm(frags)
    gc()

    # MACS3 parallelization
    if (parallel) {
        # detect maximum number of cores
        num_fragments <- length(allfragpaths)
        actual_cpus <- parallelly::availableCores()
        plan(future::multisession, workers = min(num_fragments, actual_cpus))
        
        # run CallPeaks in parallel
        gr_list <- future_lapply(n_frags, function(i) {
            CallPeaks(
                object = allfragpaths[i],
                macs3.path = macs3.path,
                mode = mode,
                outdir = outdir,
                broad = broad,
                barcodes = barcode_paths[i],
                genome = genome,
                gsize = gsize,
                additional.args = additional.args,
                name = paste0(i, "_", name),
                cleanup = cleanup,
                verbose = verbose
            )
        }, future.seed = TRUE)
        
    } else {
        # If parallel = FALSE, run CallPeaks sequentially
        gr_list <- lapply(n_frags, function(i) {
            CallPeaks(
                object = allfragpaths[i],
                macs3.path = macs3.path,
                mode = mode,
                outdir = outdir,
                broad = broad,
                barcodes = barcode_paths[i],
                genome = genome,
                gsize = gsize,
                additional.args = additional.args,
                name = paste0(i, "_", name),
                cleanup = cleanup,
                verbose = verbose
            )
        })
    }

    return(gr_list)
}

#' @method CallPeaks Fragment
#' @rdname CallPeaks
#' @concept quantification
#' @export
CallPeaks.Fragment <- function(
    object,
    macs3.path = NULL,
    mode = "callpeak", # add callpeak/hmmratac command option
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
    # get fragment file paths
    fragpath <- sapply(X = object, FUN = GetFragmentData, slot = "file.path") 

    # write cell barcodes
    cell_barcodes <- object[[1]]@`cells`
    barcode_path <- paste0(outdir, .Platform$file.sep, paste0("barcodes.txt"))
    writeLines(cell_barcodes, con = barcode_path)

    # clean objects
    rm(object)
    gc()

    gr <- CallPeaks(
        object = fragpath,
        macs3.path = macs3.path, 
        mode = mode,
        outdir = outdir,
        broad = broad,
        barcodes = barcode_path,
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
    # find macs3
    macs3.path <- macs3.path %||% unname(obj = Sys.which(names = "macs3"))
    if (nchar(x = macs3.path) == 0 && file.access(macs3.path, 1) == 0) {
        stop(
            "MACS3 not found. Please install MACS:",
            "https://macs3-project.github.io/MACS/"
        )
      }
    name <- gsub(pattern = " ", replacement = "_", x = name)
    name <- gsub(pattern = .Platform$file.sep, replacement = "_", x = name)

    # check mode
    if (!mode %in% c("callpeak")) {
        stop("Invalid macs3 command")
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

    # macs3 command
    cmd <- paste0(
        macs3.path, " ",
        mode, 
        " -t ",  object,
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
    files.to.remove <- paste0(outdir, .Platform$file.sep, files.to.remove)
    for (i in files.to.remove) {
      if (file.exists(i)) {
        file.remove(i)
      }
    }
  }
  return(gr)
}

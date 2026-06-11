#' @include generics.R
#'
NULL

#' Export bigwig files for groups of cells
#'
#' Export Tn5 insertion coverage tracks as bigwig files, one per group of cells.
#' Fragments are split by group, the genome is tiled, and the number of fragment
#' ends (Tn5 insertion sites) falling in each tile is counted and (optionally)
#' normalized before writing the bigwig files. Note that each fragment
#' contributes two insertion events (one at each end), so the tracks represent
#' insertion-site coverage rather than fragment pileup.
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param group.by The metadata variable used to group the cells
#' @param idents Identities to include (defined by group.by parameter)
#' @param normMethod Normalization method for the bigwig files. Default `RC`.
#' `RC` will divide the number of insertions in a tile by the total number of
#' fragments in the group. A scaling factor of 10^4 will be applied.
#' `ncells` will divide the number of insertions in a tile by the number of
#' cells in the group. `none` (or `NULL`) will apply no normalization.
#' The name of a metadata column can also be passed, in which case insertions
#' will be divided by the sum of that column over the cells in the group, with a
#' scaling factor of 10^4 applied.
#' @param tileSize The size of the tiles in the bigwig file
#' @param minCells The minimum number of cells in a group for it to be exported.
#' Groups with fewer than `minCells` cells are skipped.
#' @param cutoff The maximum number of insertions for a single cell in a given
#' genomic tile. Counts above this value are capped before summing across cells.
#' Note that cells are identified by their fragment-file barcode, so when an
#' assay contains multiple fragment files this cap is shared between any cells
#' that have the same barcode in different files.
#' @param chromosome A vector of chromosomes to export. If `NULL`, use all
#' chromosomes present in `seqlengths`.
#' @param seqlengths Chromosome lengths used to define the genomic tiles. Can be
#' a named numeric vector of chromosome lengths, or any object with a
#' `seqlengths` method such as a `BSgenome` or [Seqinfo::Seqinfo()] object. If
#' `NULL`, the chromosome lengths stored in the object are used; note that these
#' are frequently unset, in which case `seqlengths` must be supplied.
#' @param outdir Directory to write output files (split bed files and bigwigs).
#' Defaults to the current working directory.
#' @param cleanup Remove the intermediate per-group bed files after writing the
#' bigwig files. Default `TRUE`.
#' @param verbose Display messages
#'
#' @importFrom GenomicRanges GRanges slidingWindows
#' @importFrom IRanges IRanges
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pblapply
#' @importFrom SeuratObject DefaultAssay
#'
#' @export
#' @concept bigwig
#'
#' @return Returns a list of paths to the bigwig files that were created
#'
#' @examples
#' \dontrun{
#' ExportBigwig(object, assay = "peaks")
#' }
ExportBigwig <- function(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  normMethod = "RC",
  tileSize = 100,
  minCells = 5,
  cutoff = NULL,
  chromosome = NULL,
  seqlengths = NULL,
  outdir = getwd(),
  cleanup = TRUE,
  verbose = TRUE
) {
  # Check if output directory exists
  if (!dir.exists(paths = outdir)) {
    dir.create(path = outdir, recursive = TRUE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Please install rtracklayer: ",
      "BiocManager::install('rtracklayer')"
    )
  }
  if (length(x = Fragments(object = object)) == 0) {
    stop("This object does not have Fragments, cannot generate bigwig.")
  }
  assay <- assay %||% DefaultAssay(object = object)
  # a NULL normMethod is equivalent to applying no normalization
  normMethod <- normMethod %||% "none"

  # Resolve chromosome lengths, falling back to the object's seqinfo
  chromLengths <- seqlengths %||% seqlengths(x = object)
  if (!is.numeric(x = chromLengths)) {
    # extract a named length vector from a Seqinfo, BSgenome, etc.
    chromLengths <- Seqinfo::seqlengths(x = chromLengths)
  }
  if (is.null(x = chromLengths) || all(is.na(x = chromLengths))) {
    stop(
      "No chromosome lengths available. Supply them via the `seqlengths` ",
      "argument (a named numeric vector, or a Seqinfo or BSgenome object)."
    )
  }
  # drop chromosomes with unknown length
  chromLengths <- chromLengths[!is.na(x = chromLengths)]

  # Subset to requested chromosomes
  if (!is.null(x = chromosome)) {
    missing.chr <- setdiff(x = chromosome, y = names(x = chromLengths))
    if (length(x = missing.chr) > 0) {
      stop(
        "Requested chromosome(s) not found in seqlengths: ",
        paste(missing.chr, collapse = ", ")
      )
    }
    chromLengths <- chromLengths[chromosome]
  }
  availableChr <- names(x = chromLengths)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(
      start = rep(x = 1, times = length(x = availableChr)),
      end = as.numeric(x = chromLengths)
    )
  )
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  # match the group name sanitization done by SplitFragments so that the bed
  # file names we read back line up with the files that were written
  group.cells <- names(x = obj.groups)
  obj.groups <- gsub(pattern = " ", replacement = "_", x = obj.groups)
  obj.groups <- gsub(
    pattern = .Platform$file.sep, replacement = "_", x = obj.groups
  )
  names(x = obj.groups) <- group.cells
  # true number of cells per group (used for 'ncells' normalization; deriving
  # this from the bed files would undercount cells when fragment-file barcodes
  # collide across multiple fragment files)
  group.counts <- table(obj.groups)
  GroupsNames <- names(x = group.counts[group.counts >= minCells])
  if (length(x = GroupsNames) == 0) {
    warning(
      "No groups contain at least minCells (", minCells, ") cells; ",
      "nothing to export."
    )
    return(list())
  }
  # Check if output files already exist
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(
        sprintf(
          paste0(
            "The group \"%s\" is already present in the destination folder ",
            "and will be overwritten !"
          ),
          x
        )
      )
      file.remove(fn)
    }
  })
  # Splitting fragments file for each ident in group.by
  SplitFragments(
    object = object,
    assay = assay,
    group.by = group.by,
    idents = idents,
    outdir = outdir,
    file.suffix = "",
    append = TRUE,
    buffer_length = 256L,
    verbose = verbose
  )
  # Determine the per-group normalization factor. For a metadata column we sum
  # the values over the cells in each group here (keyed by group name), since
  # the object cell names do not necessarily match the cell barcodes written to
  # the split bed files.
  if (tolower(x = normMethod) %in% c("rc", "ncells", "none")) {
    normBy <- NULL
  } else {
    if (!normMethod %in% colnames(x = object[[]])) {
      stop(
        "normMethod must be one of 'RC', 'ncells', 'none', or the name of a ",
        "metadata column. '", normMethod, "' is not a metadata column."
      )
    }
    md <- object[[normMethod]]
    if (!is.numeric(x = md[[1]])) {
      stop(
        "The '", normMethod, "' metadata column must be numeric to be used ",
        "for normalization."
      )
    }
    normBy <- tapply(
      X = md[names(x = obj.groups), 1],
      INDEX = obj.groups,
      FUN = sum
    )
  }

  if (verbose) {
    message("Creating tiles")
  }
  # Create tiles for each chromosome, from GenomicRanges
  tiles <- unlist(
    x = slidingWindows(x = chromSizes, width = tileSize, step = tileSize)
  )
  if (verbose) {
    message("Creating bigwig files at ", outdir)
  }
  # Run the creation of bigwig for each group of cells
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  covFiles <- mylapply(
    GroupsNames,
    FUN = CreateBWGroup,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    group.counts,
    tileSize,
    normMethod,
    cutoff,
    outdir
  )
  # remove the intermediate split bed files (written for every group, including
  # those below minCells)
  if (cleanup) {
    bedfiles <- file.path(
      outdir, paste0(unique(x = obj.groups), ".bed")
    )
    file.remove(bedfiles[file.exists(bedfiles)])
  }
  return(covFiles)
}

# Helper function for ExportBigwig
#
# @param groupNamei The group to be exported
# @param availableChr Chromosomes to be processed
# @param chromLengths Chromosome lengths
# @param tiles The tiles object
# @param normBy Per-group normalization factor (a named vector keyed by group)
# used when normMethod is a metadata column
# @param nCells Per-group cell counts (a named vector keyed by group) used for
# 'ncells' normalization. If NULL, the number of unique cell barcodes in the
# group's bed file is used.
# @param tileSize The size of the tiles in the bigwig file
# @param normMethod Normalization method for the bigwig files
# 'RC' will divide the number of insertions in a tile by the number of fragments
# in the group. A scaling factor of 10^4 will be applied
# 'ncells' will divide the number of insertions in a tile by the number of cells
# in the group. 'none' will apply no normalization method. A meta.data column
# name can also be passed. A scaling factor of 10^4 will be applied
# @param cutoff The maximum number of insertions for a cell in a given tile
# @param outdir The output directory for bigwig file
#
#' @importFrom GenomicRanges seqnames GRanges
#' @importFrom IRanges coverage
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end
#' @importFrom Matrix sparseMatrix rowSums
CreateBWGroup <- function(
  groupNamei,
  availableChr,
  chromLengths,
  tiles,
  normBy,
  nCells = NULL,
  tileSize,
  normMethod,
  cutoff,
  outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message(
      "Please install rtracklayer. ",
      "http://www.bioconductor.org/packages/rtracklayer/"
    )
    return(NULL)
  }
  normMethod <- tolower(x = normMethod)
  # Read the fragments file associated with the group
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"),
    format = "bed"
  )
  cellGroupi <- unique(x = fragi$name)
  # Open the writing bigwig file
  covFile <- file.path(
    outdir,
    paste0(
      groupNamei, "-TileSize-", tileSize, "-normMethod-", normMethod, ".bw"
    )
  )

  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(x = fragi) == availableChr[k], ]
    tilesk <- tiles[
      BiocGenerics::which(
        S4Vectors::match(seqnames(x = tiles), availableChr[k], nomatch = 0) > 0
      )
    ]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
      # If fragments
    } else {
      # N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      # Add one tile if there are extra bases
      if (nTiles %% 1 != 0) {
        nTiles <- trunc(x = nTiles) + 1
      }
      # Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(x = fragik)$name, cellGroupi)

      # For each tile of this chromosome, create start tile and end tile row,
      # set the associated counts matching with the fragments
      # This changes compared to ArchR version 1.0.2
      # See https://github.com/GreenleafLab/ArchR/issues/2214
      # Clamp tile indices so fragment ends sitting beyond the chromosome
      # boundary (e.g. from read extension) fall in the last tile rather than
      # producing an out-of-bounds matrix index.
      startTile <- pmin(
        trunc(x = (start(x = fragik) - 1) / tileSize) + 1, nTiles
      )
      endTile <- pmin(
        trunc(x = (end(x = fragik) - 1) / tileSize) + 1, nTiles
      )
      mat <- sparseMatrix(
        i = c(startTile, endTile),
        j = as.vector(x = c(matchID, matchID)),
        x = rep(x = 1, times = 2 * length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )

      # Max count for a cell in a tile is set to cutoff
      if (!is.null(x = cutoff)) {
        mat@x[mat@x > cutoff] <- cutoff
      }
      # Sum the cells
      mat <- rowSums(x = mat)
      tilesk$reads <- mat
      # Normalization
      if (normMethod == "rc") {
        tilesk$reads <- tilesk$reads * 10^4 / length(x = fragi$name)
      } else if (normMethod == "ncells") {
        # use the true per-group cell count when supplied, otherwise fall back
        # to the number of unique barcodes in the bed file
        n.cells <- if (is.null(x = nCells)) {
          length(x = cellGroupi)
        } else {
          nCells[[groupNamei]]
        }
        tilesk$reads <- tilesk$reads / n.cells
      } else if (normMethod == "none") {
        # no normalization
      } else {
        # normBy holds the per-group sum of the requested metadata column
        tilesk$reads <- tilesk$reads * 10^4 / normBy[[groupNamei]]
      }
    }
    tilesk <- coverage(x = tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })

  names(x = covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}

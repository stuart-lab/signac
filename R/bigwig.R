#' Export bigwig files for groups of cells
#'
#' @param object A Seurat object
#' @param assay Name of assay to use
#' @param group.by The metadata variable used to group the cells
#' @param idents Identities to include (defined by group.by parameter)
#' @param normMethod Normalization method for the bigwig files. Deafult 'RC'.
#' 'RC' will divide the number of fragments in a tile by the total number of
#' fragments in the group. A scaling factor of 10^4 will be applied. 
#' 'ncells' will divide the number of fragments in a tile by the number of cell
#' in the group. 'none' will apply no normalization method. 
#' A vector of values for each cell can also be passed as a metadata column
#' name. A scaling factor of 10^4 will be applied
#' @param tileSize The size of the tiles in the bigwig file
#' @param minCells The minimum of cells in a group to be exported
#' @param cutoff The maximum number of fragments in a given genomic tile
#' @param chromosome A chromosomes vector to be exported
#' @param outdir Directory to write output files (splitted bed and bigwigs)
#' @param verbose Display messages
#'
#' @importFrom GenomicRanges GRanges slidingWindows
#' @importFrom future plan nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pbapply
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ExportGroupBW(object, assay = "peaks")
#' }
ExportGroupBW  <- function(
    object,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    chromosome = NULL,
    outdir = getwd(),
    verbose = TRUE
) {
  # Check if temporary directory exist
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  if (length(Fragments(object)) == 0) {
    stop("This object does not have Fragments, cannot generate bigwig.")
  }
  # Get chromosome information
  if(!is.null(x = chromosome)){
    seqlevels(object) <- chromosome
  }
  chromLengths <- seqlengths(object)
  if (is.null(chromLengths)) {
    stop("Object has no seqlength, bigwig coverages cannot be evaluated.")
  }
  availableChr <- names(chromLengths)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(
      start = rep(1, length(x = availableChr)),
      end = as.numeric(x = chromLengths)
    )
  )
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  obj.groups <- GetGroups(
    object = object,
    group.by = group.by,
    idents = idents
  )
  GroupsNames <- names(x = table(obj.groups)[table(obj.groups) > minCells])
  # Check if output files already exist
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(sprintf("The group \"%s\" is already present in the destination folder and will be overwritten !",x))
      file.remove(fn)
    }
  })      
  # Splitting fragments file for each idents in group.by
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
  # Column to normalized by
  if(!is.null(x = normMethod)) {
    if (tolower(x = normMethod) %in% c('rc', 'ncells', 'none')){
      normBy <- normMethod
    } else{
      normBy <- object[[normMethod, drop = FALSE]]
    }
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
  # Run the creation of bigwig for each cellgroups
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
    tileSize,
    normMethod,
    cutoff,
    outdir
  )
  return(covFiles)
}

# Helper function for ExportGroupBW
#
# @param groupNamei The group to be exported
# @param availableChr Chromosomes to be processed
# @param chromLengths Chromosomes lengths
# @param tiles The tiles object
# @param normBy A vector of values to normalized the cells by
# @param tileSize The size of the tiles in the bigwig file
# @param normMethod Normalization method for the bigwig files
# 'RC' will divide the number of fragments in a tile by the number of fragments
# in the group. A scaling factor of 10^4 will be applied
# 'ncells' will divide the number of fragments in a tile by the number of cell
# in the group. 'none' will apply no normalization method. A vector of values
# for each cell can also be passed as a meta.data column name. A scaling factor
# of 10^4 will be applied
# @param cutoff The maximum number of fragment in a given tile
# @param outdir The output directory for bigwig file
#
#' @importFrom GenomicRanges seqnames GRanges coverage
#' @importFrom BiocGenerics which
#' @importFrom S4Vectors match
#' @importFrom Matrix sparseMatrix rowSums
CreateBWGroup <- function(
    groupNamei,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  normMethod <- tolower(x = normMethod)
  # Read the fragments file associated to the group
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"),
    format = "bed"
  )
  cellGroupi <- unique(x = fragi$name)
  # Open the writing bigwig file
  covFile <- file.path(
    outdir,
    paste0(groupNamei, "-TileSize-",tileSize,"-normMethod-",normMethod,".bw")
  )
  
  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(fragi) == availableChr[k],]
    tilesk <- tiles[BiocGenerics::which(S4Vectors::match(seqnames(tiles), availableChr[k], nomatch = 0) > 0)]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
      # If fragments
    } else {
      # N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      # Add one tile if there is extra bases
      if (nTiles%%1 != 0) {
        nTiles <- trunc(x = nTiles) + 1
      }
      # Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      
      # For each tiles of this chromosome, create start tile and end tile row,
      # set the associated counts matching with the fragments
      # This changes compared to ArchR version 1.0.2
      # See https://github.com/GreenleafLab/ArchR/issues/2214
      mat <- sparseMatrix(
        i = c(trunc(x = (start(x = fragik) - 1) / tileSize),
              trunc(x = (end(x = fragik) - 1) / tileSize)) + 1,
        j = as.vector(x = c(matchID, matchID)),
        x = rep(1, 2*length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )
      
      # Max count for a cells in a tile is set to cutoff
      if (!is.null(x = cutoff)){
        mat@x[mat@x > cutoff] <- cutoff
      }
      # Sums the cells
      mat <- rowSums(x = mat)
      tilesk$reads <- mat
      # Normalization
      if (!is.null(x = normMethod)) {
        if (normMethod == "rc") {
          tilesk$reads <- tilesk$reads * 10^4 / length(fragi$name)
        } else if (normMethod == "ncells") {
          tilesk$reads <- tilesk$reads / length(cellGroupi)
        } else if (normMethod == "none") {
        } else {
          if (!is.null(x = normBy)){
            tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
          }
        }
      }
    }
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })
  
  names(covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}

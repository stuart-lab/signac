#' @include generics.R
#'
NULL

#' @param assay Name of assay to use. If NULL, use the default assay
#' @param key Key for link information to extract from the assay
#' @importFrom SeuratObject DefaultAssay
#' @method GetLinkedPeaks Seurat
#' @concept links
#' @rdname GetLinkedPeaks
#' @export
GetLinkedPeaks.Seurat <- function(
  object,
  features,
  key,
  assay = NULL,
  min.abs.score = 0,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  links <- GetLinkedPeaks(
    object = object[[assay]],
    features = features,
    key = key,
    min.abs.score = min.abs.score,
    ...
  )
  return(links)
}

#' @export
#' @method GetLinkedPeaks Assay5
#' @concept links
#' @rdname GetLinkedPeaks
GetLinkedPeaks.Assay5 <- function(
  object,
  ...
) {
  stop("GetLinkedPeaks requires a ChromatinAssay5 object")
}

#' @param features A list of genes to find linked peaks for
#' @param min.abs.score Minimum absolute value of the link score for a link to
#' be returned
#' @export
#' @method GetLinkedPeaks ChromatinAssay5
#' @importFrom InteractionSet anchors
#' @concept links
#' @rdname GetLinkedPeaks
GetLinkedPeaks.ChromatinAssay5 <- function(
  object,
  features,
  key,
  min.abs.score = 0,
  ...
) {
  if (missing(x = key)) {
    stop("Provide a key for link information to extract")
  }
  lnk <- Links(object = object)[[key]]
  if (length(x = lnk) == 0) {
    stop("No links present in assay. Run LinkPeaks first.")
  }
  lnk.keep <- lnk[(abs(x = lnk$score) > min.abs.score) &
                    lnk$anchor2.gene_name %in% features]
  return(unique(x = as.character(x = anchors(x = lnk.keep)$first)))
}

#' @param assay Name of assay to use. If NULL, use the default assay
#' @param key Key for link information to extract from the assay
#' @importFrom SeuratObject DefaultAssay
#' @method GetLinkedGenes Seurat
#' @concept links
#' @rdname GetLinkedGenes
#' @export
GetLinkedGenes.Seurat <- function(
  object,
  features,
  key,
  assay = NULL,
  min.abs.score = 0,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  links <- GetLinkedGenes(
    object = object[[assay]],
    features = features,
    key = key,
    min.abs.score = min.abs.score,
    ...
  )
  return(links)
}

#' @export
#' @method GetLinkedGenes Assay5
#' @concept links
#' @rdname GetLinkedGenes
GetLinkedGenes.Assay5 <- function(
  object,
  ...
) {
  stop("GetLinkedGenes requires a ChromatinAssay5 object")
}

#' @param features A list of peaks to find linked genes for
#' @param min.abs.score Minimum absolute value of the link score for a link to
#' be returned
#' @export
#' @method GetLinkedGenes ChromatinAssay5
#' @importFrom InteractionSet anchors
#' @concept links
#' @rdname GetLinkedGenes
GetLinkedGenes.ChromatinAssay5 <- function(
  object,
  features,
  key,
  min.abs.score = 0,
  ...
) {
  if (missing(x = key)) {
    stop("Provide a key for link information to extract")
  }
  lnk <- Links(object = object)[[key]]
  if (length(x = lnk) == 0) {
    stop("No links present in assay. Run LinkPeaks first.")
  }
  pknames <- as.character(x = anchors(x = lnk)$first)
  lnk.keep <- lnk[(abs(x = lnk$score) > min.abs.score) & pknames %in% features]
  return(unique(x = lnk.keep$anchor2.gene_name))
}

#' Cicero connections to links
#'
#' Convert the output of Cicero connections to a
#' `InteractionSet::GInteractions()` object.
#'
#' See the Cicero package for more information:
#' <https://bioconductor.org/packages/cicero/>
#'
#' @param conns A dataframe containing co-accessible elements. This would
#' usually be the output of `run_cicero` or
#' `assemble_connections`. Specifically, this should be a
#' dataframe where the first column contains the genomic coordinates of the
#' first element in the linked pair of elements, with chromosome, start, end
#' coordinates separated by "-" characters. The second column should be the
#' second element in the linked pair, formatted in the same way as the first
#' column. A third column should contain the co-accessibility scores.
#' @param ccans This is optional, but if supplied should be a dataframe
#' containing the cis-co-accessibility network (CCAN) information generated
#' by `generate_ccans`. Specifically, this should be a
#' dataframe containing the name of the peak in the first column, and the
#' CCAN that it belongs to in the second column.
#' @param threshold Threshold for retaining a coaccessible site. Links with
#' a value less than or equal to this threshold will be discarded.
#'
#' @export
#' @importFrom InteractionSet GInteractions
#'
#' @concept links
#' @return Returns a [InteractionSet::GInteractions()] object
ConnectionsToLinks <- function(
  conns,
  ccans = NULL,
  threshold = 0
) {
  # add group information
  if (!is.null(x = ccans)) {
    ccan.lookup <- ccans$CCAN
    names(x = ccan.lookup) <- ccans$Peak
    groups <- as.vector(x = ifelse(
      test = is.na(x = ccan.lookup[conns$Peak1]),
      yes = ccan.lookup[conns$Peak2],
      no = ccan.lookup[conns$Peak1]
    ))
    conns$group <- groups
  } else {
    conns$group <- NA
  }

  # filter based on threshold
  conns <- conns[!is.na(conns$coaccess), ]
  conns <- conns[conns$coaccess > threshold, ]

  # create ginteractions
  gi <- GInteractions(GRanges(conns$Peak1), GRanges(conns$Peak2))
  gi$score <- conns$coaccess
  gi$group <- conns$group

  return(gi)
}



#' Find candidate peak-gene links
#'
#' Identify candidate peak-gene pairs within a specified distance from a gene
#' transcription start site (TSS). Candidates can be matched either to the
#' collapsed most-5-prime gene model or to the nearest transcript TSS for each
#' peak-gene pair.
#'
#' @param object A Seurat object.
#' @param peak.assay Name of assay containing peak information.
#' @param expression.assay Name of assay containing gene expression information.
#' @param peak.layer Name of layer to pull chromatin data from.
#' @param expression.layer Name of layer to pull expression data from.
#' @param gene.coords A [GenomicRanges::GRanges] object containing gene or
#' transcript coordinates. If `NULL`, gene annotations are extracted from the
#' peak assay.
#' @param distance Maximum distance from a TSS for a peak-gene pair to be
#' considered a candidate link.
#' @param min.distance Optional minimum distance from a TSS. Candidate links
#' closer than this distance are excluded.
#' @param min.cells Minimum number of cells positive for a peak or gene for it
#' to be retained.
#' @param genes.use Optional vector of genes to test. If `NULL`, genes are
#' selected from the expression assay after `min.cells` filtering.
#' @param gene.id Set to `TRUE` if genes in the expression assay are named by
#' gene ID rather than gene name.
#' @param verbose Display messages.
#' @param tss.match.strategy Strategy used to assign candidate TSSs. Use
#' `"most5"` to use one collapsed gene-level TSS, or `"transcript_nearest"` to
#' choose the nearest transcript TSS for each peak-gene pair.
#'
#' @return A list containing objects used downstream by [LinkPeaks()]:
#' `candidate.matrix`, `candidate.links`, `peaks`, `gene.coords`,
#' `candidate.link.table`, `peak.data`, and `expression.data`.
#'
#' @importFrom GenomicRanges GRanges granges seqnames start end width strand resize findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom InteractionSet GInteractions
#' @importFrom Matrix rowSums sparseMatrix drop0 summary
#' @importFrom SeuratObject LayerData Layers
#' @importFrom stats ave
#'
#' @export
#' @concept links
FindCandidateLinks <- function(
  object,
  peak.assay,
  expression.assay,
  peak.layer = "counts",
  expression.layer = "data",
  gene.coords = NULL,
  distance = 5e5,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  gene.id = FALSE,
  verbose = TRUE,
  tss.match.strategy = c("transcript_nearest", "most5")
) {
  tss.match.strategy <- match.arg(tss.match.strategy)
  gene.key <- if (gene.id) "gene_id" else "gene_name"

  .peak_to_tss_distance <- function(peaks, tss, peak_index, tss_index) {
    peak.starts <- start(peaks)[peak_index]
    peak.ends <- end(peaks)[peak_index]
    tss.pos <- start(tss)[tss_index]

    ifelse(
      tss.pos >= peak.starts & tss.pos <= peak.ends,
      0,
      pmin(abs(tss.pos - peak.starts), abs(tss.pos - peak.ends))
    )
  }

  .add_closest_gene_info <- function(df) {
    df$closest_tss <- as.integer(
      df$tss_distance == ave(df$tss_distance, df$peak_index, FUN = min)
    )
    df$closest_gene <- df$closest_tss
    df$closest_tss_id <- NA_character_

    if ("transcript_id" %in% colnames(df)) {
      by.peak <- split(seq_len(nrow(df)), df$peak_index)

      for (idx in by.peak) {
        closest.idx <- idx[df$closest_tss[idx] == 1L]
        df$closest_tss_id[idx] <- paste(
          unique(as.character(df$transcript_id[closest.idx])),
          collapse = ";"
        )
      }
    }

    df
  }

  .add_peak_gene_overlap_info <- function(df, peaks, gene.body.coords, gene.key) {
    df$peak_overlap_gene <- 0L
    df$frac_peak_in_gene <- 0

    if (nrow(df) == 0 || length(gene.body.coords) == 0) {
      return(df)
    }

    gene.body.values <- as.character(
      mcols(gene.body.coords)[[gene.key]]
    )

    body.index <- match(as.character(df$gene), gene.body.values)
    ok <- !is.na(body.index)

    if (any(ok)) {
      peak.ranges <- peaks[df$peak_index[ok]]
      body.ranges <- gene.body.coords[body.index[ok]]

      same.seq <- as.character(seqnames(peak.ranges)) ==
        as.character(seqnames(body.ranges))

      overlap.start <- pmax(
        start(peak.ranges),
        start(body.ranges)
      )
      overlap.end <- pmin(
        end(peak.ranges),
        end(body.ranges)
      )

      overlap.width <- ifelse(
        same.seq,
        pmax(0L, overlap.end - overlap.start + 1L),
        0L
      )

      df$peak_overlap_gene[ok] <- as.integer(overlap.width > 0L)
      df$frac_peak_in_gene[ok] <- overlap.width / width(peak.ranges)
    }

    df
  }

  if (!inherits(object[[peak.assay]], "GRangesAssay")) {
    stop("The requested peak assay is not a GRangesAssay")
  }

  if (!is.null(min.distance)) {
    if (!is.numeric(min.distance)) stop("min.distance should be numeric")
    if (min.distance <= 0) min.distance <- NULL
    if (!is.null(min.distance) && min.distance >= distance) {
      stop("min.distance should be smaller than distance")
    }
  }

  gene.coords.input <- gene.coords
  if (is.null(gene.coords.input)) {
    gene.coords.input <- Annotation(object[[peak.assay]])
    if (is.null(gene.coords.input)) stop("Gene annotations not found")
  }

  if (!("gene_id" %in% colnames(mcols(gene.coords.input)))) {
    stop("gene.coords / annotation must contain 'gene_id'")
  }

  if (!("gene_name" %in% colnames(mcols(gene.coords.input)))) {
    if (gene.id) {
      mcols(gene.coords.input)$gene_name <-
        mcols(gene.coords.input)$gene_id
    } else {
      stop("gene.coords / annotation must contain 'gene_name' when gene.id = FALSE")
    }
  }

  if (!(peak.layer %in% Layers(object[[peak.assay]]))) {
    stop("Requested peak layer not found")
  }

  if (!(expression.layer %in% Layers(object[[expression.assay]]))) {
    stop("Requested expression layer not found")
  }

  peak.data <- LayerData(object[[peak.assay]], layer = peak.layer)
  expression.data <- LayerData(object[[expression.assay]], layer = expression.layer)

  peaks.keep <- rowSums(peak.data > 0) > min.cells
  genes.keep <- names(which(rowSums(expression.data > 0) > min.cells))

  if (!is.null(genes.use)) {
    genes.keep <- intersect(genes.keep, genes.use)
  }

  peak.data <- peak.data[peaks.keep, , drop = FALSE]
  expression.data <- expression.data[genes.keep, , drop = FALSE]

  if (nrow(peak.data) == 0) stop("No peaks pass min.cells filtering")
  if (nrow(expression.data) == 0) stop("No genes pass min.cells / genes.use filtering")

  genes <- rownames(expression.data)
  peaks <- granges(object[[peak.assay]])[peaks.keep]
  peak.names <- rownames(peak.data)
  if (is.null(peak.names)) peak.names <- as.character(peaks)

  gene.body.coords <- CollapseToLongestTranscript(gene.coords.input)

  if (gene.id) {
    gene.body.coords <- gene.body.coords[gene.body.coords$gene_id %in% genes]
    gene.body.coords$gene_name <- gene.body.coords$gene_id
  } else {
    gene.body.coords <- gene.body.coords[gene.body.coords$gene_name %in% genes]
  }

  gene.body.values <- as.character(
    mcols(gene.body.coords)[[gene.key]]
  )
  keep.body.genes <- genes[genes %in% gene.body.values]
  gene.body.coords <- gene.body.coords[match(keep.body.genes, gene.body.values)]

  if (tss.match.strategy == "most5") {
    gene.coords.use <- gene.body.coords
    gene.values <- mcols(gene.coords.use)[[gene.key]]

    keep.genes <- genes[genes %in% gene.values]
    if (length(keep.genes) == 0) {
      stop("No expressed genes matched gene coordinates")
    }

    gene.coords.use <- gene.coords.use[match(keep.genes, gene.values)]
    gene.values <- mcols(gene.coords.use)[[gene.key]]

    candidate.matrix <- DistanceToTSS(
      peaks = peaks,
      genes = gene.coords.use,
      distance = distance
    )

    rownames(candidate.matrix) <- peak.names
    colnames(candidate.matrix) <- gene.values

    if (!is.null(min.distance)) {
      inner.matrix <- DistanceToTSS(
        peaks = peaks,
        genes = gene.coords.use,
        distance = min.distance
      )

      rownames(inner.matrix) <- rownames(candidate.matrix)
      colnames(inner.matrix) <- colnames(candidate.matrix)

      candidate.matrix <- candidate.matrix - inner.matrix
      candidate.matrix@x[candidate.matrix@x < 0] <- 0
      candidate.matrix <- drop0(candidate.matrix)
    }

    if (sum(candidate.matrix) == 0) {
      stop("No candidate peak-gene links found")
    }

    candidate.tss <- resize(
      gene.coords.use,
      width = 1,
      fix = "start"
    )
    mcols(candidate.tss) <- NULL

    candidate.link.table <- as.data.frame(summary(candidate.matrix))
    colnames(candidate.link.table) <- c("peak_index", "gene_index", "x")

    candidate.link.table$peak <- rownames(candidate.matrix)[candidate.link.table$peak_index]
    candidate.link.table$gene <- gene.values[candidate.link.table$gene_index]
    candidate.link.table$transcript_id <- NA_character_

    candidate.link.table$tss_distance <- .peak_to_tss_distance(
      peaks = peaks,
      tss = candidate.tss,
      peak_index = candidate.link.table$peak_index,
      tss_index = candidate.link.table$gene_index
    )

    candidate.link.table <- .add_closest_gene_info(candidate.link.table)
    candidate.link.table <- .add_peak_gene_overlap_info(
      df = candidate.link.table,
      peaks = peaks,
      gene.body.coords = gene.body.coords,
      gene.key = gene.key
    )

    candidate.links <- GInteractions(
      peaks[candidate.link.table$peak_index],
      candidate.tss[candidate.link.table$gene_index]
    )

    candidate.links$peak <- candidate.link.table$peak
    candidate.links$gene <- candidate.link.table$gene
    candidate.links$transcript_id <- candidate.link.table$transcript_id
    candidate.links$tss_distance <- candidate.link.table$tss_distance
    candidate.links$peak_overlap_gene <- candidate.link.table$peak_overlap_gene
    candidate.links$frac_peak_in_gene <- candidate.link.table$frac_peak_in_gene
    candidate.links$closest_gene <- candidate.link.table$closest_gene
    candidate.links$closest_tss_id <- candidate.link.table$closest_tss_id
  }

  if (tss.match.strategy == "transcript_nearest") {
    transcript.coords.use <- gene.coords.input

    if ("type" %in% colnames(mcols(transcript.coords.use))) {
      tx.rows <- transcript.coords.use$type %in% c("transcript", "mRNA")
      if (any(tx.rows, na.rm = TRUE)) {
        transcript.coords.use <- transcript.coords.use[tx.rows]
      }
    }

    if (gene.id) {
      transcript.coords.use <- transcript.coords.use[
        transcript.coords.use$gene_id %in% genes
      ]
      transcript.coords.use$gene_name <- transcript.coords.use$gene_id
    } else {
      transcript.coords.use <- transcript.coords.use[
        transcript.coords.use$gene_name %in% genes
      ]
    }

    if (length(transcript.coords.use) == 0) {
      stop("No transcript/gene coordinates found")
    }

    transcript.gene.values <- mcols(transcript.coords.use)[[gene.key]]
    gene.values <- genes[genes %in% unique(transcript.gene.values)]

    if (length(gene.values) == 0) {
      stop("No expressed genes matched transcript coordinates")
    }

    gene.coords.use <- transcript.coords.use[match(gene.values, transcript.gene.values)]

    transcript.tss <- resize(
      transcript.coords.use,
      width = 1,
      fix = "start"
    )

    outer.hits <- findOverlaps(
      query = peaks,
      subject = Extend(
        transcript.tss,
        upstream = distance,
        downstream = distance
      ),
      type = "any"
    )

    if (length(outer.hits) == 0) {
      stop("No candidate peak-transcript links found")
    }

    outer.df <- data.frame(
      peak_index = queryHits(outer.hits),
      transcript_index = subjectHits(outer.hits)
    )

    outer.df$gene <- transcript.gene.values[outer.df$transcript_index]
    outer.df <- outer.df[outer.df$gene %in% gene.values, , drop = FALSE]

    if (nrow(outer.df) == 0) {
      stop("No candidate peak-transcript links matched expressed genes")
    }

    if (!is.null(min.distance)) {
      inner.hits <- findOverlaps(
        query = peaks,
        subject = Extend(
          transcript.tss,
          upstream = min.distance,
          downstream = min.distance
        ),
        type = "any"
      )

      if (length(inner.hits) > 0) {
        inner.df <- data.frame(
          peak_index = queryHits(inner.hits),
          transcript_index = subjectHits(inner.hits)
        )

        inner.df$gene <- transcript.gene.values[inner.df$transcript_index]
        inner.df <- inner.df[inner.df$gene %in% gene.values, , drop = FALSE]

        outer.df <- outer.df[
          !(paste(outer.df$peak_index, outer.df$gene, sep = "\r") %in%
              paste(inner.df$peak_index, inner.df$gene, sep = "\r")),
          ,
          drop = FALSE
        ]
      }
    }

    if (nrow(outer.df) == 0) {
      stop("No candidate links remain after min.distance filtering")
    }

    outer.df$gene_index <- match(outer.df$gene, gene.values)

    outer.df$tss_distance <- .peak_to_tss_distance(
      peaks = peaks,
      tss = transcript.tss,
      peak_index = outer.df$peak_index,
      tss_index = outer.df$transcript_index
    )

    outer.df$tss_seqname <- as.character(
      seqnames(transcript.tss)
    )[outer.df$transcript_index]

    outer.df$tss_position <- start(transcript.tss)[
      outer.df$transcript_index
    ]

    outer.df$tss_strand <- as.character(
      strand(transcript.tss)
    )[outer.df$transcript_index]

    outer.df$transcript_id <- if (
      "transcript_id" %in% colnames(mcols(transcript.coords.use))
    ) {
      as.character(transcript.coords.use$transcript_id[outer.df$transcript_index])
    } else {
      paste0("transcript_row_", outer.df$transcript_index)
    }

    outer.df <- outer.df[order(
      outer.df$peak_index,
      outer.df$gene_index,
      outer.df$tss_distance,
      outer.df$transcript_index
    ), , drop = FALSE]

    outer.df <- outer.df[
      !duplicated(paste(outer.df$peak_index, outer.df$gene_index, sep = "\r")),
      ,
      drop = FALSE
    ]

    outer.df$peak <- peak.names[outer.df$peak_index]
    outer.df <- .add_closest_gene_info(outer.df)
    outer.df <- .add_peak_gene_overlap_info(
      df = outer.df,
      peaks = peaks,
      gene.body.coords = gene.body.coords,
      gene.key = gene.key
    )

    candidate.matrix <- sparseMatrix(
      i = outer.df$peak_index,
      j = outer.df$gene_index,
      x = 1,
      dims = c(length(peaks), length(gene.values))
    )

    rownames(candidate.matrix) <- peak.names
    colnames(candidate.matrix) <- gene.values

    candidate.tss <- GRanges(
      seqnames = outer.df$tss_seqname,
      ranges = IRanges(
        start = outer.df$tss_position,
        end = outer.df$tss_position
      ),
      strand = outer.df$tss_strand
    )
    mcols(candidate.tss) <- NULL

    candidate.links <- GInteractions(
      peaks[outer.df$peak_index],
      candidate.tss
    )

    candidate.links$peak <- outer.df$peak
    candidate.links$gene <- outer.df$gene
    candidate.links$transcript_id <- outer.df$transcript_id
    candidate.links$tss_distance <- outer.df$tss_distance
    candidate.links$peak_overlap_gene <- outer.df$peak_overlap_gene
    candidate.links$frac_peak_in_gene <- outer.df$frac_peak_in_gene
    candidate.links$closest_gene <- outer.df$closest_gene
    candidate.links$closest_tss_id <- outer.df$closest_tss_id

    candidate.link.table <- outer.df
  }

  if (verbose) {
    message(
      "Identified ", length(candidate.links),
      " candidate peak-gene links across ", ncol(candidate.matrix),
      " genes and ", sum(rowSums(candidate.matrix) > 0),
      " peaks using strategy '", tss.match.strategy, "'"
    )
  }

  candidate.link.table <- candidate.link.table[, c("peak", "gene"), drop = FALSE]

  # Return only objects used downstream by LinkPeaks().
  list(
    candidate.matrix = candidate.matrix,
    candidate.links = candidate.links,
    peaks = peaks,
    gene.coords = gene.coords.use,
    candidate.link.table = candidate.link.table,
    peak.data = peak.data,
    expression.data = expression.data
  )
}

#' Link peaks to genes
#'
#' Find peaks that are correlated with the expression of nearby genes. For each
#' gene, this function computes the correlation coefficient between gene
#' expression and accessibility of candidate peaks within a specified distance
#' from the gene TSS. By default, links are retained using the observed
#' correlation coefficient only. If `calculate.zscore = TRUE`, the function also
#' computes an expected correlation coefficient for each peak using background
#' peaks matched on GC content, accessibility, and sequence length, and uses this
#' background distribution to compute a z-score and p-value.
#'
#' This function was inspired by the method originally described by SHARE-seq
#' (Sai Ma et al. 2020, Cell). Please consider citing the original SHARE-seq
#' work if using this function:
#' [doi: 10.1016/j.cell.2020.09.056](https://pubmed.ncbi.nlm.nih.gov/33098772/)
#'
#' @param object A Seurat object.
#' @param peak.assay Name of assay containing peak information.
#' @param expression.assay Name of assay containing gene expression information.
#' @param peak.layer Name of layer to pull chromatin data from.
#' @param expression.layer Name of layer to pull expression data from.
#' @param method Correlation method to use. One of `"pearson"` or `"spearman"`.
#' @param key Key to use when storing link information in the assay.
#' @param gene.coords A [GenomicRanges::GRanges] object containing gene or
#' transcript coordinates. If `NULL`, gene annotations are extracted from the
#' peak assay.
#' @param distance Maximum distance from a TSS for peaks to include as candidate
#' links.
#' @param min.distance Optional minimum distance from a TSS. Candidate links
#' closer than this distance are excluded.
#' @param min.cells Minimum number of cells positive for the peak and gene
#' needed to include them in the analysis.
#' @param genes.use Optional vector of genes to test. If `NULL`, genes are
#' selected from the expression assay after `min.cells` filtering.
#' @param n_sample Number of background peaks to sample when computing the null
#' distribution. Only used when `calculate.zscore = TRUE`.
#' @param pvalue_cutoff Maximum p-value for retaining a link. Links with a
#' p-value greater than or equal to this value are removed from the output. Only
#' used when `calculate.zscore = TRUE`.
#' @param score_cutoff Minimum absolute correlation coefficient for a link to be
#' retained.
#' @param calculate.zscore Compute background-matched z-scores and p-values. If
#' `FALSE` (default), background peak matching is skipped and links are retained
#' using `score_cutoff` only.
#' @param gene.id Set to `TRUE` if genes in the expression assay are named by
#' gene ID rather than gene name.
#' @param verbose Display messages.
#' @param tss.match.strategy Strategy used to assign candidate TSSs. Use
#' `"most5"` to use one collapsed gene-level TSS, or `"transcript_nearest"` to
#' choose the nearest transcript TSS for each peak-gene pair.
#' @param peak.slot Deprecated; use `peak.layer`.
#' @param expression.slot Deprecated; use `expression.layer`.
#'
#' @return Returns a Seurat object with results added to the links slot in the
#' peak assay, stored under `key`. The results are stored as an
#' [InteractionSet::GInteractions] object accessible via [Links()]. The metadata
#' stored on each link includes:
#' * `peak`: peak identifier
#' * `gene`: linked gene identifier
#' * `transcript_id`: selected transcript identifier, when available
#' * `tss_distance`: distance from the peak to the selected TSS
#' * `peak_overlap_gene`: binary indicator for any peak overlap with the linked
#' gene body
#' * `frac_peak_in_gene`: fraction of the peak width overlapping the linked gene
#' body
#' * `closest_gene`: binary indicator for whether this gene is closest to the
#' peak among candidate genes
#' * `closest_tss_id`: transcript IDs for the closest TSSs to the peak
#' * `score`: observed correlation coefficient
#' * `zscore`: z-score of the observed correlation coefficient; `NA` when
#' `calculate.zscore = FALSE`
#' * `pvalue`: p-value associated with the z-score; `NA` when
#' `calculate.zscore = FALSE`
#'
#' @importFrom SeuratObject LayerData Layers as.sparse
#' @importFrom Matrix sparseMatrix rowSums drop0
#' @importMethodsFrom Matrix t
#' @importFrom GenomicRanges seqnames
#' @importFrom S4Vectors mcols DataFrame
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pblapply
#' @importFrom lifecycle is_present deprecated deprecate_warn
#' @importFrom stats pnorm sd
#'
#' @export
#' @concept links
LinkPeaks <- function(
  object,
  peak.assay,
  expression.assay,
  peak.layer = "counts",
  expression.layer = "data",
  method = "pearson",
  key = "linkpeaks",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  calculate.zscore = FALSE,
  gene.id = FALSE,
  verbose = TRUE,
  tss.match.strategy = c("transcript_nearest", "most5"),
  peak.slot = deprecated(),
  expression.slot = deprecated()
) {
  tss.match.strategy <- match.arg(tss.match.strategy)

  if (!inherits(x = object[[peak.assay]], what = "GRangesAssay")) {
    stop("The requested assay is not a GRangesAssay")
  }

  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("min.distance should be a numeric value")
    }
    if (min.distance < 0) {
      warning("Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    } else if (min.distance == 0) {
      min.distance <- NULL
    }
  }

  if (is_present(arg = expression.slot)) {
    deprecate_warn(
      when = "1.16.0",
      what = "LinkPeaks(expression.slot)",
      with = "LinkPeaks(expression.layer)"
    )
    expression.layer <- expression.slot
  }

  if (is_present(arg = peak.slot)) {
    deprecate_warn(
      when = "1.16.0",
      what = "LinkPeaks(peak.slot)",
      with = "LinkPeaks(peak.layer)"
    )
    peak.layer <- peak.slot
  }

  features.match <- c("GC.percent", "count", "sequence.length")

  if (method == "pearson") {
    cor_method <- corSparse
  } else if (method == "spearman") {
    cor_method <- SparseSpearmanCor
  } else {
    stop("method can be one of 'pearson' or 'spearman'.")
  }

  meta.features <- NULL
  if (isTRUE(calculate.zscore)) {
    meta.features <- object[[peak.assay]][[]]
    if (!(all(c("GC.percent", "sequence.length") %in% colnames(x = meta.features)))) {
      stop(
        "DNA sequence information for each peak has not been computed.\n",
        "Run RegionStats before calling this function."
      )
    }

    if (!("count" %in% colnames(x = meta.features))) {
      hvf.info <- FindTopFeatures(
        object = LayerData(object = object[[peak.assay]], layer = peak.layer),
        verbose = FALSE
      )
      hvf.info <- hvf.info[
        rownames(meta.features), c("count", "percentile"), drop = FALSE
      ]
      meta.features <- cbind(meta.features, hvf.info)
    }
  }

  candidates <- tryCatch(
    FindCandidateLinks(
      object = object,
      peak.assay = peak.assay,
      expression.assay = expression.assay,
      peak.layer = peak.layer,
      expression.layer = expression.layer,
      gene.coords = gene.coords,
      distance = distance,
      min.distance = min.distance,
      min.cells = min.cells,
      genes.use = genes.use,
      gene.id = gene.id,
      verbose = verbose,
      tss.match.strategy = tss.match.strategy
    ),
    error = function(e) {
      if (verbose) message("Candidate-finding skipped: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(candidates)) {
    Links(object = object[[peak.assay]], key = key) <- make_empty_links()
    return(object)
  }

  peak_distance_matrix <- candidates$candidate.matrix
  candidate.links <- candidates$candidate.links
  candidate.link.table <- candidates$candidate.link.table
  peak.data <- candidates$peak.data
  expression.data <- candidates$expression.data
  gene.coords.use <- candidates$gene.coords
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)
  all.peak.chroms <- if (isTRUE(calculate.zscore)) {
    as.character(seqnames(candidates$peaks))
  } else {
    NULL
  }

  if (sum(peak_distance_matrix) == 0) {
    if (verbose) message("No peaks fall within distance threshold.")
    Links(object = object[[peak.assay]], key = key) <- make_empty_links()
    return(object)
  }

  if (verbose) {
    message(
      "Testing ", nrow(x = expression.data),
      " genes and ", sum(rowSums(x = peak_distance_matrix) > 0),
      " peaks"
    )
  }

  peak.data <- t(x = peak.data)

  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- if (verbose) pblapply else lapply
  }

  res <- mylapply(
    X = seq_along(along.with = genes.use),
    FUN = function(i) {
      peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
      gene.expression <- t(x = expression.data[genes.use[[i]], , drop = FALSE])

      if (sum(peak.use) < 2) {
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      }

      peak.access <- peak.data[, peak.use, drop = FALSE]
      if (inherits(x = peak.access, what = "IterableMatrix")) {
        peak.access <- as.sparse(x = peak.access)
      }
      if (inherits(x = gene.expression, what = "IterableMatrix")) {
        gene.expression <- as.sparse(x = gene.expression)
      }

      coef.result <- cor_method(
        X = peak.access,
        Y = gene.expression
      )
      rownames(x = coef.result) <- colnames(x = peak.access)
      coef.result <- coef.result[
        abs(x = coef.result) > score_cutoff,
        ,
        drop = FALSE
      ]

      if (nrow(x = coef.result) == 0) {
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      }

      peaks.test <- rownames(x = coef.result)
      coef.vals <- as.vector(x = coef.result)
      names(x = coef.vals) <- peaks.test

      if (!isTRUE(x = calculate.zscore)) {
        return(list(
          "gene" = rep(x = i, length(x = coef.vals)),
          "coef" = coef.vals,
          "zscore" = rep(NA_real_, length(x = coef.vals))
        ))
      }

      # select peaks at random with matching GC content and accessibility,
      # sampling from peaks on a different chromosome to the gene. Use the
      # actual peak chromosomes rather than assuming a chr-start-end peak name
      gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
      trans.peaks <- all.peaks[all.peak.chroms != gene.chrom]
      if (length(x = trans.peaks) < 2) {
        # no background available on another chromosome; skip rather than
        # computing a z-score against peaks that are cis to the gene
        warning(
          "Fewer than two peaks are available on a chromosome other than ",
          gene.chrom, "; skipping genes on this chromosome",
          call. = FALSE
        )
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      }
      meta.use <- meta.features[trans.peaks, , drop = FALSE]
      pk.use <- meta.features[peaks.test, , drop = FALSE]

      bg.peaks <- lapply(
        X = seq_len(length.out = nrow(x = pk.use)),
        FUN = function(x) {
          MatchRegionStats(
            meta.feature = meta.use,
            query.feature = pk.use[x, , drop = FALSE],
            features.match = features.match,
            n = n_sample,
            verbose = FALSE
          )
        }
      )

      # run background correlations
      unique.bg <- unique(x = unlist(x = bg.peaks))
      if (length(x = unique.bg) == 0) {
        return(list("gene" = NULL, "coef" = NULL, "zscore" = NULL))
      }

      bg.access <- peak.data[, unique.bg, drop = FALSE]
      if (inherits(x = bg.access, what = "IterableMatrix")) {
        bg.access <- as.sparse(x = bg.access)
      }

      bg.coef <- cor_method(
        X = bg.access,
        Y = gene.expression
      )
      rownames(x = bg.coef) <- unique.bg

      zscores <- vector(mode = "numeric", length = length(x = peaks.test))
      for (j in seq_along(along.with = peaks.test)) {
        coef.use <- bg.coef[bg.peaks[[j]], , drop = FALSE]
        bg.sd <- sd(x = coef.use)
        if (bg.sd == 0 || !is.finite(x = bg.sd)) {
          zscores[[j]] <- 0
        } else {
          zscores[[j]] <- (coef.result[j] - mean(x = coef.use)) / bg.sd
        }
      }
      names(x = zscores) <- peaks.test
      gc(verbose = FALSE)

      # p-values are filtered once on the assembled links below, so no
      # filtering is applied here
      list(
        "gene" = rep(x = i, length(x = coef.vals)),
        "coef" = coef.vals,
        "zscore" = zscores
      )
    }
  )

  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 3))

  if (length(x = coef.vec) == 0) {
    if (verbose) message("No links pass score_cutoff")
    Links(object = object[[peak.assay]], key = key) <- make_empty_links()
    return(object)
  }

  sig.df <- data.frame(
    peak = names(x = coef.vec),
    gene_index = as.integer(gene.vec),
    gene = genes.use[as.integer(gene.vec)],
    score = as.numeric(coef.vec),
    zscore = if (isTRUE(calculate.zscore)) {
      as.numeric(zscore.vec)
    } else {
      rep(NA_real_, length(x = coef.vec))
    },
    stringsAsFactors = FALSE
  )

  if (isTRUE(calculate.zscore)) {
    sig.df$pvalue <- 2 * pnorm(q = -abs(x = sig.df$zscore))
    sig.df <- sig.df[sig.df$pvalue < pvalue_cutoff, , drop = FALSE]

    if (nrow(sig.df) == 0) {
      if (verbose) message("No significant links after p-value filtering")
      Links(object = object[[peak.assay]], key = key) <- make_empty_links()
      return(object)
    }
  } else {
    sig.df$pvalue <- NA_real_
  }

  candidate.key <- paste(candidate.link.table$peak, candidate.link.table$gene, sep = "\r")
  sig.key <- paste(sig.df$peak, sig.df$gene, sep = "\r")
  link.idx <- match(sig.key, candidate.key)

  if (anyNA(link.idx)) {
    missing <- unique(sig.key[is.na(link.idx)])
    stop(
      "Internal error: significant links could not be matched back to candidate links. ",
      "First missing key: ", missing[[1]]
    )
  }

  links <- candidate.links[link.idx]
  links$peak <- sig.df$peak
  links$gene <- sig.df$gene

  if (!("transcript_id" %in% colnames(mcols(links)))) {
    links$transcript_id <- rep(NA_character_, length(links))
  }
  if (!("closest_tss_id" %in% colnames(mcols(links)))) {
    links$closest_tss_id <- rep(NA_character_, length(links))
  }
  if (!("peak_overlap_gene" %in% colnames(mcols(links)))) {
    links$peak_overlap_gene <- rep(NA_integer_, length(links))
  }
  if (!("frac_peak_in_gene" %in% colnames(mcols(links)))) {
    links$frac_peak_in_gene <- rep(NA_real_, length(links))
  }
  if (!("closest_gene" %in% colnames(mcols(links)))) {
    if ("closest_tss" %in% colnames(mcols(links))) {
      links$closest_gene <- as.integer(links$closest_tss)
    } else {
      links$closest_gene <- rep(NA_integer_, length(links))
    }
  }

  links$score <- sig.df$score
  links$zscore <- sig.df$zscore
  links$pvalue <- sig.df$pvalue

  .link_mcol <- function(x, col) {
    value <- mcols(x)[[col]]

    if (length(value) != length(x)) {
      stop(
        "Internal error: links$", col, " has length ", length(value),
        " but length(links) is ", length(x), "."
      )
    }

    value
  }

  # Keep only fields required downstream for export and saved link objects.
  mcols(links) <- DataFrame(
    peak = as.character(.link_mcol(links, "peak")),
    gene = as.character(.link_mcol(links, "gene")),
    transcript_id = as.character(.link_mcol(links, "transcript_id")),
    tss_distance = as.numeric(.link_mcol(links, "tss_distance")),
    peak_overlap_gene = as.integer(.link_mcol(links, "peak_overlap_gene")),
    frac_peak_in_gene = as.numeric(.link_mcol(links, "frac_peak_in_gene")),
    closest_gene = as.integer(.link_mcol(links, "closest_gene")),
    closest_tss_id = as.character(.link_mcol(links, "closest_tss_id")),
    score = as.numeric(.link_mcol(links, "score")),
    zscore = as.numeric(.link_mcol(links, "zscore")),
    pvalue = as.numeric(.link_mcol(links, "pvalue"))
  )

  Links(object = object[[peak.assay]], key = key) <- links
  object
}

### Not exported ###

# Create an empty set of links
#
# Returns a zero-length GInteractions object carrying the same metadata columns
# as the links created by LinkPeaks(), so that downstream code can rely on the
# same schema whether or not any links were found.
#
# @return Returns a GInteractions object of length zero
#' @importFrom GenomicRanges GRanges
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors mcols DataFrame
make_empty_links <- function() {
  empty.gr <- GRanges()
  gi <- GInteractions(anchor1 = empty.gr, anchor2 = empty.gr)
  mcols(x = gi) <- DataFrame(
    peak = character(),
    gene = character(),
    transcript_id = character(),
    tss_distance = numeric(),
    peak_overlap_gene = integer(),
    frac_peak_in_gene = numeric(),
    closest_gene = integer(),
    closest_tss_id = character(),
    score = numeric(),
    zscore = numeric(),
    pvalue = numeric()
  )
  return(gi)
}

#' @importFrom GenomicRanges GRanges
#' @importFrom InteractionSet GInteractions
LinksToGInteractions <- function(linkmat, gene.coords) {
  x <- as(object = linkmat, Class = "TsparseMatrix")
  peak.coords <- GRanges(colnames(x = linkmat))
  region1 <- peak.coords[x@j + 1]
  region2 <- gene.coords[x@i + 1]
  gi <- GInteractions(region1, region2)
  gi$score <- x@x
  return(gi)
}

# Find peaks near genes
#
# Find peaks that are within a given distance threshold to each gene
#
# @param peaks A GRanges object containing peak coordinates
# @param genes A GRanges object containing gene coordinates
# @param distance Distance threshold. Peaks within this distance from the gene
# will be recorded.
#
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges resize
#
# @return Returns a sparse matrix
DistanceToTSS <- function(
  peaks,
  genes,
  distance = 200000
) {
  tss <- resize(x = genes, width = 1, fix = "start")
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = "any",
    select = "all"
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- as.character(x = peaks)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

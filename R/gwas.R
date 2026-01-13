# R/gwas.R
# GWAS Visualization for Signac

#' @importFrom utils globalVariables
#'
NULL
globalVariables(names = c("in_credset", "ld_category", "log10p", "pos"), package = "Signac")

#' Load GWAS summary statistics for a genomic region
#'
#' @param gwas.file Path to GWAS/QTL summary statistics file (tabix-indexed TSV)
#' @param region Genomic region (GRanges object or string like "chr1:1000000-2000000")
#' @return data.frame with columns: chr, pos, pval, snp, and optionally beta, se
#'
#' @details
#' Input format based on GWAS-SSF (NHGRI-EBI GWAS Catalog standard).
#'
#' Required columns: \code{chromosome}, \code{base_pair_location}, \code{p_value}
#'
#' Optional columns: \code{variant_id}, \code{beta}, \code{standard_error}
#'
#' Column matching is case-insensitive. Also accepts QTL data in the same format.
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @export
LoadGWASRegion <- function(gwas.file, region) {
  
  # Convert region to GRanges if string
  if (!inherits(region, "GRanges")) {
    region_str <- gsub(":", "-", region)
    parts <- strsplit(region_str, "-")[[1]]
    if (length(parts) != 3) {
      stop("Region must be in format 'chr:start-end' or 'chr-start-end'")
    }
    region <- GRanges(
      seqnames = parts[1],
      ranges = IRanges(start = as.numeric(parts[2]), end = as.numeric(parts[3]))
    )
  }
  
  # Read file
  gwas_data <- fread(gwas.file, data.table = FALSE)
  colnames_lower <- tolower(colnames(gwas_data))
  
  # Validate required columns
  required <- c("chromosome", "base_pair_location", "p_value")
  missing <- required[!required %in% colnames_lower]
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "), "\n",
      "Found: ", paste(colnames(gwas_data), collapse = ", "), "\n",
      "Format: SSF-based (chromosome, base_pair_location, p_value)"
    )
  }
  
  # Column lookup helper
  get_col <- function(name) {
    idx <- which(colnames_lower == name)
    if (length(idx) == 0) return(NULL)
    colnames(gwas_data)[idx[1]]
  }
  
  # Map to internal names
  gwas_data$chr <- gwas_data[[get_col("chromosome")]]
  gwas_data$pos <- as.integer(gwas_data[[get_col("base_pair_location")]])
  gwas_data$pval <- as.numeric(gwas_data[[get_col("p_value")]])
  
  # Optional columns
  snp_col <- get_col("variant_id")
  if (!is.null(snp_col)) {
    gwas_data$snp <- gwas_data[[snp_col]]
  } else {
    gwas_data$snp <- paste0("chr", gwas_data$chr, ":", gwas_data$pos)
  }
  
  beta_col <- get_col("beta")
  if (!is.null(beta_col)) gwas_data$beta <- as.numeric(gwas_data[[beta_col]])
  
  se_col <- get_col("standard_error")
  if (!is.null(se_col)) gwas_data$se <- as.numeric(gwas_data[[se_col]])
  
  # Filter to region
  gwas_data$chr <- gsub("^chr", "", gwas_data$chr)
  region_chr <- gsub("^chr", "", as.character(seqnames(region)))
  
  gwas_data <- gwas_data[
    gwas_data$chr == region_chr &
      gwas_data$pos >= start(region) &
      gwas_data$pos <= end(region),
  ]
  gwas_data <- gwas_data[!is.na(gwas_data$pval), ]
  
  return(gwas_data)
}

#' Load LD data for a lead SNP
#'
#' @param ld.file Path to LD file
#' @return data.frame with columns: chr, pos, r2
#'
#' @details
#' Required columns: \code{chromosome}, \code{position}, \code{r2}
#'
#' Column matching is case-insensitive. Values should be pairwise r-squared to a lead SNP.
#'
#' @export
LoadLDData <- function(ld.file) {
  ld_data <- fread(ld.file, data.table = FALSE)
  colnames_lower <- tolower(colnames(ld_data))
  
  # Validate required columns
  required <- c("chromosome", "position", "r2")
  missing <- required[!required %in% colnames_lower]
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "), "\n",
      "Found: ", paste(colnames(ld_data), collapse = ", "), "\n",
      "Format: chromosome, position, r2"
    )
  }
  
  # Column lookup helper
  get_col <- function(name) {
    idx <- which(colnames_lower == name)
    colnames(ld_data)[idx[1]]
  }
  
  # Extract and standardize
  result <- data.frame(
    chr = gsub("^chr", "", as.character(ld_data[[get_col("chromosome")]])),
    pos = as.integer(ld_data[[get_col("position")]]),
    r2 = as.numeric(ld_data[[get_col("r2")]])
  )
  
  return(result)
}

#' Load fine-mapping credible sets
#'
#' @param credset.file Path to credible set file
#' @param credset.threshold Posterior inclusion probability threshold (default: 0.01)
#' @return data.frame with columns: chr, pos, pip, credset_id
#'
#' @details
#' Required columns: \code{chromosome}, \code{position}, \code{pip}, \code{cs}
#'
#' Column matching is case-insensitive. Variants with cs = -1 are not in a credible set.
#'
#' @export
LoadCredibleSets <- function(credset.file, credset.threshold = 0.01) {
  cs_data <- fread(credset.file, data.table = FALSE)
  colnames_lower <- tolower(colnames(cs_data))
  
  # Check for header
  if (all(grepl("^V[0-9]+$", colnames(cs_data)))) {
    stop("File has no header row")
  }
  
  # Validate required columns
  required <- c("chromosome", "position", "pip", "cs")
  missing <- required[!required %in% colnames_lower]
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "), "\n",
      "Found: ", paste(colnames(cs_data), collapse = ", "), "\n",
      "Format: chromosome, position, pip, cs"
    )
  }
  
  # Column lookup helper
  get_col <- function(name) {
    idx <- which(colnames_lower == name)
    colnames(cs_data)[idx[1]]
  }
  
  # Extract and standardize
  result <- data.frame(
    chr = gsub("^chr", "", as.character(cs_data[[get_col("chromosome")]])),
    pos = as.integer(cs_data[[get_col("position")]]),
    pip = as.numeric(cs_data[[get_col("pip")]]),
    credset_id = cs_data[[get_col("cs")]]
  )
  
  # Filter: PIP >= threshold AND in a credible set (cs != -1)
  result <- result[result$pip >= credset.threshold & result$credset_id != -1, ]
  
  if (nrow(result) == 0) {
    warning("No variants passed filters (pip >= ", credset.threshold, " and cs != -1)")
  }
  
  return(result)
}

#' Create GWAS Manhattan plot track
#'
#' @param region Genomic region (GRanges or string like "chr10-112900000-113100000")
#' @param gwas.file Path to GWAS summary statistics file
#' @param ld.file Path to LD file. Optional.
#' @param ld.lead.snp Lead SNP for LD (required if ld.file provided)
#' @param credset.file Path to credible set file. Optional.
#' @param credset.threshold PIP threshold for credible sets (default: 0.01)
#' @param p.threshold Genome-wide significance threshold (default: 5e-8)
#' @param ymax Maximum y-axis value (default: auto)
#' @param point.size Point size (default: 1)
#' @param point.color Point color when no LD (default: "steelblue")
#' @param show.axis Show x-axis (default: TRUE)
#' @return ggplot2 object
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_hline scale_y_continuous
#' @importFrom ggplot2 theme_classic labs theme element_blank element_line element_text
#' @importFrom ggplot2 scale_shape_manual scale_size_manual scale_color_manual
#' @export
#' @concept visualization
GWASTrack <- function(
    region,
    gwas.file,
    ld.file = NULL,
    ld.lead.snp = NULL,
    credset.file = NULL,
    credset.threshold = 0.01,
    p.threshold = 5e-8,
    ymax = NULL,
    point.size = 1,
    point.color = "steelblue",
    show.axis = TRUE
) {
  
  # Load GWAS data
  gwas_data <- LoadGWASRegion(gwas.file, region)
  if (nrow(gwas_data) == 0) {
    stop("No GWAS data found in region")
  }
  gwas_data$log10p <- -log10(gwas_data$pval)
  
  # Validate LD parameters
  if (!is.null(ld.file) && is.null(ld.lead.snp)) {
    stop("ld.lead.snp required when ld.file provided")
  }
  
  # LocusZoom colors
  ld_colors <- c(
    "r2_0-0.2" = "#0000CD",
    "r2_0.2-0.4" = "#00CED1",
    "r2_0.4-0.6" = "#32CD32",
    "r2_0.6-0.8" = "#FFA500",
    "r2_0.8-1.0" = "#FF0000"
  )
  
  # Merge LD data
  if (!is.null(ld.file)) {
    ld_data <- LoadLDData(ld.file)
    gwas_data <- merge(gwas_data, ld_data, by = c("chr", "pos"), all.x = TRUE)
    gwas_data$r2 <- as.numeric(gwas_data$r2)
    gwas_data$ld_category <- cut(
      gwas_data$r2,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c("r2_0-0.2", "r2_0.2-0.4", "r2_0.4-0.6", "r2_0.6-0.8", "r2_0.8-1.0"),
      include.lowest = TRUE
    )
  } else {
    gwas_data$r2 <- NA
    gwas_data$ld_category <- NA
  }
  
  # Merge credible set data
  if (!is.null(credset.file)) {
    credset_data <- LoadCredibleSets(credset.file, credset.threshold)
    gwas_data <- merge(gwas_data, credset_data, by = c("chr", "pos"), all.x = TRUE)
  } else {
    gwas_data$pip <- NA
    gwas_data$credset_id <- NA
  }
  
  # Y-axis limit
  if (is.null(ymax)) {
    ymax <- max(gwas_data$log10p, na.rm = TRUE) * 1.1
  }
  
  # Build plot
  if (!is.null(credset.file)) {
    gwas_data$in_credset <- !is.na(gwas_data$pip)
    
    if (!is.null(ld.file)) {
      # LD + credible sets
      p <- ggplot(gwas_data, aes(x = pos, y = log10p, color = ld_category)) +
        geom_point(aes(shape = in_credset, size = in_credset), alpha = 0.6) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
        scale_size_manual(values = c("FALSE" = point.size, "TRUE" = point.size * 2), guide = "none") +
        scale_color_manual(values = ld_colors, name = expression(LD~(r^2)), na.value = "grey50")
    } else {
      # Credible sets only
      p <- ggplot(gwas_data, aes(x = pos, y = log10p)) +
        geom_point(aes(shape = in_credset, size = in_credset, color = in_credset), alpha = 0.6) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
        scale_size_manual(values = c("FALSE" = point.size, "TRUE" = point.size * 2), guide = "none") +
        scale_color_manual(
          values = c("FALSE" = point.color, "TRUE" = "#FF0000"),
          name = "Credible set",
          labels = c("FALSE" = "No", "TRUE" = "Yes")
        )
    }
    
  } else if (!is.null(ld.file)) {
    # LD only
    p <- ggplot(gwas_data, aes(x = pos, y = log10p, color = ld_category)) +
      geom_point(size = point.size, alpha = 0.6) +
      scale_color_manual(values = ld_colors, name = expression(LD~(r^2)), na.value = "grey50")
    
  } else {
    # Basic plot
    p <- ggplot(gwas_data, aes(x = pos, y = log10p)) +
      geom_point(color = point.color, size = point.size, alpha = 0.6)
  }
  
  # Common elements
  p <- p +
    geom_hline(yintercept = -log10(p.threshold), linetype = "dashed", color = "red", alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    theme_classic() +
    labs(x = "Position (bp)", y = expression(-log[10](italic(P)))) +
    theme(
      axis.title.x = if (show.axis) element_text() else element_blank(),
      axis.text.x = if (show.axis) element_text() else element_blank(),
      axis.line.x = if (show.axis) element_line() else element_blank(),
      axis.ticks.x = if (show.axis) element_line() else element_blank()
    )
  
  return(p)
}
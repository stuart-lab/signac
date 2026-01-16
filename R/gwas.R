# R/gwas.R
# GWAS Visualization for Signac

#' @importFrom utils globalVariables
#'
NULL
globalVariables(names = c("in_credset", "ld_category", "log10p", "base_pair_location"), package = "Signac")

#' Load GWAS-SSF file
#'
#' @param gwas.file Path to GWAS/QTL summary statistics file (TSV)
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
#' 
#' @references \doi{10.1101/2022.07.15.500230}
#' @export
LoadGWAS <- function(gwas.file) {
  
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
  
  # # Map to internal names
  # gwas_data$chr <- gwas_data[[get_col("chromosome")]]
  # gwas_data$pos <- as.integer(gwas_data[[get_col("base_pair_location")]])
  # gwas_data$pval <- as.numeric(gwas_data[[get_col("p_value")]])
  
  # Optional columns
  snp_col <- get_col("variant_id")
  effect_allele_col <- get_col("effect_allele")
  other_allele_col <- get_col("other_allele")
  
  if (is.null(snp_col) & !(is.null(x = effect_allele_col) | is.null(x = other_allele_col))) {
    gwas_data$variant_id <- paste(
      gwas_data$chromosome, gwas_data$base_pair_location,
      gwas_data[[effect_allele_col]], gwas_data[[other_allele_col]], sep = "_"
    )
  }
  
  # beta_col <- get_col("beta")
  # if (!is.null(beta_col)) gwas_data$beta <- as.numeric(gwas_data[[beta_col]])
  # 
  # se_col <- get_col("standard_error")
  # if (!is.null(se_col)) gwas_data$se <- as.numeric(gwas_data[[se_col]])
  
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
#' @param gwas Path to GWAS summary statistics file, or a dataframe
#' containing the gwas data in the GWAS-SSF format
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
#' @importFrom Seqinfo seqnames
#' @importFrom GenomicRanges start end
#' @export
#' @concept visualization
GWASTrack <- function(
    gwas,
    region,
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
  if (is.character(x = gwas)) {
    gwas <- LoadGWAS(gwas)
  } 
  
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  
  # subset to region
  chromosome <- as.character(x = seqnames(x = region))
  gwas <- gwas[
    gwas$chromosome == chromosome &
      gwas$base_pair_location >= start(x = region) &
      gwas$base_pair_location <= end(x = region),
  ]
  gwas <- gwas[!is.na(gwas$p_value), ]
  
  
  if (nrow(x = gwas) == 0) {
    stop("No GWAS data found in region")
  }
  gwas$log10p <- -log10(gwas$p_value)
  
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
    # TODO need to update column names
    gwas <- merge(gwas, ld_data, by = c("chr", "pos"), all.x = TRUE)
    gwas$r2 <- as.numeric(gwas$r2)
    gwas$ld_category <- cut(
      gwas$r2,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c("r2_0-0.2", "r2_0.2-0.4", "r2_0.4-0.6", "r2_0.6-0.8", "r2_0.8-1.0"),
      include.lowest = TRUE
    )
  } else {
    gwas$r2 <- NA
    gwas$ld_category <- NA
  }
  
  # Merge credible set data
  if (!is.null(credset.file)) {
    credset_data <- LoadCredibleSets(credset.file, credset.threshold)
    # TODO need to update column names
    gwas <- merge(gwas, credset_data, by = c("chr", "pos"), all.x = TRUE)
    gwas$in_credset <- !is.na(gwas$pip)
  }
  
  # Y-axis limit
  if (is.null(ymax)) {
    ymax <- max(gwas$log10p, na.rm = TRUE) * 1.1
  }
  
  # Build plot
  if ("in_credset" %in% colnames(x = gwas)) {
    if (!is.null(ld.file)) {
      # LD + credible sets
      p <- ggplot(gwas, aes(x = base_pair_location, y = log10p, color = ld_category)) +
        geom_point(aes(shape = in_credset, size = in_credset), alpha = 0.6) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
        scale_size_manual(values = c("FALSE" = point.size, "TRUE" = point.size * 2), guide = "none") +
        scale_color_manual(values = ld_colors, name = expression(LD~(r^2)), na.value = "grey50")
    } else {
      # Credible sets only
      p <- ggplot(gwas, aes(x = base_pair_location, y = log10p)) +
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
    p <- ggplot(gwas, aes(x = base_pair_location, y = log10p, color = ld_category)) +
      geom_point(size = point.size, alpha = 0.6) +
      scale_color_manual(values = ld_colors, name = expression(LD~(r^2)), na.value = "grey50")
    
  } else {
    # Basic plot
    p <- ggplot(gwas, aes(x = base_pair_location, y = log10p)) +
      geom_point(color = point.color, size = point.size, alpha = 0.6)
  }
  
  # Common elements
  p <- p +
    geom_hline(yintercept = -log10(p.threshold), linetype = "dashed", color = "red", alpha = 0.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
    theme_classic() +
    labs(x = paste0(chromosome, " position (bp)"),
         y = expression(-log[10](italic(P)))) +
    theme(
      axis.title.x = if (show.axis) element_text() else element_blank(),
      axis.text.x = if (show.axis) element_text() else element_blank(),
      axis.line.x = if (show.axis) element_line() else element_blank(),
      axis.ticks.x = if (show.axis) element_line() else element_blank()
    ) +
    xlim(c(start(x = region), end(x = region)))
  
  return(p)
}
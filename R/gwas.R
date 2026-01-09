# R/gwas.R
# GWAS Visualization for Signac
#
# GWAS data must be in GWAS-SSF (Summary Statistics Format) as defined by
# NHGRI-EBI GWAS Catalog. See: https://github.com/EBISPOT/gwas-summary-statistics-standard
#
# Required columns: chromosome, base_pair_location, p_value
# Optional columns: effect_allele, other_allele, beta, standard_error,
#                   effect_allele_frequency, variant_id

#' @importFrom utils globalVariables
#'
NULL
globalVariables(names = c("in_credset", "ld_category", "log10p", "pos"), package = "Signac")

#' Load GWAS summary statistics for a genomic region
#'
#' Loads GWAS summary statistics in SSF (Summary Statistics Format) for a
#' specified genomic region. SSF is the standard format defined by NHGRI-EBI
#' GWAS Catalog.
#'
#' @param gwas.file Path to GWAS summary statistics file (tabix-indexed TSV in SSF format)
#' @param region Genomic region (GRanges object or string like "chr1:1000000-2000000")
#' @return data.frame with standardized columns: chr, pos, pval, snp, and optional
#'   beta, se, effect_allele, other_allele, eaf
#'
#' @details
#' GWAS data must be in SSF format with the following columns:
#' \itemize{
#'   \item Required: \code{chromosome}, \code{base_pair_location}, \code{p_value}
#'   \item Optional: \code{effect_allele}, \code{other_allele}, \code{beta},
#'         \code{standard_error}, \code{effect_allele_frequency}, \code{variant_id}
#' }
#'
#' Column matching is case-insensitive but column names must match exactly
#' (no fuzzy matching).
#'
#' To convert other formats to SSF:
#' \itemize{
#'   \item VCF (OpenGWAS): \code{bcftools query -f '\%CHROM\\t\%POS\\t\%ALT\\t\%REF[\\t\%ES\\t\%SE\\t\%AF\\t\%LP]\\n' gwas.vcf.gz | awk ...}
#'   \item FinnGen: \code{awk} to rename columns
#'   \item PLINK2: use \code{--gwas-ssf} flag
#' }
#' See vignette("gwas-formats") for detailed conversion instructions.
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges
#' @export
LoadGWASRegion <- function(gwas.file, region) {
  # Convert region to GRanges if needed
  if (!inherits(region, "GRanges")) {
    # Parse region string manually
    # Format: "chr10:112950000-113050000" or "chr10-112950000-113050000"
    region_str <- gsub(":", "-", region)
    parts <- strsplit(region_str, "-")[[1]]
    
    if (length(parts) != 3) {
      stop("Region must be in format 'chr:start-end' or 'chr-start-end'")
    }
    
    region <- GRanges(
      seqnames = parts[1],
      ranges = IRanges(
        start = as.numeric(parts[2]),
        end = as.numeric(parts[3])
      )
    )
  }
  
  # Read GWAS file
  gwas_data <- fread(gwas.file, data.table = FALSE)
  
  # Strict SSF column validation (case-insensitive matching)
  colnames_lower <- tolower(colnames(gwas_data))
  
  # Required SSF columns
  required_cols <- c("chromosome", "base_pair_location", "p_value")
  missing <- required_cols[!required_cols %in% colnames_lower]
  
  if (length(missing) > 0) {
    stop(
      "Missing required SSF columns: ", paste(missing, collapse = ", "), "\n",
      "Found columns: ", paste(colnames(gwas_data), collapse = ", "), "\n\n",
      "GWAS data must be in SSF format with columns: chromosome, base_pair_location, p_value\n",
      "See vignette('gwas-formats') for conversion instructions from VCF, FinnGen, PLINK, etc."
    )
  }
  
  # Helper function for case-insensitive column lookup
  get_col <- function(name) {
    idx <- which(colnames_lower == name)
    if (length(idx) == 0) return(NULL)
    colnames(gwas_data)[idx[1]]
  }
  
  # Get required columns
  
  chr_col <- get_col("chromosome")
  pos_col <- get_col("base_pair_location")
  p_col <- get_col("p_value")
  
  # Get optional SSF columns
  snp_col <- get_col("variant_id")
  beta_col <- get_col("beta")
  se_col <- get_col("standard_error")
  ea_col <- get_col("effect_allele")
  oa_col <- get_col("other_allele")
  eaf_col <- get_col("effect_allele_frequency")
  
  # Standardize to internal column names
  gwas_data$chr <- gwas_data[[chr_col]]
  gwas_data$pos <- as.integer(gwas_data[[pos_col]])
  gwas_data$pval <- as.numeric(gwas_data[[p_col]])
  
  # Handle variant_id (create position-based ID if not present)
  if (!is.null(snp_col)) {
    gwas_data$snp <- gwas_data[[snp_col]]
  } else {
    gwas_data$snp <- paste0("chr", gwas_data$chr, ":", gwas_data$pos)
  }
  
  # Map optional columns
  if (!is.null(beta_col)) gwas_data$beta <- as.numeric(gwas_data[[beta_col]])
  if (!is.null(se_col)) gwas_data$se <- as.numeric(gwas_data[[se_col]])
  if (!is.null(ea_col)) gwas_data$effect_allele <- gwas_data[[ea_col]]
  if (!is.null(oa_col)) gwas_data$other_allele <- gwas_data[[oa_col]]
  if (!is.null(eaf_col)) gwas_data$eaf <- as.numeric(gwas_data[[eaf_col]])
  
  # Handle chromosome naming (chr1 vs 1)
  gwas_data$chr <- gsub("^chr", "", gwas_data$chr)
  region_chr <- gsub("^chr", "", as.character(seqnames(region)))
  
  # Filter to region
  gwas_data <- gwas_data[
    gwas_data$chr == region_chr &
      gwas_data$pos >= start(region) &
      gwas_data$pos <= end(region),
  ]
  
  # Remove rows with missing p-values
  gwas_data <- gwas_data[!is.na(gwas_data$pval), ]
  
  return(gwas_data)
}

#' Load LD data from LDlink or PLINK
#'
#' @param ld.file Path to LD file (LDlink or PLINK pairwise format)
#' @return data.frame with columns: chr, pos, r2
#' @details LD must be pairwise r² relative to one lead SNP.
#'   Matched to GWAS by position (CHR + POS).
#'
#'   Supported formats:
#'   \itemize{
#'     \item LDlink: requires \code{Coord} and \code{R2} columns
#'     \item PLINK pairwise: requires \code{CHR_B}, \code{BP_B}, and \code{R2} columns
#'   }
#' @export
LoadLDData <- function(ld.file) {
  ld_data <- fread(ld.file, data.table = FALSE)
  
  # LDlink format
  if ("Coord" %in% colnames(ld_data) && "R2" %in% colnames(ld_data)) {
    coord_split <- strsplit(ld_data$Coord, ":", fixed = TRUE)
    ld_data$chr <- sapply(coord_split, function(x) gsub("chr", "", x[1]))
    ld_data$pos <- as.integer(sapply(coord_split, function(x) x[2]))
    ld_data$r2 <- ld_data$R2
    return(ld_data[, c("chr", "pos", "r2")])
  }
  
  # PLINK pairwise format
  if ("CHR_B" %in% colnames(ld_data) && "BP_B" %in% colnames(ld_data) && "R2" %in% colnames(ld_data)) {
    ld_data$chr <- gsub("chr", "", as.character(ld_data$CHR_B))
    ld_data$pos <- as.integer(ld_data$BP_B)
    ld_data$r2 <- as.numeric(ld_data$R2)
    return(ld_data[, c("chr", "pos", "r2")])
  }
  
  stop(
    "Unrecognized LD format.\n",
    "Supported formats:\n",
    "  - LDlink: requires 'Coord' and 'R2' columns\n",
    "  - PLINK pairwise: requires 'CHR_B', 'BP_B', and 'R2' columns"
  )
}

#' Create GWAS Manhattan plot track
#'
#' Creates a Manhattan plot for GWAS summary statistics, optionally colored by
#' LD (r²) to a lead SNP and/or highlighting fine-mapped credible set variants.
#'
#' @param region Genomic region to plot (GRanges or string like "chr10:112900000-113100000")
#' @param gwas.file Path to GWAS summary statistics file (SSF format, tabix-indexed)
#' @param ld.file Path to LD data file (LDlink or PLINK pairwise format). Optional.
#' @param ld.lead.snp rsID that LD was calculated relative to (required if ld.file provided)
#' @param credset.file Path to credible set file (SuSiE or FINEMAP format). Optional.
#' @param credset.threshold Posterior inclusion probability threshold for credible sets (default: 0.01)
#' @param p.threshold Genome-wide significance threshold (default: 5e-8)
#' @param ymax Maximum y-axis value (default: auto-scale)
#' @param point.size Size of points (default: 1)
#' @param point.color Color for points when no LD coloring (default: "steelblue")
#' @param show.axis Show x-axis labels and ticks (default: TRUE)
#' @return ggplot2 object
#'
#' @details
#' LD coloring uses standard LocusZoom colors:
#' \itemize{
#'   \item r² 0.8-1.0: Red
#'   \item r² 0.6-0.8: Orange
#'   \item r² 0.4-0.6: Green
#'   \item r² 0.2-0.4: Cyan
#'   \item r² 0-0.2: Dark blue
#' }
#'
#' Credible set variants are shown as larger diamonds.
#'
#' @importFrom ggplot2 ggplot geom_point aes geom_hline scale_y_continuous
#' @importFrom ggplot2 theme_classic labs theme element_blank element_line element_text
#' @importFrom ggplot2 scale_shape_manual scale_size_manual scale_color_manual ggsave
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
  # Load data
  gwas_data <- LoadGWASRegion(gwas.file, region)
  
  if (nrow(gwas_data) == 0) {
    stop("No GWAS data found in region")
  }
  
  # Calculate -log10(p)
  gwas_data$log10p <- -log10(gwas_data$pval)
  
  # Validate LD parameters
  if (!is.null(ld.file) && is.null(ld.lead.snp)) {
    stop("ld.lead.snp required when ld.file provided")
  }
  
  # Standard LocusZoom colors (defined once, used in multiple places)
  ld_colors <- c(
    "r2_0-0.2" = "#0000CD",
    "r2_0.2-0.4" = "#00CED1",
    "r2_0.4-0.6" = "#32CD32",
    "r2_0.6-0.8" = "#FFA500",
    "r2_0.8-1.0" = "#FF0000"
  )
  
  # Merge LD data if provided
  if (!is.null(ld.file)) {
    ld_data <- LoadLDData(ld.file)
    
    # Position-based matching
    if ("chr" %in% colnames(ld_data) && "pos" %in% colnames(ld_data)) {
      gwas_data <- merge(gwas_data, ld_data, by = c("chr", "pos"), all.x = TRUE)
    } else if ("snp" %in% colnames(ld_data)) {
      gwas_data <- merge(gwas_data, ld_data, by = "snp", all.x = TRUE)
    }
    
    # Ensure r2 is numeric
    gwas_data$r2 <- as.numeric(gwas_data$r2)
    
    # Bin r² values into standard LD categories
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
  
  # Merge credible set data if provided
  if (!is.null(credset.file)) {
    credset_data <- LoadCredibleSets(credset.file, credset.threshold)
    
    # Position-based matching
    if ("chr" %in% colnames(credset_data)) {
      gwas_data <- merge(gwas_data, credset_data, by = c("chr", "pos"), all.x = TRUE)
    } else {
      gwas_data <- merge(gwas_data, credset_data, by = "pos", all.x = TRUE)
    }
  } else {
    gwas_data$pip <- NA
    gwas_data$credset_id <- NA
  }
  
  # Set y-axis limit
  if (is.null(ymax)) {
    ymax <- max(gwas_data$log10p, na.rm = TRUE) * 1.1
  }
  
  # Create plot based on what data is available
  if (!is.null(credset.file)) {
    # Highlight credible set variants
    gwas_data$in_credset <- !is.na(gwas_data$pip)
    
    if (!is.null(ld.file)) {
      # Both LD and credible sets
      p <- ggplot(gwas_data, aes(x = pos, y = log10p, color = ld_category)) +
        geom_point(aes(shape = in_credset, size = in_credset), alpha = 0.6) +
        scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 18),
          guide = "none"
        ) +
        scale_size_manual(
          values = c("FALSE" = point.size, "TRUE" = point.size * 2),
          guide = "none"
        ) +
        scale_color_manual(
          values = ld_colors,
          name = expression(LD~(r^2)),
          na.value = "grey50"
        )
    } else {
      # Credible sets without LD
      p <- ggplot(gwas_data, aes(x = pos, y = log10p)) +
        geom_point(aes(shape = in_credset, size = in_credset,
                       color = in_credset), alpha = 0.6) +
        scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 18),
          guide = "none"
        ) +
        scale_size_manual(
          values = c("FALSE" = point.size, "TRUE" = point.size * 2),
          guide = "none"
        ) +
        scale_color_manual(
          values = c("FALSE" = point.color, "TRUE" = "#FF0000"),
          name = "In credible set",
          labels = c("FALSE" = "No", "TRUE" = "Yes")
        )
    }
    
  } else if (!is.null(ld.file)) {
    # LD coloring only (no credible sets)
    p <- ggplot(gwas_data, aes(x = pos, y = log10p, color = ld_category)) +
      geom_point(size = point.size, alpha = 0.6) +
      scale_color_manual(
        values = ld_colors,
        name = expression(LD~(r^2)),
        na.value = "grey50"
      )
    
  } else {
    # Basic plot (no LD, no credible sets)
    p <- ggplot(gwas_data, aes(x = pos, y = log10p)) +
      geom_point(color = point.color, size = point.size, alpha = 0.6)
  }
  
  # Add common elements
  p <- p +
    geom_hline(
      yintercept = -log10(p.threshold),
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
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

#' Load fine-mapping credible sets from SuSiE or FINEMAP
#'
#' @param credset.file Path to credible set file
#' @param credset.threshold Posterior inclusion probability threshold (default: 0.01)
#' @return data.frame with columns: chr, pos, pip, credset_id
#'
#' @details
#' Supports SuSiE and FINEMAP output formats. Column detection is flexible to
#' handle variations in output from different tools.
#'
#' Required columns (case-insensitive, multiple names accepted):
#' \itemize{
#'   \item Position: \code{pos}, \code{position}, or \code{bp}
#'   \item Chromosome: \code{chr}, \code{chromosome}, \code{chrom}, or \code{seqnames}
#'   \item Posterior probability: \code{prob}, \code{pip}, or \code{pp}
#'   \item Credible set ID: \code{cs}, \code{credset}, or \code{credible_set}
#' }
#'
#' Variants are filtered to those with PIP >= threshold AND assigned to a

#' credible set (cs != -1).
#'
#' @export
LoadCredibleSets <- function(credset.file, credset.threshold = 0.01) {
  cs_data <- fread(credset.file, data.table = FALSE)
  
  # Check for required columns
  colnames_lower <- tolower(colnames(cs_data))
  
  if (all(grepl("^V[0-9]+$", colnames(cs_data)))) {
    stop("Credible set file has no header. Please add column names or use the full file with headers.")
  }
  
  # Find position
  pos_col <- colnames(cs_data)[which(colnames_lower %in% c("pos", "position", "bp"))[1]]
  if (is.na(pos_col)) stop("Could not find position column (pos/position/bp)")
  
  # Find chromosome
  chr_col <- colnames(cs_data)[which(colnames_lower %in% c("chr", "chromosome", "chrom", "seqnames"))[1]]
  if (is.na(chr_col)) stop("Could not find chromosome column (chr/chromosome)")
  
  # Find posterior probability
  pp_col <- colnames(cs_data)[which(colnames_lower %in% c("prob", "pip", "pp"))[1]]
  if (is.na(pp_col)) stop("Could not find posterior probability column (prob/pip/pp)")
  
  # Find credible set ID
  cs_col <- colnames(cs_data)[which(colnames_lower %in% c("cs", "credset", "credible_set"))[1]]
  if (is.na(cs_col)) stop("Could not find credible set column (cs/credset)")
  
  # Extract and standardize
  result <- data.frame(
    chr = gsub("chr", "", as.character(cs_data[[chr_col]])),
    pos = as.integer(cs_data[[pos_col]]),
    pip = as.numeric(cs_data[[pp_col]]),
    credset_id = cs_data[[cs_col]]
  )
  
  # Filter: PIP >= threshold AND in a credible set (cs != -1)
  result <- result[result$pip >= credset.threshold & result$credset_id != -1, ]
  
  if (nrow(result) == 0) {
    warning("No variants in credible sets above threshold")
  }
  
  return(result)
}
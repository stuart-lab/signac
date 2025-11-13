# R/gwas.R
# GWAS Visualization for Signac
# Version 1: Basic Manhattan plot only

#' Load GWAS summary statistics for a genomic region
#' 
#' @param gwas.file Path to GWAS summary statistics file
#' @param region Genomic region (GRanges object or string like "chr1:1000000-2000000")
#' @return data.frame with standardized columns
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
  
  # Auto-detect column names (handle variations)
  colnames_lower <- tolower(colnames(gwas_data))
  
  # Find CHR column
  chr_col <- colnames(gwas_data)[which(colnames_lower %in% c("chr", "#chr", "chrom", "#chrom", "chromosome", "#chromosome"))[1]]
  if (is.na(chr_col)) stop("Could not find chromosome column")
  
  # Find POS column  
  pos_col <- colnames(gwas_data)[which(colnames_lower %in% c("pos", "bp", "position"))[1]]
  if (is.na(pos_col)) stop("Could not find position column")
  
  # Find P column
  p_col <- colnames(gwas_data)[which(colnames_lower %in% c("p", "pval", "pvalue", "p_value"))[1]]
  if (is.na(p_col)) stop("Could not find p-value column")
  
  # Optional columns
  snp_col <- colnames(gwas_data)[which(colnames_lower %in% c("snp", "rsid", "rs", "variant_id", "id"))[1]]
  beta_col <- colnames(gwas_data)[which(colnames_lower %in% c("beta", "b", "effect"))[1]]
  se_col <- colnames(gwas_data)[which(colnames_lower %in% c("se", "stderr", "standard_error"))[1]]
  
  # Standardize column names
  gwas_data$chr <- gwas_data[[chr_col]]
  gwas_data$pos <- gwas_data[[pos_col]]
  gwas_data$pval <- gwas_data[[p_col]]
  
  if (!is.na(snp_col)) gwas_data$snp <- gwas_data[[snp_col]]
  if (!is.na(beta_col)) gwas_data$beta <- gwas_data[[beta_col]]
  if (!is.na(se_col)) gwas_data$se <- gwas_data[[se_col]]
  
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

#' Create GWAS Manhattan plot track
#'
#' @param region Genomic region to plot
#' @param gwas.file Path to GWAS summary statistics
#' @param p.threshold Significance threshold (default: 5e-8)
#' @param ymax Maximum y-axis value (default: auto-scale)
#' @param point.size Size of points (default: 1)
#' @param point.color Color for points (default: "steelblue")
#' @return ggplot2 object
#' @importFrom ggplot2 ggplot geom_point aes geom_hline scale_y_continuous
#' @importFrom ggplot2 theme_classic labs theme element_blank ggsave
#' @export
#' @concept visualization
GWASTrack <- function(
    region,
    gwas.file,
    p.threshold = 5e-8,
    ymax = NULL,
    point.size = 1,
    point.color = "steelblue"
) {
  # Load data
  gwas_data <- LoadGWASRegion(gwas.file, region)
  
  if (nrow(gwas_data) == 0) {
    stop("No GWAS data found in region")
  }
  
  # Calculate -log10(p)
  gwas_data$log10p <- -log10(gwas_data$pval)
  
  # Set y-axis limit
  if (is.null(ymax)) {
    ymax <- max(gwas_data$log10p, na.rm = TRUE) * 1.1
  }
  
  # Create plot
  p <- ggplot(gwas_data, aes(x = pos, y = log10p)) +
    geom_point(color = point.color, size = point.size, alpha = 0.6) +
    geom_hline(
      yintercept = -log10(p.threshold),
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, ymax)
    ) +
    theme_classic() +
    labs(x = "Position (bp)", y = expression(-log[10](italic(P)))) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  return(p)
}
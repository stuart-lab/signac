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
  snp_col <- colnames(gwas_data)[which(colnames_lower %in% c("snp", "rsid", "rsids", "rs", "variant_id", "id"))[1]]
  beta_col <- colnames(gwas_data)[which(colnames_lower %in% c("beta", "b", "effect"))[1]]
  se_col <- colnames(gwas_data)[which(colnames_lower %in% c("se", "stderr", "standard_error"))[1]]
  
  # Standardize column names
  gwas_data$chr <- gwas_data[[chr_col]]
  gwas_data$pos <- gwas_data[[pos_col]]
  gwas_data$pval <- gwas_data[[p_col]]
  
  if (!is.na(snp_col)) {
    gwas_data$snp <- gwas_data[[snp_col]]
  } else {
    # If no SNP ID, create position-based ID
    gwas_data$snp <- paste0("chr", gwas_data$chr, ":", gwas_data$pos)
  }
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

#' Load LD data from LDlink or PLINK
#' 
#' @param ld.file Path to LD file (LDlink or PLINK pairwise format)
#' @details LD must be pairwise r² relative to one lead SNP.
#'   Matched to GWAS by position (CHR + POS).
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
  
  stop("Unrecognized LD format. Needs LDlink (Coord, R2) or PLINK (CHR_B, BP_B, R2)")
}

#' Create GWAS Manhattan plot track
#'
#' @param region Genomic region to plot
#' @param gwas.file Path to GWAS summary statistics
#' @param p.threshold Significance threshold (default: 5e-8)
#' @param ymax Maximum y-axis value (default: auto-scale)
#' @param point.size Size of points (default: 1)
#' @param point.color Color for points (default: "steelblue")
#' @param ld.lead.snp rsID that LD was calculated relative to (required if ld.file provided)
#' @return ggplot2 object
#' @importFrom ggplot2 ggplot geom_point aes geom_hline scale_y_continuous
#' @importFrom ggplot2 theme_classic labs theme element_blank ggsave
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
    
  } else {
    gwas_data$r2 <- NA
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
  
  # Create plot
  if (!is.null(credset.file)) {
    # Highlight credible set variants
    gwas_data$in_credset <- !is.na(gwas_data$pip)
    
    if (!is.null(ld.file)) {
      # Both LD and credible sets
      p <- ggplot(gwas_data, aes(x = pos, y = log10p, color = ld_category)) +
        geom_point(aes(shape = in_credset, size = in_credset), alpha = 0.6) +
        scale_shape_manual(
          values = c("FALSE" = 16, "TRUE" = 18),  # circle vs diamond
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
        axis.title.x = if(show.axis) element_text() else element_blank(),
        axis.text.x = if(show.axis) element_text() else element_blank(),
        axis.line.x = if(show.axis) element_line() else element_blank(),
        axis.ticks.x = if(show.axis) element_line() else element_blank()
      )
  
  } else if (!is.null(ld.file)) {
    # Bin r² values into standard LD categories
    gwas_data$ld_category <- cut(
      gwas_data$r2,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c("r2_0-0.2", "r2_0.2-0.4", "r2_0.4-0.6", "r2_0.6-0.8", "r2_0.8-1.0"),
      include.lowest = TRUE
    )
    
    # Standard LocusZoom colors (blue to red gradient)
    ld_colors <- c(
      "r2_0-0.2" = "#0000CD",      # Dark blue
      "r2_0.2-0.4" = "#00CED1",    # Cyan
      "r2_0.4-0.6" = "#32CD32",    # Green  
      "r2_0.6-0.8" = "#FFA500",    # Orange
      "r2_0.8-1.0" = "#FF0000"     # Red
    )
    
    # Create plot with LD coloring
    p <- ggplot(gwas_data, aes(x = pos, y = log10p, color = ld_category)) +
      geom_point(size = point.size, alpha = 0.6) +
      scale_color_manual(
        values = ld_colors,
        name = expression(LD~(r^2)),
        na.value = "grey50"  # SNPs without LD info
      ) +
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
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      theme(
        axis.title.x = if(show.axis) element_text() else element_blank(),
        axis.text.x = if(show.axis) element_text() else element_blank(),
        axis.line.x = if(show.axis) element_line() else element_blank(),
        axis.ticks.x = if(show.axis) element_line() else element_blank()
      )
    
  } else {
    # Original plot without LD coloring
    p <- ggplot(gwas_data, aes(x = pos, y = log10p)) +
      geom_point(color = point.color, size = point.size, alpha = 0.6) +
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
        axis.title.x = if(show.axis) element_text() else element_blank(),
        axis.text.x = if(show.axis) element_text() else element_blank(),
        axis.line.x = if(show.axis) element_line() else element_blank(),
        axis.ticks.x = if(show.axis) element_line() else element_blank()
      )
  }
  
  return(p)
}

#' Load fine-mapping credible sets from SuSiE or FINEMAP
#' 
#' @param credset.file Path to credible set file
#' @param credset.threshold Posterior probability threshold (default: 0.01)
#' @details Supports SuSiE and FINEMAP output formats. Matched by position.
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
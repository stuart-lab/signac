# R/gwas.R
# GWAS Visualization for Signac

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
  gwas_data <- fread(file = gwas.file, data.table = FALSE)
  colnames_lower <- tolower(x = colnames(x = gwas_data))
  
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

  # Optional columns
  variant_present <- "variant_id" %in% colnames(x = gwas_data)
  effect_allele_present <- "effect_allele" %in% colnames(x = gwas_data)
  other_allele_present <- "other_allele" %in% colnames(x = gwas_data)
  
  if (!variant_present & (effect_allele_present & other_allele_present)) {
    # fill in the variant id
    gwas_data[['variant_id']] <- paste(
      gwas_data[['chromosome']], gwas_data[['base_pair_location']],
      gwas_data[['effect_allele']], gwas_data[['other_allele']], sep = "_"
    )
  }
  
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
#' @importFrom data.table fread
#' @export
LoadLDData <- function(ld.file) {
  ld_data <- fread(input = ld.file, data.table = FALSE)
  colnames_lower <- tolower(x = colnames(x = ld_data))
  
  # Validate required columns
  required <- c("chromosome", "position", "r2")
  missing <- required[!required %in% colnames_lower]
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "), "\n",
      "Found: ", paste(colnames(x = ld_data), collapse = ", "), "\n",
      "Format: chromosome, position, r2"
    )
  }
  
  # Extract and standardize
  result <- data.frame(
    chromosome = as.character(x = ld_data[["chromosome"]]),
    base_pair_location = as.integer(x = ld_data[["position"]]),
    r2 = as.numeric(x = ld_data[["r2"]])
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
#' @importFrom data.table fread
#'
#' @export
LoadCredibleSets <- function(credset.file, credset.threshold = 0.01) {
  cs_data <- fread(input = credset.file, data.table = FALSE)
  colnames_lower <- tolower(x = colnames(x = cs_data))
  
  # Check for header
  if (all(grepl("^V[0-9]+$", colnames(x = cs_data)))) {
    stop("File has no header row")
  }
  
  # Validate required columns
  required <- c("chromosome", "position", "pip", "cs")
  missing <- required[!required %in% colnames_lower]
  if (length(missing) > 0) {
    stop(
      "Missing required columns: ", paste(missing, collapse = ", "), "\n",
      "Found: ", paste(colnames(x = cs_data), collapse = ", "), "\n",
      "Format: chromosome, position, pip, cs"
    )
  }
  
  # Extract and standardize
  result <- data.frame(
    chromosome = as.character(x = cs_data[["chromosome"]]),
    base_pair_location = as.integer(x = cs_data[["position"]]),
    pip = as.numeric(x = cs_data[["pip"]]),
    credset_id = cs_data[["cs"]]
  )
  
  # Filter: PIP >= threshold AND in a credible set (cs != -1)
  result <- result[result[['pip']] >= 
                     credset.threshold & result[['credset_id']] != -1, ]
  
  if (nrow(x = result) == 0) {
    warning("No variants passed filters (pip >= ", credset.threshold, " and cs != -1)")
  }
  
  return(result)
}

#' Create GWAS locus zoom track
#' 
#' This function will plot the p-values associated with variants in a given
#' region of the genome (genome position: x-axis; -log10(p): y-axis).
#' 
#' If an LD file is provided using the \code{ld.file} parameter, or a column
#' named \code{r2} is present in the input GWAS data, the variants will be
#' colored according to their LD value.
#' 
#' If a fine mapping file is provided using the \code{credset.file} parameter,
#' or a column named \code{in_credset} is present in the input GWAS data,
#' variants in the credible set will be denoted by shape. If LD information is
#' not provided, credible variants will also be denoted by color.
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
#' @importFrom ggplot2 ggplot geom_point aes_string geom_hline scale_y_continuous
#' theme_classic labs theme element_blank element_line element_text
#' scale_shape_manual scale_size_manual scale_color_manual
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
    gwas <- LoadGWAS(gwas.file = gwas)
  } 
  
  if (!inherits(x = region, what = "GRanges")) {
    region <- StringToGRanges(regions = region)
  }
  
  # subset to region
  chromosome <- as.character(x = seqnames(x = region))
  gwas <- gwas[
    gwas[['chromosome']] == chromosome &
      gwas[['base_pair_location']] >= start(x = region) &
      gwas[['base_pair_location']] <= end(x = region),
  ]
  gwas <- gwas[!is.na(x = gwas[['p_value']]), ]
  
  if (nrow(x = gwas) == 0) {
    stop("No GWAS data found in region")
  }
  gwas[['log10p']] <- -log10(x = gwas[['p_value']])
  
  # Validate LD parameters
  if (!is.null(x = ld.file) && is.null(x = ld.lead.snp)) {
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
  if (!is.null(x = ld.file)) {
    ld_data <- LoadLDData(ld.file = ld.file)
    gwas <- merge(
      x = gwas,
      y = ld_data,
      by = c("chromosome", "base_pair_location"),
      all.x = TRUE
    )
    gwas[['r2']] <- as.numeric(x = gwas[['r2']])
  }
  
  if ('r2' %in% colnames(x = gwas)) {
    gwas[['ld_category']] <- cut(
      gwas[['r2']],
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
      labels = c("r2_0-0.2", "r2_0.2-0.4", "r2_0.4-0.6", "r2_0.6-0.8", "r2_0.8-1.0"),
      include.lowest = TRUE
    )
  }
  
  # Merge credible set data
  if (!is.null(x = credset.file)) {
    credset_data <- LoadCredibleSets(
      credset.file = credset.file, credset.threshold = credset.threshold
    )
    gwas <- merge(
      x = gwas,
      y = credset_data,
      by = c("chromosome", "base_pair_location"),
      all.x = TRUE
    )
    gwas[['in_credset']] <- !is.na(x = gwas[['pip']])
  }
  
  # Y-axis limit
  if (is.null(x = ymax)) {
    ymax <- max(gwas[['log10p']], na.rm = TRUE) * 1.1
  }
  
  # Build plot
  if ("in_credset" %in% colnames(x = gwas)) {
    if ("ld_category" %in% colnames(x = gwas)) {
      # LD + credible sets
      p <- ggplot(data = gwas, mapping = aes_string(
        x = 'base_pair_location',
        y = 'log10p',
        color = 'ld_category'
      )) +
        geom_point(
          aes_string(shape = 'in_credset', size = 'in_credset'), alpha = 0.6
        ) +
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
      # Credible sets only
      p <- ggplot(data = gwas, mapping = aes_string(
        x = 'base_pair_location',
        y = 'log10p')) +
        geom_point(
          mapping = aes_string(
            shape = 'in_credset',
            size = 'in_credset',
            color = 'in_credset'),
          alpha = 0.6
        ) +
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
          name = "Credible set",
          labels = c("FALSE" = "No", "TRUE" = "Yes")
        )
    }
    
  } else if ("ld_category" %in% colnames(x = gwas)) {
    # LD only
    p <- ggplot(
      data = gwas, mapping = aes_string(
        x = 'base_pair_location',
        y = 'log10p',
        color = 'ld_category'
        )
      ) +
      geom_point(size = point.size, alpha = 0.6) +
      scale_color_manual(
        values = ld_colors,
        name = expression(LD~(r^2)),
        na.value = "grey50"
      )
  } else {
    # Basic plot
    p <- ggplot(
      data = gwas, mapping = aes_string(
        x = 'base_pair_location',
        y = 'log10p'
        )
      ) +
      geom_point(color = point.color, size = point.size, alpha = 0.6)
  }
  
  # Common elements
  p <- p +
    geom_hline(
      yintercept = -log10(x = p.threshold),
      linetype = "dashed",
      color = "red",
      alpha = 0.5
    ) +
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
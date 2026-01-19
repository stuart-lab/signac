#' Load GWAS-SSF file
#'
#' @param gwas.file Path to GWAS/QTL summary statistics file (TSV)
#' @return data.frame with following columns present: chromosome, base_pair_location, p_value
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
#' @return data.frame with columns: chromosome, base_pair_location, r2
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
#' @return data.frame with columns: chromosome, base_pair_location, pip, credset_id
#'
#' @details
#' Required columns: \code{chromosome}, \code{position}, \code{pip}, \code{credset_id}
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
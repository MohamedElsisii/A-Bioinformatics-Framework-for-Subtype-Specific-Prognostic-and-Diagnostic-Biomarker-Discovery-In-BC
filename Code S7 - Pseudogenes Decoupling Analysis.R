# ==================================
# Script: Code S7 - Pseudogenes Decoupling Analysis.R
#
# Purpose:
# This script investigates the regulatory relationship between diagnostic
# pseudogenes and their protein-coding parent genes. It calculates the
# correlation between each pseudogene-parent pair in normal tissue vs.
# tumor subtypes and tests for a significant change in this correlation
# (regulatory decoupling or reinforcement).
#
# Key Outputs:
# 1. A CSV file containing correlation stats (normal & tumor), the difference
#    (Delta_r), and statistical significance for each pair.
# 2. A bar plot visualizing the Delta_r for all tested pairs, with significant
#    changes highlighted.
# ==================================

# Load necessary libraries, suppressing startup messages
suppressPackageStartupMessages({
  library(dplyr)     # For data manipulation (filter, bind_rows, etc.)
  library(ggplot2)   # For creating the visualization plot
})

# -------------------------------------------------------
# 1. Load Data
# -------------------------------------------------------
# Load the full tumor expression matrix
brca_expr <- readRDS("BRCA.rda")      
# Load the normal tissue expression matrix (control group)
normal_expr <- readRDS("Normal.rds") 
# Load sample IDs for the subtypes we are analyzing
her2_samples <- readRDS("Her2_Samples.rda")
basal_samples <- readRDS("Basal_Samples.rda")

# Clean TCGA sample names to standard format (replace dots with dashes)
her2_samples <- gsub("\\.", "-", her2_samples)
basal_samples <- gsub("\\.", "-", basal_samples)

# -------------------------------------------------------
# 2. Define Pseudogene-Parent Pairs
# -------------------------------------------------------
# This table defines which pseudogenes to test and their corresponding parent genes.
# In a production environment, this might be loaded from an external file.
pseudo_parents <- data.frame(
  Pseudogene = c("HMGN2P15", "TPTE2P2", "TMEM161BP1", "VN1R53P", "TNRC18P1", 
                 "AIDAP2", "CPHL1P", "RAB6C-AS1", "GTF2IP7", "RAD17P1", 
                 "POTEKP", "AZGP1P1", "PLAC9P1", "CYP4Z2P", "RN7SL314P"),
  ParentGene = c("HMGN2", "TPTE2", "TMEM161B", "VN1R53", "TNRC18", 
                 "AIDA", "CPHL1", "RAB6C", "GTF2I", "RAD17", 
                 "POTEK", "AZGP1", "PLAC9", "CYP4Z1", "7SL")
)

# -------------------------------------------------------
# 3. Core Functions
# -------------------------------------------------------

#' @title Fisher's Z-test for Correlation Difference
#' @description Tests if two correlation coefficients (r1, r2) from independent
#' samples (sizes n1, n2) are statistically different.
#' @param r1 Correlation coefficient 1 (e.g., normal tissue)
#' @param n1 Sample size 1
#' @param r2 Correlation coefficient 2 (e.g., tumor tissue)
#' @param n2 Sample size 2
#' @return Two-tailed p-value
fisher_z_test <- function(r1, n1, r2, n2) {
  # Handle cases where correlation could not be computed (NA)
  if (is.na(r1) || is.na(r2)) return(NA)
  
  # Fisher r-to-z transformation
  z1 <- 0.5 * log((1 + r1) / (1 - r1))
  z2 <- 0.5 * log((1 + r2) / (1 - r2))
  
  # Standard error of the difference
  se_diff <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  
  # Z-score for the difference
  z_score <- (z1 - z2) / se_diff
  
  # Calculate two-tailed p-value from normal distribution
  p_val <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  return(p_val)
}

#' @title Calculate Decoupling Statistics for One Pair
#' @description Computes Spearman correlations in normal and subtype tumor samples,
#' calculates the difference (Delta_r), and tests for significance.
#' @param pseudogene Symbol of the pseudogene
#' @param parent Symbol of the parent gene
#' @param normal_expr Expression matrix for normal samples
#' @param subtype_expr Expression matrix for tumor subtype samples
#' @param subtype_name Name of the subtype (for labeling results)
#' @return A data frame with one row containing all statistics, or NULL if genes are missing.
calc_decoupling <- function(pseudogene, parent, normal_expr, subtype_expr, subtype_name) {
  
  # 1. Check if both genes exist in the expression matrices
  # We check subtype_expr, assuming normal_expr has the same gene rows.
  if (!all(c(pseudogene, parent) %in% rownames(subtype_expr))) {
    warning(paste("Skipping pair:", pseudogene, "-", parent, "for subtype", subtype_name, "(one or both not found)."))
    return(NULL) 
  }
  
  # 2. Extract expression vectors
  ps_norm <- as.numeric(normal_expr[pseudogene, ])
  pa_norm <- as.numeric(normal_expr[parent, ])
  ps_sub  <- as.numeric(subtype_expr[pseudogene, ])
  pa_sub  <- as.numeric(subtype_expr[parent, ])
  
  # 3. Calculate Spearman correlations
  # suppressWarnings is used to handle potential ties or near-constant values gracefully
  r_norm <- suppressWarnings(cor(ps_norm, pa_norm, method = "spearman", use = "pairwise.complete.obs"))
  r_sub  <- suppressWarnings(cor(ps_sub, pa_sub, method = "spearman", use = "pairwise.complete.obs"))
  
  # 4. Get effective sample sizes (number of non-NA pairs)
  n1 <- sum(complete.cases(ps_norm, pa_norm))
  n2 <- sum(complete.cases(ps_sub, pa_sub))
  
  # 5. Perform statistical test
  p_val <- fisher_z_test(r_norm, n1, r_sub, n2)
  
  # 6. Calculate Delta_r (Change in correlation: Tumor - Normal)
  # Positive Delta_r = Reinforced correlation in tumor
  # Negative Delta_r = Decoupled (lost) correlation in tumor
  delta_r <- r_sub - r_norm 
  
  # 7. Return results as a data frame
  data.frame(
    Subtype = subtype_name,
    Pseudogene = pseudogene,
    Parent = parent,
    R_normal = r_norm,
    R_subtype = r_sub,
    Delta_r = delta_r,
    P_value = p_val
  )
}

# -------------------------------------------------------
# 4. Analysis Execution
# -------------------------------------------------------

# --- Analyze HER2-enriched Subtype ---
# Subset HER2 samples
her2_expr <- brca_expr[, her2_samples]
# Analyze the specific HMGN2P15 pair for this subtype
her2_results <- calc_decoupling("HMGN2P15", "HMGN2", normal_expr, her2_expr, "HER2-enriched")

# --- Analyze Basal-like Subtype ---
# Subset Basal-like samples
basal_expr <- brca_expr[, basal_samples]
# Filter the master pair list to exclude the HER2-specific one
basal_pseudo_parents <- pseudo_parents %>%
  filter(Pseudogene != "HMGN2P15")

# Loop through all remaining pairs for the Basal-like subtype
basal_results_list <- lapply(1:nrow(basal_pseudo_parents), function(i) {
  calc_decoupling(basal_pseudo_parents$Pseudogene[i], basal_pseudo_parents$ParentGene[i],
                  normal_expr, basal_expr, "Basal-like")
})
# Combine the list of results into a single data frame
basal_results <- bind_rows(basal_results_list)

# --- Combine and Adjust Results ---
# Merge HER2 and Basal results
all_results <- bind_rows(her2_results, basal_results)
# Apply Benjamini-Hochberg (FDR) correction to the p-values
all_results$Adj_P_value <- p.adjust(all_results$P_value, method = "BH")

# Save full results to CSV
write.csv(all_results, "Pseudogene_Decoupling_Analysis_Results.csv", row.names = FALSE)
# Print results to console for quick review
print(all_results)

# -------------------------------------------------------
# 5. Visualization
# -------------------------------------------------------
# Create a bar plot to visualize the change in correlation (Delta_r)
p <- ggplot(all_results, aes(x = reorder(Pseudogene, Delta_r), y = Delta_r, fill = Subtype)) +
  geom_col() + # Draw bars
  # Add an asterisk (*) for statistically significant pairs (Adj_P < 0.05)
  geom_text(aes(label = ifelse(Adj_P_value < 0.05, "*", "")), 
            vjust = 0.7, hjust = -0.3, size = 6, fontface = "bold") +
  coord_flip() + # Flip to make it a horizontal bar plot for easier reading of labels
  labs(
    title = "Pseudogene–Parent Decoupling Analysis (Δr)",
    subtitle = "Significant decoupling or reinforcement marked with *",
    x = "Pseudogene",
    y = "Δr (Correlation in Tumor - Correlation in Normal)"
  ) +
  # Add a dashed line at 0 to separate reinforced (+) from decoupled (-)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

# Save the plot
ggsave("Pseudogene_Decoupling_Plot.png", plot = p, width = 8, height = 6, dpi = 300)
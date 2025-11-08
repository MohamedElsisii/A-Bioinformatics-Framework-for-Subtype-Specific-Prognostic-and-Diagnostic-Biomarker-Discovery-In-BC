# ==================================
# Script: Code S5 - Receiver Operating Characteristic Curve and Area Under the Curve Analysis.R
#
# Purpose:
# This script evaluates the diagnostic performance of candidate biomarkers
# for distinguishing a specific breast cancer subtype from all other subtypes.
# It calculates the Area Under the Receiver Operating Characteristic Curve (AUC)
# using a rigorous repeated k-fold cross-validation (CV) approach.
#
# Key Outputs:
# 1. Per-subtype summary CSVs containing diagnostic metrics for all tested genes.
# 2. A combined master summary CSV for all subtypes.
# 3. Individual ROC plots (PNG format) for each gene that passes the CV threshold.
# ==================================

# Load necessary libraries, suppressing startup messages for cleaner output
suppressPackageStartupMessages({
  library(pROC)      # For ROC curve analysis and AUC calculation
  library(dplyr)     # For data manipulation (filter, mutate, arrange, etc.)
  library(readr)     # For efficient reading/writing of CSV files
  library(ggplot2)   # For creating and saving ROC plots
})

# -------------------------
# User params
# -------------------------
# Directory where individual ROC plots will be saved
plot_dir <- "roc_cv10_plots"
# Create the directory if it doesn't exist (recursive = TRUE handles parent dirs)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Minimum number of positive (subtype of interest) and negative (all other subtypes)
# samples required to run cross-validation. This prevents unstable results from
# too few samples. Adjusted based on typical subtype prevalence in TCGA.
thresholds <- list(
  LumA  = list(min_pos = 15, min_neg = 15),
  LumB  = list(min_pos = 12, min_neg = 12),
  Basal = list(min_pos = 10, min_neg = 10),
  Her2  = list(min_pos = 7,  min_neg = 7)
)

# Cross-validation settings:
k_folds    <- 10  # Number of folds for CV (standard is 10)
cv_repeats <- 3   # Number of times to repeat the entire 10-fold CV process
set.seed(42)      # Set seed for reproducibility of random splits
log_transform <- FALSE # Whether to log2-transform expression data before analysis

# -------------------------
# Helpers
# -------------------------
# Helper to create safe file names by replacing non-alphanumeric chars with underscores
safe_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# Helper to standardize sample IDs (e.g., TCGA.A1.A0SD -> TCGA-A1-A0SD)
norm_samples <- function(x) gsub("\\.", "-", as.character(x))

#' @title Robust File Loader
#' @description Tries to load an R object from a file, handling both .rds
#' (single object) and .rda/.RData (environment) formats.
#' @param path Path to the file
#' @return The loaded R object
load_rds_get <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  # First, try reading as an RDS file
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(obj)) return(obj)
  
  # If that fails, try loading as an RDA file into a new environment
  env <- new.env()
  nm <- tryCatch(load(path, envir = env), error = function(e) character(0))
  # Error handling for empty or multiple-object RDA files
  if (length(nm) == 0) stop("No readable object found in ", path)
  if (length(nm) > 1) warning("Multiple objects in ", path, " — returning first: ", nm[1])
  # Return the first object loaded into the environment
  get(nm[1], envir = env)
}

#' @title Plot and Save ROC Curve
#' @description Generates a standard ROC curve plot using ggplot2 and saves it.
#' Note: This plots the ROC for the *entire* dataset, not one specific CV fold,
#' to give an overall visual representation.
#' @param labels True class labels (0 or 1)
#' @param scores Predicted scores (gene expression values)
#' @param out_file File path for the saved plot
#' @param title Plot title
#' @param subtitle Plot subtitle (usually contains AUC stats)
#' @return TRUE if successful, FALSE otherwise
plot_roc_curve <- function(labels, scores, out_file, title = NULL, subtitle = NULL) {
  # Calculate ROC object quietly
  roc_obj <- tryCatch(pROC::roc(labels, scores, quiet = TRUE), error = function(e) NULL)
  if (is.null(roc_obj)) return(FALSE)
  
  # Create the plot using ggroc
  p <- ggroc(roc_obj) +
    # Add a diagonal dashed line representing random guessing (AUC = 0.5)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(title = title, subtitle = subtitle, x = "False Positive Rate", y = "True Positive Rate") +
    theme_classic() +
    theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9))
  
  # Save the plot to file
  ggsave(out_file, plot = p, width = 6, height = 6, dpi = 300)
  TRUE
}

#' @title Create Stratified Folds
#' @description Splits positive and negative samples separately to ensure
#' each fold has a representative ratio of classes (stratification).
#' @param pos_idx Indices of positive samples
#' @param neg_idx Indices of negative samples
#' @param k Number of folds
#' @return A list of k vectors, each containing indices for one fold
create_stratified_folds <- function(pos_idx, neg_idx, k) {
  # Randomly assign each positive/negative sample to one of k folds
  folds_pos <- split(sample(pos_idx), rep(1:k, length.out = length(pos_idx)))
  folds_neg <- split(sample(neg_idx), rep(1:k, length.out = length(neg_idx)))
  # Combine positive and negative indices for each fold
  lapply(seq_len(k), function(i) c(folds_pos[[i]], folds_neg[[i]]))
}

#' @title Calculate Cross-Validated AUC
#' @description Performs repeated, stratified k-fold CV to estimate AUC
#' mean and standard deviation.
#' @param scores Gene expression values
#' @param labels True class labels (0 for other, 1 for subtype of interest)
#' @param k Number of folds (default 10)
#' @param repeats Number of full CV repetitions (default 3)
#' @return A list containing mean AUC, SD of AUC, number of successful folds, and all raw AUCs
cv_auc <- function(scores, labels, k = 10, repeats = 3) {
  labels <- as.integer(labels)
  pos_idx <- which(labels == 1)
  neg_idx <- which(labels == 0)
  
  # Basic check: need at least 2 samples of each class to do ANYTHING meaningful
  if (length(pos_idx) < 2 || length(neg_idx) < 2) {
    return(list(mean = NA_real_, sd = NA_real_, n_folds = 0L, aucs = numeric(0)))
  }
  
  all_aucs <- vector("numeric", 0)
  # Outer loop for repeats
  for (rep_i in seq_len(repeats)) {
    # Create new random stratified folds for this repeat
    folds <- create_stratified_folds(pos_idx, neg_idx, k)
    # Inner loop for k folds
    for (fold_idx in seq_along(folds)) {
      test_idx <- folds[[fold_idx]]
      # Skip fold if it's too small or only has one class (can happen with very small datasets)
      if (length(test_idx) < 4) next
      if (length(unique(labels[test_idx])) < 2) next
      
      # Calculate AUC for this specific test fold
      roc_obj <- tryCatch(pROC::roc(labels[test_idx], scores[test_idx], quiet = TRUE), error = function(e) NULL)
      if (!is.null(roc_obj)) {
        all_aucs <- c(all_aucs, as.numeric(pROC::auc(roc_obj)))
      }
    }
  }
  
  # If no folds were successful, return NAs
  if (length(all_aucs) == 0) {
    return(list(mean = NA_real_, sd = NA_real_, n_folds = 0L, aucs = numeric(0)))
  }
  
  # Return summary statistics
  list(mean = mean(all_aucs), sd = if (length(all_aucs) > 1) sd(all_aucs) else NA_real_,
       n_folds = length(all_aucs), aucs = all_aucs)
}

# -------------------------
# Load BRCA expression matrix
# -------------------------
message("Loading BRCA expression matrix...")
BRCA_obj <- load_rds_get("BRCA.rda")
# Handle cases where the loaded object might be wrapped in a list
if (is.list(BRCA_obj) && length(BRCA_obj) == 1 && (is.matrix(BRCA_obj[[1]]) || is.data.frame(BRCA_obj[[1]]))) {
  BRCA <- BRCA_obj[[1]]
} else {
  BRCA <- BRCA_obj
}
# Ensure it's a numeric matrix for fast indexing
if (!is.matrix(BRCA)) BRCA <- as.matrix(BRCA)
mode(BRCA) <- "numeric"
message("BRCA dims: ", nrow(BRCA), " genes x ", ncol(BRCA), " samples")

# -------------------------
# Load candidate lists (all *_Candidates.rda)
# -------------------------
# Automatically find all candidate gene files in the working directory
cand_files <- list.files(pattern = "_Candidates\\.rda$", ignore.case = TRUE)
candidate_gene_list <- list()

if (length(cand_files) == 0) {
  message("No candidate files found (pattern: _Candidates.rda).")
}

# Load each candidate file and store in a named list (key = subtype prefix)
for (f in cand_files) {
  subname <- sub("_Candidates\\.rda$", "", f, ignore.case = TRUE)
  obj <- load_rds_get(f)
  # Ensure genes are unique character vectors
  candidate_gene_list[[subname]] <- unique(as.character(unlist(obj)))
}
message("Loaded candidate lists for subtypes: ", paste(names(candidate_gene_list), collapse = ", "))

# -------------------------
# Load sample lists (all *_Samples.rda) and build metadata
# -------------------------
# Automatically find all sample list files to determine true class labels
sample_files <- list.files(pattern = "_Samples\\.rda$", ignore.case = TRUE)
meta_list <- list()

for (f in sample_files) {
  subname <- sub("_Samples\\.rda$", "", f, ignore.case = TRUE)
  samp <- load_rds_get(f)
  # Standardize sample IDs to match the expression matrix
  samp <- norm_samples(as.character(unlist(samp)))
  # Create a mini-metadata frame for this subtype
  meta_list[[subname]] <- data.frame(SampleID = samp, Subtype = subname, stringsAsFactors = FALSE)
}

if (length(meta_list) == 0) stop("No sample files found matching *_Samples.rda")
# Combine all subtype sample lists into one master metadata table
metadata <- bind_rows(meta_list)

# --- Align Metadata with Expression Matrix ---
# Keep only samples that are actually present in the BRCA matrix
count_cols_trim <- norm_samples(colnames(BRCA))
metadata <- metadata %>% filter(SampleID %in% count_cols_trim)
metadata$Subtype <- factor(metadata$Subtype, levels = unique(metadata$Subtype))

# Reorder the BRCA matrix columns to match the filtered metadata exactly
keep_cols <- unique(metadata$SampleID)
BRCA <- BRCA[, match(keep_cols, count_cols_trim), drop = FALSE]
colnames(BRCA) <- keep_cols
message("Filtered BRCA to ", ncol(BRCA), " samples available in sample lists/metadata")

# -------------------------
# Main loop: subtype -> genes
# -------------------------
results <- list()

# Iterate through each subtype that has a candidate gene list
for (sub in names(candidate_gene_list)) {
  message("\nProcessing subtype: ", sub)
  genes <- unique(candidate_gene_list[[sub]])
  
  # Skip if this subtype has no sample metadata
  if (!sub %in% metadata$Subtype) {
    warning("No sample metadata for subtype: ", sub, " — skipping.")
    next
  }
  
  # Define Positive samples (this subtype) and Negative samples (ALL other subtypes)
  pos_samps <- metadata$SampleID[metadata$Subtype == sub]
  neg_samps <- setdiff(colnames(BRCA), pos_samps)
  
  if (length(pos_samps) == 0) {
    warning("No positive samples found for subtype ", sub, " — skipping.")
    next
  }
  
  # Get sample size thresholds for this subtype, use default if not defined
  thr <- thresholds[[sub]]
  if (is.null(thr)) {
    thr <- list(min_pos = 10, min_neg = 10)
    message("  No thresholds defined for ", sub, " — using defaults min_pos=10,min_neg=10")
  }
  
  # Pre-calculate column indices for speed inside the gene loop
  pos_idx_all <- match(pos_samps, colnames(BRCA))
  neg_idx_all <- match(neg_samps, colnames(BRCA))
  
  subtype_results <- vector("list", length(genes))
  ii <- 0L
  
  # --- Gene Loop ---
  for (g in genes) {
    ii <- ii + 1L
    
    # Handle case where a candidate gene is missing from the expression matrix
    if (!g %in% rownames(BRCA)) {
      subtype_results[[ii]] <- data.frame(
        Subtype = sub, Gene = g, InMatrix = FALSE, n_pos = NA_integer_, n_neg = NA_integer_,
        Mean_Pos = NA_real_, Mean_Neg = NA_real_, Mean_Diff = NA_real_,
        Wilcox_p = NA_real_, AUC_CV = NA_real_, AUC_CV_SD = NA_real_, AUC_nfolds = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Extract expression values for positive and negative classes
    expr_pos <- as.numeric(BRCA[g, pos_idx_all])
    expr_neg <- as.numeric(BRCA[g, neg_idx_all])
    
    # Optional log transformation
    if (log_transform) {
      expr_pos <- log2(expr_pos + 1)
      expr_neg <- log2(expr_neg + 1)
    }
    
    # Basic descriptive statistics
    npos <- sum(!is.na(expr_pos))
    nneg <- sum(!is.na(expr_neg))
    mean_pos <- if (npos > 0) mean(expr_pos, na.rm = TRUE) else NA_real_
    mean_neg <- if (nneg > 0) mean(expr_neg, na.rm = TRUE) else NA_real_
    mean_diff <- if (!is.na(mean_pos) && !is.na(mean_neg)) mean_pos - mean_neg else NA_real_
    
    # Wilcoxon rank-sum test (non-parametric comparison of means)
    w_p <- tryCatch({
      if (npos >= 1 && nneg >= 1) wilcox.test(expr_pos, expr_neg, alternative = "two.sided")$p.value else NA_real_
    }, error = function(e) NA_real_)
    
    # Prepare data for ROC analysis
    labels <- c(rep(1L, length(expr_pos)), rep(0L, length(expr_neg)))
    scores <- c(expr_pos, expr_neg)
    ok <- !is.na(scores) & !is.na(labels)
    
    auc_cv_mean <- NA_real_; auc_cv_sd <- NA_real_; auc_nfolds <- 0L
    
    # Check if we have enough samples to run rigorous CV
    if (sum(ok & labels == 1L) >= thr$min_pos && sum(ok & labels == 0L) >= thr$min_neg) {
      # Run repeated, stratified k-fold CV
      cv <- cv_auc(scores[ok], labels[ok], k = k_folds, repeats = cv_repeats)
      auc_cv_mean <- cv$mean; auc_cv_sd <- cv$sd; auc_nfolds <- cv$n_folds
      
      # If CV was successful (at least one fold worked), save a representative ROC plot
      if (auc_nfolds > 0) {
        fname <- file.path(plot_dir, paste0("ROC_", safe_name(sub), "_", safe_name(g), ".png"))
        # Create subtitle with CV statistics
        subtitle <- paste0("AUC=", ifelse(is.na(auc_cv_mean), "NA", formatC(auc_cv_mean, digits = 3)),
                           ", sd=", ifelse(is.na(auc_cv_sd), "NA", formatC(auc_cv_sd, digits = 3)),
                           ", folds=", auc_nfolds)
        # Plot using the full dataset for visualization
        tryCatch(plot_roc_curve(labels[ok], scores[ok], fname, title = paste0(g, " | ", sub), subtitle = subtitle),
                 error = function(e) NULL)
      }
    } else {
      message("  Skipping CV for ", g, " in ", sub, " (pos=", sum(ok & labels==1), " neg=", sum(ok & labels==0), ")")
    }
    
    # Store results for this gene
    subtype_results[[ii]] <- data.frame(
      Subtype = sub, Gene = g, InMatrix = TRUE, n_pos = npos, n_neg = nneg,
      Mean_Pos = mean_pos, Mean_Neg = mean_neg, Mean_Diff = mean_diff,
      Wilcox_p = w_p, AUC_CV = auc_cv_mean, AUC_CV_SD = auc_cv_sd, AUC_nfolds = auc_nfolds,
      stringsAsFactors = FALSE
    )
  } # End gene loop
  
  # Combine all gene results for this subtype
  res_df_sub <- bind_rows(subtype_results)
  
  # Apply Benjamini-Hochberg (FDR) correction to Wilcoxon p-values
  res_df_sub <- res_df_sub %>%
    mutate(Wilcox_p_num = ifelse(is.na(Wilcox_p), 1, as.numeric(Wilcox_p)),
           Wilcox_p_adj = p.adjust(Wilcox_p_num, method = "BH")) %>%
    arrange(desc(AUC_CV)) # Sort by best AUC first
  
  # Save subtype-specific summary CSV
  outname <- paste0("ROC_", safe_name(sub), "_summary.csv")
  readr::write_csv(res_df_sub, outname)
  message("Wrote ", outname)
  
  # Add to master results list
  results[[sub]] <- res_df_sub
} # End subtype loop

# -------------------------
# Combined summary
# -------------------------
# Combine all subtype results into one master CSV file
if (length(results) > 0) {
  res_df_all <- bind_rows(results)
  res_df_all <- res_df_all %>% arrange(Subtype, desc(AUC_CV))
  readr::write_csv(res_df_all, "ROC_summary_by_subtype.csv")
  message("Wrote ROC_summary_by_subtype.csv")
} else {
  message("No results to save.")
}
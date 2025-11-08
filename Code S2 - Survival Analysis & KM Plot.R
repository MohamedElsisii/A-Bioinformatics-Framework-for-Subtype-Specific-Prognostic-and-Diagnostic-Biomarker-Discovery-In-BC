# Survival_analysis_LumA_optimized.R
# Optimized Kaplan-Meier survival loop for a list of candidate genes (LumA).
#  - Pre-validates gene presence and expression to avoid wasted iterations
#  - Matches sample IDs reliably (uses TCGA submitter trimming)
#  - Robust error handling and informative logging
#  - Saves per-gene KM plots and a consolidated CSV with BH-adjusted p-values

suppressPackageStartupMessages({
  library(survminer)
  library(survival)
  library(tidyverse)
  library(DESeq2)
})

# -------------------------
# User options
# -------------------------
genes_file      <- "LumA_Significant.rda"
samples_file    <- "LumA_Samples.rda"
counts_file     <- "BRCA.rds"   # raw counts or SummarizedExperiment
clinical_file   <- "Clinical_Data.csv"
out_dir         <- "survival_plots"
min_samples     <- 10     # minimum samples required after join
plot_width      <- 7
plot_height     <- 6

dir.create(out_dir, showWarnings = FALSE)

# -------------------------
# Helpers
# -------------------------
log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
trim_submitter <- function(x) sub("^([^-]+-[^-]+-[^-]+).*$", "\\1", x)

# Safe read vector rds
read_char_rds <- function(path) {
  obj <- readRDS(path)
  as.character(unique(unlist(obj)))
}

# -------------------------
# Load inputs
# -------------------------
log_msg("Loading genes and sample lists...")
Genes <- read_char_rds(genes_file)
log_msg("Candidate genes: ", length(Genes))

LumASamples <- read_char_rds(samples_file) %>% trim_submitter()
log_msg("LumA sample IDs (trimmed): ", length(LumASamples))

log_msg("Loading counts...")
counts_obj <- readRDS(counts_file)
counts_mat <- if (inherits(counts_obj, "SummarizedExperiment")) assay(counts_obj) else as.matrix(counts_obj)
col_submitters <- trim_submitter(colnames(counts_mat))
log_msg("Counts matrix: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " samples")

# Subset to LumA samples that exist in counts
keep_cols_idx <- which(col_submitters %in% LumASamples)
if (length(keep_cols_idx) == 0) stop("No LumA samples found in counts matrix after trimming IDs.")
counts_LumA <- counts_mat[, keep_cols_idx, drop = FALSE]
colnames(counts_LumA) <- trim_submitter(colnames(counts_LumA))
log_msg("Counts subset: ", ncol(counts_LumA), " LumA columns kept")

# Load clinical and harmonize
log_msg("Loading clinical data...")
clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
if (!"submitter_id" %in% colnames(clinical)) {
  stop("Clinical file must contain a 'submitter_id' column")
}
clinical$submitter_id_trim <- trim_submitter(clinical$submitter_id)
clinical$deceased <- ifelse(tolower(clinical$vital_status) == "alive", 0, 1)
clinical$overall_survival <- as.numeric(coalesce(clinical$days_to_death, clinical$days_to_last_follow_up))

# Keep only LumA clinical records and complete cases
clinical_LumA <- clinical %>% filter(submitter_id_trim %in% colnames(counts_LumA)) %>%
  select(submitter_id_trim, deceased, overall_survival) %>% drop_na()
log_msg("Clinical records matched: ", nrow(clinical_LumA))
if (nrow(clinical_LumA) < min_samples) stop("Not enough matched clinical records for LumA (n < min_samples)")

# Make sure order of columns in counts matches clinical rows
counts_LumA <- counts_LumA[, clinical_LumA$submitter_id_trim, drop = FALSE]

# -------------------------
# VST transformation
# -------------------------
log_msg("Filtering low-count genes and applying VST...")
keep_genes <- rowSums(counts_LumA) >= 10
counts_LumA_filtered <- counts_LumA[keep_genes, , drop = FALSE]
log_msg("Genes kept for VST: ", nrow(counts_LumA_filtered))

dds <- DESeqDataSetFromMatrix(countData = round(counts_LumA_filtered),
                              colData = data.frame(row.names = clinical_LumA$submitter_id_trim),
                              design = ~ 1)
vsd <- vst(dds, blind = FALSE)
counts_vst <- assay(vsd)

# Limit Genes to those present and expressed
Genes_present <- intersect(Genes, rownames(counts_vst))
Genes_present <- Genes_present[rowSums(counts_vst[Genes_present, , drop = FALSE]) > 0]
log_msg("Candidate genes present & expressed: ", length(Genes_present))
if (length(Genes_present) == 0) stop("No candidate genes available in VST matrix")

# -------------------------
# Per-gene survival function
# -------------------------
run_gene_survival <- function(gene) {
  tryCatch({
    expr_vals <- counts_vst[gene, ]
    df <- tibble(submitter_id_trim = colnames(counts_vst), expr = as.numeric(expr_vals)) %>%
      left_join(clinical_LumA, by = "submitter_id_trim") %>% drop_na()
    if (nrow(df) < min_samples) return(NULL)
    
    # Determine optimal cutpoint
    sc <- surv_cutpoint(df, time = "overall_survival", event = "deceased", variables = "expr")
    cp <- sc$cutpoint$cutpoint
    if (is.null(cp) || is.na(cp)) return(NULL)
    
    df <- df %>% mutate(strata = ifelse(expr > cp, "HIGH", "LOW"))
    if (length(unique(df$strata)) < 2) return(NULL)
    
    fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = df)
    pval <- tryCatch({ surv_pvalue(fit)$pval }, error = function(e) NA_real_)
    
    # Save KM plot (with risk table)
    plot_obj <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                           title = paste0("Survival: ", gene),
                           subtitle = paste0("Cutoff=", round(cp, 3)),
                           legend.title = gene)
    ggsave(filename = file.path(out_dir, paste0("Survival_KM_", gene, ".png")), plot = plot_obj$plot,
           width = plot_width, height = plot_height)
    
    tibble(gene = gene, n = nrow(df), cutoff = cp, logrank_p_value = pval)
  }, error = function(e) {
    message("Error for gene ", gene, ": ", e$message)
    NULL
  })
}

# Collect results
res_df <- bind_rows(res_list)
if (nrow(res_df) == 0) {
  log_msg("No valid survival results produced.")
} else {
  res_df <- res_df %>% mutate(logrank_adj_p_value = p.adjust(logrank_p_value, method = "BH"))
  write.csv(res_df, file = file.path(out_dir, "LumA_survival_results.csv"), row.names = FALSE)
  log_msg("Saved results for ", nrow(res_df), " genes to ", file.path(out_dir, "LumA_survival_results.csv"))
}

log_msg("Finished.")

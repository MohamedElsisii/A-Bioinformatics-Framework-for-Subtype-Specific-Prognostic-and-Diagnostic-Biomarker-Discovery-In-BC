# -----------------------------------------------------------------------------
# GENE PROGNOSTIC VALIDATION SCRIPT
# -----------------------------------------------------------------------------
# This script performs a comprehensive validation of a single gene's prognostic
# value in a cancer subtype cohort. The analysis includes:
# 1. Checking for confounding with clinical variables (age, stage).
# 2. Visualizing survival differences using Kaplan-Meier plots.
# 3. Assessing independent prognostic significance with a multivariate Cox model.
# 4. Quantifying standalone predictive accuracy using ROC/AUC analysis.
# -----------------------------------------------------------------------------

# --- 0. SETUP: LOAD LIBRARIES ---

# Suppress startup messages for a cleaner console output
suppressPackageStartupMessages({
  library(survival)    # Core package for survival analysis (Surv, coxph, survfit)
  library(survminer)   # For creating publication-ready survival plots (ggsurvplot)
  library(DESeq2)      # Used for expression data handling (if expr_vst is a DESeq object)
  library(tidyverse)   # A collection of R packages for data science (dplyr, ggplot2)
  library(broom)       # Helps to tidy model output into data frames
  library(pROC)        # For Receiver Operating Characteristic (ROC) curve analysis
})

# --- 1. INITIALIZATION: DEFINE GENE AND PREPARE DATA ---

# Define the gene symbol to be validated in this run.
gene_to_validate <- "DLL3"

# Create the master data frame for the analysis.
# This combines the gene's expression data with relevant clinical information.
df_gene <- data.frame(
  expr = expr_vst[gene_to_validate, ],      # Gene expression values (VST normalized)
  age = clinical_subtype_multi$age_at_diagnosis,  # Patient age at diagnosis
  stage = clinical_subtype_multi$ajcc_pathologic_stage, # Pathological stage
  OS_time = clinical_subtype_multi$overall_survival,   # Overall survival time (e.g., in days or months)
  OS_status = clinical_subtype_multi$deceased         # Survival status (e.g., 0 = alive, 1 = deceased)
) %>%
  # IMPORTANT: Filter out any rows with missing data (NA) to ensure model compatibility.
  filter(!is.na(expr), !is.na(age), !is.na(stage), !is.na(OS_time), !is.na(OS_status))


# --- 2. CONFOUNDING ANALYSIS: CHECK ASSOCIATION WITH CLINICAL FACTORS ---
message("Step 1: Checking for potential confounding variables...")

# Test the correlation between gene expression and patient age.
# Using Spearman's rank correlation because expression data is often not normally distributed.
cor_age_test <- cor.test(df_gene$expr, df_gene$age, method = "spearman")
print(cor_age_test)

# Test for an association between gene expression and pathological stage.
# Using Kruskal-Wallis test, a non-parametric alternative to ANOVA, suitable for categorical variables.
kruskal_stage_test <- kruskal.test(expr ~ stage, data = df_gene)
print(kruskal_stage_test)


# --- 3. SURVIVAL VISUALIZATION: KAPLAN-MEIER PLOT ---
message("Step 2: Generating Kaplan-Meier plot to visualize survival...")

# Determine the optimal expression cutpoint to stratify patients into "high" and "low" groups.
# This maximizes the survival difference between the two groups.
km_cut <- surv_cutpoint(df_gene, time = "OS_time", event = "OS_status", variables = "expr")

# Apply the calculated cutpoint to categorize the patients.
km_cat <- surv_categorize(km_cut)

# Fit the Kaplan-Meier survival model.
fit <- survfit(Surv(OS_time, OS_status) ~ expr, data = km_cat)

# Create and print the Kaplan-Meier plot.
# This plot visually represents the survival probability over time for the high vs. low expression groups.
km_plot <- ggsurvplot(
  fit,
  data = km_cat,
  pval = TRUE,             # Display the log-rank p-value on the plot
  conf.int = TRUE,         # Show confidence intervals for the survival curves
  risk.table = TRUE,       # Add a table showing the number of patients at risk over time
  title = paste("Kaplan-Meier Plot for", gene_to_validate)
)
print(km_plot)


# --- 4. MULTIVARIATE ANALYSIS: COX PROPORTIONAL HAZARDS MODEL ---
message("Step 3: Running multivariate Cox model to assess independent prognostic value...")

# Merge the dichotomized expression categories ("high"/"low") back into the main data frame.
df_gene_final <- cbind(df_gene, Expression = km_cat$expr)

# Fit the multivariate Cox proportional hazards model.
# This model evaluates the gene's prognostic significance while adjusting for age and stage.
cox_model <- coxph(
  Surv(OS_time, OS_status) ~ Expression + age + stage,
  data = df_gene_final
)
print(summary(cox_model))

# Visualize the Cox model results with a forest plot.
# This plot shows the Hazard Ratio (HR) and 95% Confidence Interval for each variable.
ggforest(cox_model, data = df_gene_final,
         main = paste(gene_to_validate, "Multivariate Cox Model"),
         cpositions = c(0.02, 0.22, 0.4),
         fontsize = 1.0)

# Check the proportional hazards (PH) assumption, a key requirement for the Cox model.
# A p-value < 0.05 for any variable suggests its effect on survival changes over time, violating the assumption.
ph_test <- cox.zph(cox_model)
print(ph_test)

# --- 4a. MODEL CORRECTION (if PH assumption is violated) ---
# If the `ph_test` showed a violation for 'stage' (a categorical variable), a stratified Cox model is used.
# Stratification adjusts for the violating variable without estimating its HR, thus satisfying the assumption.
cox_model_stratified <- coxph(
  Surv(OS_time, OS_status) ~ Expression + age + strata(stage),
  data = df_gene_final
)

# Print the summary of the new, corrected stratified model.
print(summary(cox_model_stratified))

# Re-check the PH assumption on the stratified model to confirm the violation is resolved.
# The 'stage' variable is no longer tested as it is now a stratum.
ph_test_stratified <- cox.zph(cox_model_stratified)
print(ph_test_stratified)


# --- 5. PREDICTIVE PERFORMANCE: ROC/AUC ANALYSIS ---
message("Step 4: Quantifying standalone predictive performance with ROC analysis...")

# Create a ROC object using the continuous expression values to predict patient status.
# This assesses how well the gene alone can discriminate between cases (deceased) and controls (alive).
roc_obj <- roc(response = df_gene$OS_status, predictor = df_gene$expr)

# Print the Area Under the Curve (AUC) and its 95% Confidence Interval (CI).
# AUC = 0.5 -> no predictive ability (random chance)
# AUC = 1.0 -> perfect predictive ability
print(roc_obj)
ci(roc_obj)

# Plot the ROC curve.
# This visualizes the trade-off between the true positive rate and the false positive rate.
plot(roc_obj, main = paste("ROC Curve for", gene_to_validate), print.auc = TRUE, col = "blue")


# --- 6. CONCLUSION ---
message("--- Validation for ", gene_to_validate, " is complete. ---")
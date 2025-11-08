# A-Bioinformatics-Framework-for-Subtype-Specific-Prognostic-and-Diagnostic-Biomarker-Discovery-In-BC
![Copy of Graphical Abstract - Bioinformatic Analysis Reveals Novel Subtype-Specific Prognostic and Diagnostic Biomarkers in Breast Cancer _1](https://github.com/user-attachments/assets/f2f5a393-951d-4927-bfd7-4c97cdd50aea)

# Overview
This repository contains the complete bioinformatics workflow associated with our study, "Resolving Breast Cancer Heterogeneity: A Bioinformatics Framework for Subtype-Specific Prognostic and Diagnostic Biomarker Discovery".

Breast cancer (BC) is a highly heterogeneous disease. While molecular subtyping using PAM50 classifier has improved clinical precision, significant variability in patient outcomes persists within each subtype. This framework was developed to bridge this gap by identifying robust, subtype-specific biomarkers that can serve as precise tools for both prognosis and early diagnosis.

Our multi-step analysis integrates differential expression, survival modeling, multivariate Cox regression, and machine learning-based diagnostic evaluation to propose a refined panel of 194 candidate biomarkers across four major molecular subtypes (Luminal A, Luminal B, Basal-like, HER2-enriched).

# Workflow

The analysis pipeline implemented in the provided scripts follows this systematic framework:

1- **Data Acquisition & Subtyping**: Processing TCGA-BRCA RNA-seq data and classifying tumors using the PAM50 classifier.

2- **Differential Expression Analysis**: Identifying unique Differentially Expressed Genes (DEGs) for each subtype against normal tissue using edgeR and limma.

3- **Prognostic Screening**: Kaplan-Meier survival analysis to identify initial prognostic candidates.

4- **Robust Validation**: Univariate and multivariate Cox proportional hazards regression, adjusted for clinical covariates (age, stage).

5- **Diagnostic Evaluation**: ROC curve analysis with 10-fold cross-validation to assess discriminatory power.

6- **Advanced Functional Characterization**:

  - Correlation with PAM50 signature genes.

  - Pseudogene–parent decoupling analysis.

  - Pathway-level functional rewiring analysis using MSigDB Hallmark pathways.

# Repository Structure

This repository hosts the 9 core R scripts used to execute the entire analysis pipeline.

## Script Descriptions

### Code S1 - Differentially Expressed Genes.R
Performs DEG analysis using **edgeR** and **limma** to find subtype-specific signatures.

### Code S2 - Survival Analysis & KM Plot.R
Screens DEGs for prognostic potential using **Kaplan-Meier estimates** and generates KM plots.

### Code S3 - Univariate & Multivariate Cox.R
Builds univariate and multivariate **Cox models** to assess independent prognostic value.

### Code S4 - Prognostic Biomarkers Validation.R
Validates prognostic markers, checking for **proportional hazards** and adjusting models if necessary.

### Code S5 - Receiver Operating Characteristic Curve and Area Under the Curve Analysis.R
Evaluates diagnostic performance using **ROC curves** and 10-fold cross-validation.

### Code S6 - Correlation with PAM50 Genes.R
Assesses correlation between candidate markers and standard **PAM50 signature genes**.

### Code S7 - Pseudogenes Decoupling Analysis.R
Analyzes regulatory divergence between **pseudogenes** and their **parent genes**.

### Code S8 - Functional Rewiring of lncRNAs.R
Performs **GSVA** to identify subtype-specific pathway functional rewiring for **lncRNAs**.

### Code S9 - Functional Rewiring of Protein-Coding Genes.R
Performs **GSVA** to identify subtype-specific pathway functional rewiring for **protein-coding genes**.

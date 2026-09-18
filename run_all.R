################################################################################
################ MASTER RUN SCRIPT: run_all.R ###################################
#################################################################################
#* by: Isabel Santonja
#* last modified: 18/09/2026
#* description: This script reproduces the entire analysis pipeline for the
#* manuscript on shift schedules and sleep using the data from the online
#* survey of the SHIFT2HEALTH project
# ==========================================

# 1. Load required library for robust relative paths ####
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
library(rmarkdown)

cat("Starting full analysis pipeline...\n")

# 2. Run data processing scripts ####
cat("\nStep 1: Cleaning raw data...\n")
## 2.1 Edit long schedules dataset
cat("\nStep 1.1: Edits in long dataset...\n")
source(here("scripts", "01_schedules_cleaning.R"))

cat("\nStep 1.2: Initial data analysis and wrangling...\n")
render(input = here("reports", "02_initial_data_analysis_and_wrangling.Rmd"), 
       output_format = "html_document")

cat("\nStep 2: Fit VMF cluster models...\n")
render(input = here("reports", "03_vmf_shift_clusters.Rmd"),
       output_format = "html_document")

cat("\nStep 3: Shift-specific sleep duration and quality (Figure 1)...\n")
source(here("scripts", "04_shift_specific_sleep.R"))

cat("\nStep 4: Prepare participant-level datasets (level2_data, brm_data)...\n")
source(here("scripts", "05_prepare_model_data.R"))

cat("\nStep 5: Sample description (Sup. Table S2)...\n")
source(here("scripts", "06_sample_description.R"))

cat("\nStep 6: DAG and derivation of adjustment sets (Sup. Figure S3)...\n")
source(here("scripts", "07_dag.R"))

cat("\nStep 7: Isotemporal substitution analysis (Figure 2, Sup. Fig. 5-6)...\n")
source(here("scripts", "08_isotemporal_substitution.R"))

cat("\nStep 8: Classic shift work metrics and sensitivity analyses (Tables 2-3, Sup. Tables S3-S8)...\n")
source(here("scripts", "09_classic_metrics_and_sensitivity.R"))

cat("\nPipeline complete! All outputs generated successfully.\n")

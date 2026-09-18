################################################################################
############### CLASSIC SHIFT WORK METRICS: TABLES 2-3 & SENSITIVITY ###########
################################################################################

#* by: Isabel Santonja
#* description: This script fits the Bayesian logistic regression models for
#* Tables 2-3 (classic night shift metrics vs chronic insomnia and short
#* average sleep duration; Models 1, 2 and the mutually-adjusted Model 4), and
#* the sensitivity analyses reported in Supplementary Tables S3-S8:
#* complete-case analysis (S3-S4), non-informative priors (S5-S6), and
#* stratification by chronotype (S7-S8).

# This script reads "data/brm_data.RDS", built by 05_prepare_model_data.R -
# see that script and the project README for details.

# Load necessary libraries and functions
library(tidyverse)
library(brms)
library(broom.mixed)
library(writexl)

source(here::here("functions", "tables.R"))
source(here::here("functions", "model_functions.R"))

# Load the data
data <- readRDS(here::here("data", "brm_data.RDS")) |>
  mutate(
    schedule_type = as.factor(case_when(
      number_of_shifts == 1 & !night_shift_reported
      ~ "Permanent day",
      number_of_shifts == 1 & night_shift_reported
      ~ "Permanent night",
      number_of_shifts > 1
      ~ "Rotating shifts"
    )),
    night_shift_freq_rel = 10 * oq_night_shift_freq_rel, # Scale by 10 for interpretability
    number_of_shifts_cat3 = factor(case_when(
      number_of_shifts >= 3 ~ "3+", # Merge 4th category for stratified analyses
      .default = number_of_shifts_cat
    )),
    chronotype_3cat = as.factor(chronotype_3cat),
    years_worked_ns_grouped = case_when(
      years_worked_ns_cat %in% c("less than one year", "1-5 years")
      ~ "up to 5 years",
      is.na(years_worked_ns_cat)
      ~ "0 years",
      years_worked_ns_cat %in% c(
        "16-25 years", "26-35 years", "36-40 years", "longer than 40 years"
      )
      ~ "16+ years",
      .default = years_worked_ns_cat
  )
  ) |>
  mutate(years_worked_ns_grouped = fct_relevel(
    years_worked_ns_grouped,
    "0 years", "up to 5 years", "6-15 years", "16+ years"
  ))

## Prevalences
cat_exposures <- c("night_shifts_updated", "years_worked_ns_grouped",
                    "schedule_type", "number_of_shifts_cat3",
                    "shift_rotation_speed")

for (var in cat_exposures) {
  print(prev_tab(data, chronic_insomnia, !!sym(var)))
  print(prev_tab(data, avg_short_sleep, !!sym(var)))
}


################################### TABLE 2 ####################################
# 1. Define the sets of covariates, priors and functions ####
## Define variables
### Exposures
main_exposures <- c(cat_exposures, "working_hours_scaled", "night_shift_freq_rel")

## Define the sets of covariates in a named list
#* Model 1: age and gender. Model 2: Model 1 + chronotype, country, education,
#* income satisfaction, marital status, BMI, alcohol and smoking (as in
#* Tables 2-3). Model 4 (mutually adjusted) uses the Model 2 covariates.
model_covariates_list <- list(
  "model1" = c("age_scaled", "gender_grouped"),
  "model2" = c("age_scaled", "gender_grouped", "chronotype", "country_grouped",
               "education_grouped", "income_satisfaction",
               "marital_status_grouped", "alcohol", "bmi_scaled",
               "smoking_status")
)

## Define base priors for the model parameters
#* exposure-specific priors will be defined inside the loop
base_ins_priors <- c(
  set_prior("normal(-2.2, 1)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b")
)

base_ssl_priors <- c(
  set_prior("normal(-2.7, 1)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b")
)

# 2. Fit insomnia models ####
ins_results <- run_exposure_models(
  data = data,
  outcome = "chronic_insomnia",
  main_exposures = main_exposures,
  model_covariates_list = model_covariates_list,
  base_priors = base_ins_priors,
  results_dir = "results/t2/insomnia",
  label = "t2_ins"
)

## 2.1 Mutually adjusted model (Model 4)
non_pred_exposures <- c(
  "schedule_type", "years_worked_ns_grouped",
  "working_hours_scaled")

fit_mut_adj_model(
  data = data,
  outcome = "chronic_insomnia",
  exposures = non_pred_exposures,
  covariates = model_covariates_list[["model2"]],
  priors = base_ins_priors,
  results_dir = "results/t2/insomnia",
  model_name = "ins_ma_4")

# 3. Fit short sleep models ####
ssl_results <- run_exposure_models(
  data = data,
  outcome = "avg_short_sleep",
  main_exposures = main_exposures,
  model_covariates_list = model_covariates_list,
  base_priors = base_ssl_priors,
  results_dir = "results/t2/ssl",
  label = "t2_ssl"
)

## 3.1 Mutually adjusted model (Model 4)
fit_mut_adj_model(
  data = data,
  outcome = "avg_short_sleep",
  exposures = non_pred_exposures,
  covariates = model_covariates_list[["model2"]],
  priors = base_ssl_priors,
  results_dir = "results/t2/ssl",
  model_name = "ssl_ma_4")

### Check generalized variance inflation factors in the mutually adjusted model
fit_lm <- glm(avg_short_sleep ~ schedule_type + years_worked_ns_log_scaled +
                working_hours_scaled + age_scaled + gender_grouped +
                country_grouped + education_grouped + income_satisfaction +
                marital_status_grouped,
              data = data, family = binomial)

car::vif(fit_lm)
#* GVIFs were <1.3 (no meaningful multicollinearity), i.e. the effect of years
#* worked on night shifts was largely independent of, rather than confounded
#* by, current schedule type.


############################ SENSITIVITY ANALYSES ##############################

# 1. Complete-case analysis, excluding missings in covariates (Tables S3-S4) ####
data_complete <- data |>
  drop_na(oq_bmi, country, oq_alc_per_wd, smoking_status) |>
  filter(smoking_status != "Unknown") |>
  mutate( # rescale continuous covariates
    working_hours_scaled = as.numeric(scale(working_hours)),
    years_worked_ns_log_scaled = as.numeric(scale(years_worked_ns_log)),
    age_scaled = as.numeric(scale(oq_age)),
    bmi_scaled = as.numeric(scale(oq_bmi))
  )
cat("Data after excluding missings:", nrow(data_complete), "rows\n")

ins_results_wo_missings <- run_exposure_models(
  data = data_complete,
  outcome = "chronic_insomnia",
  main_exposures = main_exposures,
  model_covariates_list = model_covariates_list,
  base_priors = base_ins_priors,
  results_dir = "results/t2_sens/insomnia_wo_missings",
  label = "t2_ins"
)
ssl_results_wo_missings <- run_exposure_models(
  data = data_complete,
  outcome = "avg_short_sleep",
  main_exposures = main_exposures,
  model_covariates_list = model_covariates_list,
  base_priors = base_ssl_priors,
  results_dir = "results/t2_sens/ssl_wo_missings",
  label = "t2_ssl"
)

## Mutually adjusted models (Model 4)
fit_mut_adj_model(
  data = data_complete,
  outcome = "chronic_insomnia",
  exposures = non_pred_exposures,
  covariates = model_covariates_list[["model2"]],
  priors = base_ins_priors,
  results_dir = "results/t2_sens/insomnia_wo_missings",
  model_name = "ins_ma_4")

fit_mut_adj_model(
  data = data_complete,
  outcome = "avg_short_sleep",
  exposures = non_pred_exposures,
  covariates = model_covariates_list[["model2"]],
  priors = base_ssl_priors,
  results_dir = "results/t2_sens/ssl_wo_missings",
  model_name = "ssl_ma_4")

rm(data_complete) # free up memory

# 2. Non-informative priors (Tables S5-S6) ####
wider_base_priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),
  set_prior("normal(0, 100)", class = "b")
  )

ins_results_wider_priors <- run_exposure_models(
  data = data,
  outcome = "chronic_insomnia",
  main_exposures = main_exposures,
  model_covariates_list = model_covariates_list,
  base_priors = wider_base_priors,
  exposure_prior_dist = "normal(0, 10)",
  m_treedepth = 15,
  results_dir = "results/t2_sens/insomnia_wider_priors",
  label = "t2_ins"
)
ssl_results_wider_priors <- run_exposure_models(
  data = data,
  outcome = "avg_short_sleep",
  main_exposures = main_exposures,
  model_covariates_list = model_covariates_list,
  base_priors = wider_base_priors,
  exposure_prior_dist = "normal(0, 10)",
  m_treedepth = 15,
  results_dir = "results/t2_sens/ssl_wider_priors",
  label = "t2_ssl"
  )

## Mutually adjusted models (Model 4)
fit_mut_adj_model(
  data = data,
  outcome = "chronic_insomnia",
  exposures = non_pred_exposures,
  covariates = model_covariates_list[["model2"]],
  priors = wider_base_priors,
  exposure_prior_dist = "normal(0, 10)",
  results_dir = "results/t2_sens/insomnia_wider_priors",
  model_name = "ins_ma_4")

fit_mut_adj_model(
  data = data,
  outcome = "avg_short_sleep",
  exposures = non_pred_exposures,
  covariates = model_covariates_list[["model2"]],
  priors = wider_base_priors,
  exposure_prior_dist = "normal(0, 10)",
  results_dir = "results/t2_sens/ssl_wider_priors",
  model_name = "ssl_ma_4")

############################ STRATIFIED ANALYSIS (Tables S7-S8) ################
## Adjustment covariates excluding chronotype itself (used as the stratifier)
strat_analysis_cov <- model_covariates_list$model2
strat_analysis_cov <- strat_analysis_cov[!strat_analysis_cov %in% c("chronotype")]

ins_results_stratified <- run_exposure_models_stratified(
  data = data,
  outcome = "chronic_insomnia",
  stratify_by = "chronotype_3cat",
  main_exposures = main_exposures,
  model_covariates = strat_analysis_cov,
  base_priors = base_ins_priors,
  m_treedepth = 10,
  results_dir = "./results",
  label = "ins_chron_strat"
)

ssl_results_stratified <- run_exposure_models_stratified(
  data = data,
  outcome = "avg_short_sleep",
  stratify_by = "chronotype_3cat",
  main_exposures = main_exposures,
  model_covariates = strat_analysis_cov,
  base_priors = base_ssl_priors,
  m_treedepth = 10,
  results_dir = "./results",
  label = "ssl_chron_strat"
)

# NOTE: the source project also contained an unfinished, in-progress block
# fitting mutually-adjusted models per chronotype stratum. It has been left
# out of this cleaned script because, as written, it fit every stratum's
# model on the full (unsplit) dataset instead of each chronotype subset -
# i.e. it did not actually produce three different stratified models. If a
# mutually-adjusted stratified estimate is needed, this should be re-written
# using purrr::map2(data_split, names(data_split), ~ fit_mut_adj_model(data = .x, ...))
# (or similar), fitting on .x rather than on `data`.

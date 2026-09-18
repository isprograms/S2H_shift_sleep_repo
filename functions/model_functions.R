################################################################################
###################### FIT BAYESIAN MODELS #####################################
################################################################################

library(tidyverse)
library(brms)
library(broom.mixed)
library(writexl)

################################################################################
#################### Logistic regression for binary outcomes ###################
################################################################################
#* Fit Bayesian logistic regression model with logit link. This function
#* returns the fitted model object, and takes as input: dataset, outcome
#* variable as a string, vector of exposure variable names as strings, vector
#* of covariates, and prior distributions.
fit_blog_model <- function(data, outcome, exposures, covariates,
                           priors, m_treedepth = 10){
  # Define formula
  predictors <- c(exposures, covariates)

  righthand_formula <- paste0(predictors, collapse = " + ")

  formula <- as.formula(paste(outcome, "~", righthand_formula))

  # Fit model
  model <- brm(
    formula = formula,
    data = data,
    family = bernoulli(link = "logit"),
    prior = priors,
    seed = 647,
    chains = 4, iter = 4000, warmup = 2000, cores = 4,
    backend = "cmdstanr",
    control = list(max_treedepth = m_treedepth),
    refresh = 0
  )

  return(model)
}

################################################################################
####### Fit several logistic regression models and store their results #########
################################################################################
#* Run Bayesian logistic regression models for multiple exposures and covariate
run_exposure_models <- function(data,
                                outcome,
                                main_exposures,
                                model_covariates_list,
                                base_priors,
                                exposure_prior_dist = "normal(0, 1)",
                                m_treedepth = 10,
                                results_dir,
                                label = "model_results") {

  plot_dir <- file.path(results_dir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  all_results_list <- list()

  for (exposure in main_exposures) {
    tibble_name <- paste0("models_", exposure)
    all_results_list[[tibble_name]] <- tibble()

    # 1. Determine coefficient names for the specific prior
    if (is.factor(data[[exposure]])) {
      full_formula <- as.formula(paste(outcome, "~", exposure))
      priors_info <- get_prior(full_formula, data = data, family = bernoulli())
      coef_names <- priors_info$coef[priors_info$class == "b" &
                                       priors_info$coef != "(Intercept)" &
                                       priors_info$coef != ""]
    } else {
      coef_names <- exposure
    }

    # 2. Build the specific prior using the custom distribution argument
    current_priors <- c(
      base_priors,
      set_prior(exposure_prior_dist, class = "b", coef = coef_names)
    )

    # 3. Loop over models
    for (m_idx in seq_along(model_covariates_list)) {
      cov_set <- model_covariates_list[[m_idx]]
      model_label <- names(model_covariates_list)[m_idx]
      model_name <- paste0(label, "_", exposure, "_", model_label)

      cat("\n--- Fitting:", model_name, "---\n")

      current_model <- fit_blog_model(
        data = data, outcome = outcome, exposures = c(exposure),
        covariates = cov_set, priors = current_priors,
        m_treedepth = m_treedepth
      )

      # Save objects
      saveRDS(current_model, file.path(results_dir, paste0(model_name, ".rds")))
      p <- pp_check(current_model, ndraws = 100) + ggtitle(model_name)
      ggsave(file.path(plot_dir, paste0("ppcheck_", model_name, ".png")), plot = p)

      # 4. Extract Results
      tidy_results <- tidy(current_model, exponentiate = TRUE, conf.level = 0.95) |>
        filter(term %in% coef_names) |>
        mutate(
          "OR (95% CI)" = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
          model = model_label,
          Exposure = term
        ) |>
        select(model, Exposure, "OR (95% CI)")

      all_results_list[[tibble_name]] <- bind_rows(all_results_list[[tibble_name]], tidy_results)
    }
  }

  writexl::write_xlsx(all_results_list, file.path(results_dir, paste0(label, "_summary.xlsx")))
  return(all_results_list)
}

################################################################################
##########################Mutually adjusted models #############################
################################################################################
#* Fit a model with several exposures entered simultaneously. `priors` holds
#* the intercept and covariate priors; the exposure coefficients get their own
#* prior (`exposure_prior_dist`), as in run_exposure_models().
fit_mut_adj_model <- function(
    data, outcome, exposures, covariates, priors, results_dir, model_name,
    exposure_prior_dist = "normal(0, 1)"){
  ## Get names of coefficients of main exposures
  exposures_formula <- paste0(exposures, collapse = " + ")
  full_formula <- as.formula(paste(outcome, "~", exposures_formula))
  priors_info <- get_prior(full_formula, data = data, family = bernoulli())
  coef_names <- priors_info$coef[priors_info$class == "b" &
                                   priors_info$coef != "(Intercept)" &
                                   priors_info$coef != ""]

  current_priors <- c(
    priors,
    set_prior(exposure_prior_dist, class = "b", coef = coef_names)
  )

  # Fit model
  current_model <- fit_blog_model(
    data = data,
    outcome = outcome,
    exposures = exposures,
    covariates = covariates,
    priors = current_priors
  )

  # Save objects
  saveRDS(current_model, file.path(results_dir, paste0(model_name, ".rds")))
  p <- pp_check(current_model, ndraws = 100) + ggtitle(model_name)
  ggsave(file.path(results_dir, paste0("ppcheck_", model_name, ".png")), plot = p)

  # Extract Results
  tidy_results <- tidy(current_model, exponentiate = TRUE, conf.level = 0.95) |>
    filter(term %in% coef_names) |>
    mutate(
      "OR (95% CI)" = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high),
      Exposure = term
    ) |>
    select(Exposure, "OR (95% CI)")

  writexl::write_xlsx(tidy_results, file.path(results_dir, paste0(model_name, "_summary.xlsx")))
  return(tidy_results)
}

################################################################################
##########################Stratified analyses #############################
################################################################################
#* Fit one model per exposure and stratum. Categorical exposures return one
#* row per level (prevalence of the outcome and OR, reference level without
#* OR); continuous exposures return one row with the OR per unit increase.
run_exposure_models_stratified <- function(data,
                                           outcome,
                                           main_exposures,
                                           model_covariates,
                                           base_priors,
                                           stratify_by,
                                           exposure_prior_dist = "normal(0, 1)",
                                           m_treedepth = 10,
                                           results_dir,
                                           label = "model_results") {

  # Setup directories
  plot_dir <- file.path(results_dir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  strata_levels <- unique(na.omit(data[[stratify_by]]))
  all_results_list <- list()

  for (exposure in main_exposures) {
    tibble_name <- paste0("models_", exposure)
    all_results_list[[tibble_name]] <- tibble()

    for (stratum in strata_levels) {
      # 1. Subset and Clean
      stratum_data <- data[data[[stratify_by]] == stratum, ]
      stratum_data <- droplevels(stratum_data)

      # Safety check: skip if empty
      if (nrow(stratum_data) == 0) next

      is_categorical <- is.factor(stratum_data[[exposure]]) ||
        is.character(stratum_data[[exposure]])

      # 2. Determine priors
      full_formula <- as.formula(paste(outcome, "~", exposure))
      priors_info <- get_prior(full_formula, data = stratum_data, family = bernoulli())
      coef_names <- priors_info$coef[priors_info$class == "b" &
                                       priors_info$coef != "(Intercept)" &
                                       priors_info$coef != ""]

      if (length(coef_names) == 0) next

      current_priors <- c(
        base_priors,
        set_prior(exposure_prior_dist, class = "b", coef = coef_names)
      )

      # Match each non-reference level to its coefficient name (same order)
      if (is_categorical) {
        exposure_levels <- levels(factor(stratum_data[[exposure]]))
        if (length(coef_names) != length(exposure_levels) - 1) {
          stop("Could not match levels of ", exposure, " to coefficients in stratum ", stratum)
        }
        level_map <- tibble(level = exposure_levels[-1], term = coef_names)
      }

      # 3. Fit models
      model_name <- paste0(label, "_", exposure, "_", stratum)

      cat("\n--- Fitting:", model_name, "---\n")

      # Call the model fitting function
      current_model <- fit_blog_model(
        data = stratum_data,
        outcome = outcome,
        exposures = c(exposure),
        covariates = model_covariates,
        priors = current_priors,
        m_treedepth = m_treedepth
      )

      # Save individual model object
      saveRDS(current_model, file.path(results_dir, paste0(model_name, ".rds")))

      # 4. Extract results
      or_tab <- tidy(
        current_model,
        exponentiate = TRUE,
        conf.level = 0.95
      ) |>
        filter(term %in% coef_names) |>
        mutate("OR (95% CI)" = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high)) |>
        select(term, "OR (95% CI)")

      # Only proceed if we actually have coefficients to report
      if (nrow(or_tab) == 0) next

      if (is_categorical) {
        # Prevalence of the outcome (N (%)) in each level of the exposure
        prev_map <- stratum_data |>
          mutate(level = trimws(as.character(.data[[exposure]]))) |>
          group_by(level) |>
          summarise(
            n_event = sum(.data[[outcome]], na.rm = TRUE),
            total = n(),
            formatted = sprintf("%d (%.1f%%)", n_event, 100 * n_event / total),
            .groups = "drop"
          )

        tidy_results <- prev_map |>
          left_join(level_map, by = "level") |>
          left_join(or_tab, by = "term") |>
          mutate(Stratum = as.character(stratum)) |>
          select(!!sym(exposure) := level, formatted, Stratum, "OR (95% CI)")
      } else {
        tidy_results <- or_tab |>
          mutate(
            !!sym(exposure) := exposure,
            formatted = NA_character_,
            Stratum = as.character(stratum)
          ) |>
          select(!!sym(exposure), formatted, Stratum, "OR (95% CI)")
      }

      all_results_list[[tibble_name]] <- bind_rows(all_results_list[[tibble_name]], tidy_results)
    }
  }

  # Export to Excel
  writexl::write_xlsx(all_results_list, file.path(results_dir, paste0(label, "_summary.xlsx")))

  return(all_results_list)
}

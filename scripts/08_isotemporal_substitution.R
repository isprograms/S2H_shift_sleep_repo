################################################################################
################ ISOTEMPORAL SUBSTITUTION: FIGURE 2 #############################
################################################################################

# by: Isabel Santonja
# description: This script fits the Bayesian logistic regression models for
# the isotemporal substitution analysis (reallocating day shifts to other
# shift clusters) and produces Figure 2 and Supplementary Figures 5-6, showing
# the association between shift cluster schedule composition and chronic
# insomnia / short average sleep duration.
# last updated: 18/09/2026

# This script reads "data/brm_data.RDS", built by 05_prepare_model_data.R: a
# participant-level dataset with one row per participant, containing the
# relative frequency of each shift cluster (day/evening/night/long day/long
# evening), the two sleep outcomes (chronic_insomnia, avg_short_sleep) and all
# covariates.

# SET-UP ####
# Load libraries
library(tidyverse)
library(brms)
library(broom.mixed)
library(patchwork)

source(here::here("functions", "model_functions.R"))
source(here::here("functions", "tables.R"))

# Load data
data <- readRDS(file = here::here("data", "brm_data.RDS")) |>
  mutate(
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
    ),
    years_worked_ns_grouped = fct_relevel(
      years_worked_ns_grouped,
      "0 years", "up to 5 years", "6-15 years", "16+ years"
    ),
    across(
      c(night_freq, evening_freq, long_day_freq, long_evening_freq),
      ~ case_when(
        . == 0 ~ "0",
        . > 0 & . < 0.15 ~ "10",
        . >= 0.15 & . < 0.30 ~ "20",
        . >= 0.30 & . < 0.50 ~ "40",
        . >= 0.50 & . < 0.70 ~ "60",
        . >= 0.70 & . < 0.90 ~ "80",
        . >= 0.900 ~ "100",
      ),
      .names = "{.col}_cat"
    )
  ) |>
  mutate(across(ends_with("_cat"), factor)) |>
  mutate(across(c("night_freq_cat", "evening_freq_cat",
                  "long_day_freq_cat", "long_evening_freq_cat"),
                ~ fct_relevel(., "0", "10", "20", "40", "60", "80", "100")))

# FUNCTIONS ####
## Dataframe of ORs
get_or_tibble <- function(model) {
  require(broom.mixed)

  or_tibble <- tidy(model, exponentiate = TRUE) |>
    select(term, estimate, conf.low, conf.high) |>
    filter(grepl("_cat", term)) |>
    mutate(
      credible = if_else(
        conf.high < 1.00 | conf.low > 1.0,
        TRUE, FALSE
      ),
      # Extract shift type
      shift_type = factor(str_extract(term, ".*_cat")),
      # Extract % (digits at the end of the string)
      freq = factor(str_extract(term, "\\d+$")),
      shift_type = factor(case_when(
        shift_type == "evening_freq_cat" ~ "Evening shifts",
        shift_type == "night_freq_cat" ~ "Night shifts",
        shift_type == "long_day_freq_cat" ~ "Long day shifts",
        shift_type == "long_evening_freq_cat" ~ "Long evening shifts"
      )),
      shift_type = fct_relevel(shift_type, "Evening shifts", "Night shifts")
    ) |>
    select(-term) |>
    select(shift_type, freq, everything()) |>
    mutate(freq = case_when(
      freq == "0" ~ "0%",
      freq == "10" ~ ">0 to 15%",
      freq == "20" ~ ">15 to 30%",
      freq == "40" ~ ">30 to 50%",
      freq == "60" ~ ">50 to 70%",
      freq == "80" ~ ">70 to 90%",
      freq == "100" ~ ">90%"
    )) |>
    mutate(freq = fct_relevel(
      freq, ">90%", ">70 to 90%", ">50 to 70%", ">30 to 50%", ">15 to 30%",
      ">0 to 15%", "0%"
    ))

  return(or_tibble)
}

## Function to create Substitution Forest Plots
draw_substitution_plot <- function(plot_data, title, outcome_name) {
  # Create x axis title
  x_title <- paste0(outcome_name, ": OR (95% CI)")

  # Draw plot
  ggplot(plot_data,
         aes(x = estimate, y = freq)) +
    # Add the coloured point for estimate
    geom_point(aes(color = shift_type), size = 3) +
    # Add the CI bars
    geom_errorbarh(aes(xmin = conf.low,
                       xmax = conf.high,
                       colour = shift_type),
                   height = 0.2,
                   linewidth = 1) +
    # Add reference line at OR = 1
    geom_vline(xintercept = 1,
               linetype = "dashed",
               color = "gray50") +
    # Facet by shift type with free y-axis scales
    facet_wrap(~ shift_type, ncol = 1, scales = "free_y") +
    # Color scale matching
    scale_colour_manual(values = shift_type_colours) +
    labs(
      title = title,
      x = x_title,
      y = "Relative frequency (%) of shift cluster added",
      colour = "Shift type"
      ) +
    theme_minimal() +
    theme(
      strip.text = element_text(
        size = 10,
        face = "bold",
        hjust = 0 # Align left
        ),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16),
      axis.text.y = element_text(size = 10),
      plot.caption = element_text(
        size = 10,
        hjust = 0,
        ),
      legend.position = "none"
  )
}

# PALETTES ####
shift_type_colours <- c("Day shifts" = "#A6D854",
                        "Evening shifts"   = "#FC8D62",
                        "Night shifts" = "#8DA0CB",
                        "Long day shifts" = "#66C2A5",
                        "Long evening shifts" = "#E78AC3"
)

# MODEL SPECIFICATIONS ####
## Variables
exposures <- c("evening_freq_cat", "night_freq_cat",
               "long_day_freq_cat", "long_evening_freq_cat")
## Define the sets of covariates in a named list
model_covariates_list <- list(
  "model1" = c("age_scaled", "gender_grouped"),
  "model2" = c("age_scaled", "gender_grouped", "chronotype_3cat", "country_grouped",
                "education_grouped", "income_satisfaction",
                "marital_status_grouped", "alcohol", "bmi_scaled",
                "smoking_status"),
  "model3" = c("age_scaled", "gender_grouped", "chronotype_3cat", "country_grouped",
               "education_grouped", "income_satisfaction",
               "marital_status_grouped", "alcohol", "bmi_scaled",
               "smoking_status", "years_worked_ns_grouped", "working_hours_scaled")
)
## Priors
coef_frequencies <- c("night_freq_cat10", "night_freq_cat20", "night_freq_cat40",
                      "night_freq_cat60", "night_freq_cat80", "night_freq_cat100",
                      "evening_freq_cat10", "evening_freq_cat20", "evening_freq_cat40",
                      "evening_freq_cat60", "evening_freq_cat80", "evening_freq_cat100",
                      "long_day_freq_cat10", "long_day_freq_cat20", "long_day_freq_cat40",
                      "long_day_freq_cat60", "long_day_freq_cat80", "long_day_freq_cat100",
                      "long_evening_freq_cat10", "long_evening_freq_cat20",
                      "long_evening_freq_cat40", "long_evening_freq_cat60",
                      "long_evening_freq_cat80", "long_evening_freq_cat100")

ins_priors <- c(
  set_prior("normal(-2.2, 1)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("normal(0, 1.0)", class = "b",
            coef = coef_frequencies
            ))

ssl_priors <- c(
  set_prior("normal(-2.7, 1)", class = "Intercept"),
  set_prior("normal(0, 0.5)", class = "b"),
  set_prior("normal(0, 1.0)", class = "b",
            coef = coef_frequencies)
)

# PLOTS ####
## Figure 2A: risk of insomnia (Model 2) ####
# Use permanent day shift as reference category
ins_m2 <- fit_blog_model(
  data = data,
  outcome = "chronic_insomnia",
  exposures = exposures,
  covariates = model_covariates_list[["model2"]],
  priors = ins_priors
)

ins_m2_tb <- get_or_tibble(ins_m2)

fig2a <- draw_substitution_plot(ins_m2_tb, "A)", "Chronic insomnia") +
  labs(
    subtitle = "Reference: permanent day shift workers (100% day shifts)"
  )

## Figure 2B: risk of short average sleep duration (Model 2) ####
ssl_m2 <- fit_blog_model(
  data = data,
  outcome = "avg_short_sleep",
  exposures = exposures,
  covariates = model_covariates_list[["model2"]],
  priors = ssl_priors
)

ssl_m2_tb <- get_or_tibble(ssl_m2)

fig2b <- draw_substitution_plot(ssl_m2_tb, "B)", "Short average sleep duration") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    subtitle = element_blank()
  )

## Figure 2 ####
fig2 <- fig2a + fig2b

ggsave(here::here("figures", "fig2.pdf"), width = 8, height = 10)

## Supplementary Figure 5 (Model 1) ####
ins_m1 <- fit_blog_model(
  data = data,
  outcome = "chronic_insomnia",
  exposures = exposures,
  covariates = model_covariates_list[["model1"]],
  priors = ins_priors
)

ins_m1_tb <- get_or_tibble(ins_m1)

supfig5a <- draw_substitution_plot(ins_m1_tb, "A)", "Chronic insomnia") +
  labs(
    subtitle = "Reference: permanent day shift workers (100% day shifts)"
  )

ssl_m1 <- fit_blog_model(
  data = data,
  outcome = "avg_short_sleep",
  exposures = exposures,
  covariates = model_covariates_list[["model1"]],
  priors = ssl_priors
)

ssl_m1_tb <- get_or_tibble(ssl_m1)

supfig5b <- draw_substitution_plot(ssl_m1_tb, "B)", "Short average sleep duration") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )

supfig5 <- supfig5a + supfig5b
ggsave(here::here("figures", "supfig5.pdf"), width = 8, height = 10)

## Supplementary Figure 6 (Model 3, mutually adjusted) ####
ins_m3 <- fit_blog_model(
  data = data,
  outcome = "chronic_insomnia",
  exposures = exposures,
  covariates = model_covariates_list[["model3"]],
  priors = ins_priors
)

ins_m3_tb <- get_or_tibble(ins_m3)

supfig6a <- draw_substitution_plot(ins_m3_tb, "A)", "Chronic insomnia") +
  labs(
    subtitle = "Reference: permanent day shift workers (100% day shifts)"
  )

ssl_m3 <- fit_blog_model(
  data = data,
  outcome = "avg_short_sleep",
  exposures = exposures,
  covariates = model_covariates_list[["model3"]],
  priors = ssl_priors
)

ssl_m3_tb <- get_or_tibble(ssl_m3)

supfig6b <- draw_substitution_plot(ssl_m3_tb, "B)", "Short average sleep duration") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.caption = element_blank()
  )

supfig6 <- supfig6a + supfig6b
ggsave(here::here("figures", "supfig6.pdf"), width = 8, height = 10)

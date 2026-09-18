################################################################################
##################### PREPARE PARTICIPANT-LEVEL MODEL DATASET ##################
################################################################################

# by: Isabel Santonja
# description: This script builds "data/brm_data.RDS", the participant-level
# dataset (one row per participant) used by the sample description
# (Sup. Table S2) and the regression scripts (Figure 2, Sup. Figures 5-6,
# Tables 2-3, Sup. Tables S3-S8). It combines the shift-level sleep data with
# participant-level covariates and shift-cluster composition, and derives
# the outcome variables (average sleep duration, short average sleep
# duration, continuous/chronic insomnia) and the scaled/regrouped versions
# of the covariates used in the regression models.
#Last updated: 18/09/2026

# This script uses two datasets produced earlier in the pipeline:
# "data/data_analysis" (one row per participant, from
# 02_initial_data_analysis_and_wrangling.Rmd) and "data/schedules_sleep_long.RDS"
# (one row per shift, from 04_shift_specific_sleep.R). It first builds the
# participant-level dataset "data/level2_data.RDS" (participant characteristics
# plus the relative frequency of time spent in each shift cluster) and then
# derives the model dataset from it.

library(tidyverse)

# Load data
wide_data <- readRDS(file = here::here("data", "data_analysis"))
long_data <- readRDS(file = here::here("data", "schedules_sleep_long.RDS"))

# 1. Harmonise shift-cluster names and levels ####
long_data <- long_data |>
  mutate(
    shift_type_4 = fct_recode(shift_type_4,
                              Early = "Very Early",
                              Day = "Early",
                              Evening = "Late"),
    cluster_k5 = factor(cluster_k5),
    cluster_k5 = factor(case_when(
      cluster_k5 == "Afternoon shifts" ~ "Evening shifts",
      cluster_k5 == "Long early shifts" ~ "Long day shifts",
      cluster_k5 == "Long late shifts" ~ "Long evening shifts",
      .default = cluster_k5
      ))
  ) |>
  mutate(
    shift_type_4 = fct_relevel(shift_type_4,
                               "Early", "Day"),
    cluster_k5 = fct_relevel(cluster_k5,
      "Day shifts", "Evening shifts", "Night shifts"
    )
  )

# 2. Participant-level dataset  ####
## Years worked on night shifts (0 for always day workers) and number of
## distinct shift schedules
l2_data <- wide_data |>
  mutate(
    years_worked_ns = as.numeric(case_when(
      is.na(years_worked_ns) & night_shift_history == "Always day worker" ~ 0,
      .default = years_worked_ns
    )),
    years_worked_ns_log = log1p(years_worked_ns),
    number_of_shifts_cat = factor(number_of_shifts)
  )

## Relative frequency (0-1) of time spent in each shift cluster, per participant
shift_freq_wide <- long_data |>
  select(id, cluster_k5, shift_freq_rel_clean) |>
  group_by(id, cluster_k5) |>
  summarise(
    shift_cluster_rel_freq = sum(shift_freq_rel_clean),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = cluster_k5,
    values_from = shift_cluster_rel_freq,
    values_fill = 0 # clusters the participant never worked
  )

l2_data <- l2_data |>
  left_join(shift_freq_wide, by = "id") |>
  rename(
    day_freq = `Day shifts`,
    evening_freq = `Evening shifts`,
    night_freq = `Night shifts`,
    long_day_freq = `Long day shifts`,
    long_evening_freq = `Long evening shifts`
  )

saveRDS(l2_data, file = here::here("data", "level2_data.RDS"))

rm(wide_data, shift_freq_wide)

# 3. Average sleep duration and short average sleep duration ####
## Average sleep duration on workdays (frequency-weighted across shifts)
sleep_summary_wd <- long_data |>
  group_by(id) |>
  summarise(
    avg_workdays_sleep_duration = sum(sleep_duration * shift_freq_rel_clean)
  ) |>
  mutate(
    workdays_short_sleep = if_else(
      avg_workdays_sleep_duration < 6,
      TRUE,
      FALSE
    )
  )

## Add sleep duration on days off, and compute the overall weighted average
## (5 workdays + 2 days off per week; see manuscript Methods)
sleep_summary <- l2_data |>
  select(id, oq_wf_dw_sleep_duration, oq_wf_sleep_duration) |>
  mutate(
    freedays_sleep_duration = as.numeric(if_else(
      is.na(oq_wf_dw_sleep_duration) & !is.na(oq_wf_sleep_duration),
      oq_wf_sleep_duration,
      oq_wf_dw_sleep_duration
    ))
  ) |>
  select(id, freedays_sleep_duration) |>
  left_join(sleep_summary_wd, by = "id") |>
  mutate(
    avg_sleep_duration =
      ((avg_workdays_sleep_duration * 5) + (freedays_sleep_duration * 2)) / 7,
    avg_short_sleep = if_else(
      avg_sleep_duration < 6,
      TRUE,
      FALSE
    )
  )

# Add sleep duration vars to analysis dataset
brm_data <- l2_data |>
  inner_join(sleep_summary, by = "id")

rm(l2_data, sleep_summary_wd, sleep_summary)

# 4. Continuous insomnia score and working hours scaling ####
brm_data <- brm_data |>
  mutate(
    # replace missings in sleep interference with 0
    # (question with branching logic)
    sp_interference = if_else(
      is.na(oq_sp_interference),
      0,
      oq_sp_interference
    ),
    isi_red = oq_sp_a + oq_sp_b + oq_sp_c + sp_interference
    ) |>
  mutate(
    working_hours_scaled = as.numeric( # standardize variable
      scale(working_hours)
      )
    )

# 5. Scale continuous covariates, regroup small categories, set reference levels ####
brm_data <- brm_data |>
  mutate(
    age_scaled = as.numeric(scale(oq_age)),
    bmi_scaled = scale(bmi_imputed),
    years_worked_ns_log_scaled = as.numeric(scale(years_worked_ns_log))
  ) |>
  mutate( # regroup small categories
    gender_grouped = if_else(
      gender %in% c(
        "Non-binary",
        "Prefer not to answer",
        "Other"),
      "Other",
      gender
    ),
    country_grouped = fct_lump_prop(
      country, prop = 0.05, other_level = "Other"
      ),
    country_grouped = if_else(
      is.na(country_grouped),
      "Unknown",
      country_grouped
    ),
    work_sector_grouped = fct_lump_prop(
      work_sector, prop = 0.05, other_level = "Other"
      ),
    chronotype_3cat = case_when(
      chronotype %in% c("Definitely morning", "Rather morning") ~ "Early",
      chronotype %in% c("Definitely evening", "Rather evening") ~ "Late",
      .default = "Intermediate"
    ),
     marital_status_grouped = fct_collapse(
      marital_status,
      "Divorced, separated, or partner has died" = c(
        "Divorced/separated", "Widowed/partner has died"
        )
    ),
    marital_status_grouped = fct_relevel(
      marital_status_grouped,
      "Single", "Married/living together with my partner",
      "Divorced, separated, or partner has died"
      ),
    alcohol = fct_collapse(
      alcohol_weekly,
      "None" = "0 drink",
      "2 - 3 drinks" = c(
        "2 drinks", "3 drinks"
          ),
      "10 or more drinks" = c(
          "10 - 15 drinks", "16 or more drinks"
          )
      ),
    alcohol = fct_relevel(alcohol, "None", "1 drink",
                          "2 - 3 drinks",
                          "4 - 5 drinks",
                          "6 - 9 drinks",
                          "10 or more drinks"
                          )
  ) |>
  mutate( # specify new reference categories
    shift_rotation_speed = fct_relevel(
      shift_rotation_speed,
      "Not applicable", "Daily change", "Change every 2-4 days",
      "Change every week", "Change every 2-3 weeks",
      "Rotation speed varies", "Other"
      ),
    night_shifts_updated = fct_relevel(
      night_shifts_updated,
      "Always day worker", "Former night shift worker"
      )
  ) |>
  mutate( # make factors unordered
    income_satisfaction = factor(
      income_satisfaction,
      ordered = FALSE
    ),
    smoking_status = factor(
      smoking_status,
      ordered = FALSE
    )
  )

# 6. Save data ####
saveRDS(brm_data, file = here::here("data", "brm_data.RDS"))

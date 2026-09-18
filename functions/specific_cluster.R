################################################################################
################## SPECIFIC CLUSTER FUNCTIONS ##################################
################################################################################
#* by: Isabel Santonja
#* created on: 06/08/2025
#* last modified: 
#* description: "Functions to use with the clusters of the online survey 
#* dataset."
#* WARNING: "Many of these functions are very specific to the dataset and
#* probably cannot be used in other contexts without modification."

library(tidyverse)
library(circular)
library(hms)

################################################################################
############ Function to compute cluster-wise summary stats ###################
################################################################################
#* This function takes two arguments: a tibble with the variables in circular 
#* format and the factor variable that defines the clusters. It returns a
#* tibble with the summary statistics for each cluster.
cluster_schedules_summary <- function(dataset, cluster_var){
  dataset |>  
    group_by({{cluster_var}}) |>
    summarise(
      n_shifts = n(),
      mean_start = mean.circular(shift_start_circular, na.rm = TRUE),
      sd_start = sd.circular(shift_start_circular, na.rm = TRUE),
    
      mean_end = mean.circular(shift_end_circular, na.rm = TRUE),
      sd_end = sd.circular(shift_end_circular, na.rm = TRUE),
    
      mean_duration = mean(shift_duration_clean, na.rm = TRUE),
      sd_duration = sd(shift_duration_clean, na.rm = TRUE),
    
      night_shift_n = sum(is_night_shift_num, na.rm = TRUE),
      night_shift_pct = mean(is_night_shift_num, na.rm = TRUE) * 100
      ) |>
    mutate(
      mean_start = ifelse(mean_start < 0, mean_start + 2 * pi, mean_start),
      mean_end = ifelse(mean_end < 0, mean_end + 2 * pi, mean_end)
      )|>
    mutate(
      mean_start_hour = (as.numeric(mean_start) * 24) / (2 * pi),
      sd_start_hour = (as.numeric(sd_start) * 24) / (2 * pi),
      mean_end_hour = (as.numeric(mean_end) * 24) / (2 * pi),
      sd_end_hour = (as.numeric(sd_end) * 24) / (2 * pi)
      ) |> 
    select(-mean_start, -mean_end, -sd_start, -sd_end)|>
    select(
      {{cluster_var}}, n_shifts, 
      mean_start_hour, sd_start_hour, 
      mean_end_hour, sd_end_hour, 
      mean_duration, sd_duration, 
      night_shift_n,
      night_shift_pct
      )|>
    mutate(
      across(mean_start_hour:night_shift_pct,
             ~ round(.x, 2))
      )
}
################################################################################
########################## SCHEDULES CLEANING  #################################
################################################################################

# by: Isabel Santonja
# last modified: 18/09/2026
# description: This script contains minimal corrections to the the long shifts 
# dataset.

# 1. Load data and libraries
library(dplyr)
library(hms)
library(magrittr)

shifts_long <- readRDS(here::here("raw_data", "shift_schedules_long.RDS"))

# 2. Fixing mistakes
## 2.1 Correct one respondent's reported shift end time
## One respondent reported shift start/end times of 07:00-07:00 (implying a
## 24h shift), but had separately reported/corrected the shift duration to
## 17h. An earlier cleaning pass had recalculated the end time to 00:00 to
## match that shorter duration; this is reverted here, keeping the originally
## reported end time (07:00), since the shorter duration likely reflects the
## respondent excluding on-shift sleep time from their estimate, rather than
## an error in the reported end time.

shifts_long_edited <- shifts_long %>%
  mutate(
    shift_end_clean = as_hms(ifelse(
      id == 103307 & shift_end_clean == as_hms("00:00:00"), 
      as_hms("07:00:00"),
      shift_end_clean)
  ))

# 3. Save the cleaned data

saveRDS(shifts_long_edited, here::here("data", "shift_schedules_long_edited.rds"))

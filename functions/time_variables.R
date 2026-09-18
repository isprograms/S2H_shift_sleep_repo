################################################################################
################## FUNCTIONS FOR TIME VARIABLES ################################
################################################################################

#* by: Isabel Santonja
#* created on: 07/07/2025
#* last modified: 
#* description: "Functions to work with time variables"


################################################################################
########################### Function to round times ############################
################################################################################
#* This function takes an hms variable and a numeric value representing the 
#* degree of rounding (in h, e.g. 1.0h or 0.25h).
#* It defaults to rounding to the nearest half hour
round_hms_to_nearest <- function(time, interval_hours = 0.5) {
  interval_seconds = 3600 * interval_hours
  as_hms(round(as.numeric(time) / interval_seconds) * interval_seconds)
}
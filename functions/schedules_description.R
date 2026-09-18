################################################################################
################ SCHEDULES DESCRIPTION FUNCTIONS ###############################
################################################################################

#* by: Isabel Santonja
#* created on: 04/07/2025
#* last modified: 16/07/2026
#* description: "Functions to describe the (long) schedules dataset"
#* WARNING: "Many of these functions are very specific to the dataset and
#* should not be used in other contexts without modification."

################################################################################
## Function to calculate summary statistics of shift duration per shift type ###
################################################################################
#* This function takes the binary variable indicating the type of shift as 
#* argument and returns a table with summary statistics of the shift duration
shift_dur_by_type <- function(type_of_shift, schedules){
  # Build string with shift type
  ## Capture the expression
  type_of_shift_quo <- enquo(type_of_shift)
  ## Get the variable name as a string
  type_name <- as_label(type_of_shift_quo)
  ## Remove "is_" prefix
  shift_type_clean <- str_remove(type_name, "^is_")
  shift_type_clean <- str_remove(shift_type_clean, "_shift")
  # Build table
  schedules %>% 
    filter({{type_of_shift}}) %>% 
    ungroup() %>% 
    summarise(
      n = sum(!is.na(shift_duration_clean)),
      mean = mean(shift_duration_clean, na.rm = TRUE),
      sd = sd(shift_duration_clean, na.rm = TRUE),
      median = median(shift_duration_clean, na.rm = TRUE),
      p25 = quantile(shift_duration_clean, probs = 0.25, na.rm = TRUE),
      p75 = quantile(shift_duration_clean, probs = 0.75, na.rm = TRUE)
    ) %>% 
    mutate("Shift type" = shift_type_clean) %>% 
    select("Shift type", everything())
}

################################################################################
######## Function to generate table with summary statistics of shift ###########
################# duration variable stratified by shift type ###################
################################################################################
#* This function takes a list of shift types variablesand returns a table with 
#* summary and one row per specified shift type
table_duration_by_type <- function(var_list, schedules) {
  # Create empty list to store results
  results_list <- list()
  # Ensure the variable list is a character vector
  if (!is.character(var_list)) {
    stop("var_list must be a character vector of variable names.")
  }
  # Loop through each variable in the list
  for (var in var_list) {
    # Call the shift_dur_by_type function for each variable
    result <- shift_dur_by_type(!!sym(var), schedules)
    results_list[[var]] <- result |> 
      mutate(across(where(is.numeric), ~ round(.x, digits = 2)))
  }
  
  # Combine all results into a single data frame
  combined_results <- bind_rows(results_list)

  return(combined_results)
}

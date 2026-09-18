################################################################################
#########################CUSTOM TABLES##########################################
################################################################################
library(tidyverse)
library(janitor)


################################################################################
############################PREVALENCE TABLE####################################
################################################################################
prev_tab <- function(data, col_var, row_var){
  # Cross-tab with counts 
  tab_counts <- tabyl(data, {{row_var}}, {{col_var}})|>
    adorn_totals("row")
  
  # Cross-tab with row-wise %
  tab_percent <- tab_counts %>%
    adorn_percentages("row") %>%
    adorn_rounding(digits = 3)
  
  # Combine counts and percentages in one table
  tab_combined <- tab_counts
  
  for (col in names(tab_counts)[-1]) {  # skip first column (row labels)
    tab_combined[[col]] <- paste0(tab_counts[[col]], 
                                  " (", 
                                  tab_percent[[col]]*100, 
                                  "%)")
  }
  
  tab_combined
}

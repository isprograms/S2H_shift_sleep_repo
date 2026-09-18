################################################################################
#########################PLOTS##################################################
################################################################################
library(tidyverse)
library(hms)
library(scales)

source(here::here("functions", "time_variables.R"))

################################################################################
###########################CREATE HISTOGRAMS####################################
################################################################################
histogram <- function(data, column, bins = 30) {
  ggplot(data, aes(x = {{ column }})) +
    geom_histogram(aes(y = after_stat(count / sum(count)) * 100),
                   colour = "grey", fill = "lightblue", bins = bins) +
    ylab("Proportion (%)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}


################################################################################
###################CREATE WRAPPED HISTOGRAMS####################################
################################################################################
wrapped_histogram <- function(data, column, wrap, bins = 30, rows = 1) {
  histogram(data, {{ column }}, bins) +
    facet_wrap(vars({{ wrap }}), nrow = rows)
}

################################################################################
##########################DUMBBELL PLOT FOR SHIFT SCHEDULES#####################
################################################################################
#* This function creates a dumbbell-style plot for shift schedules, where each 
#* line represents a distinct schedule with its start and end times. 
#* It takes 2 arguments: a tibble with variables for shift start 
#* (shift_start_clean)and end times (shift_end_clean), 
#* and a variable to colour the lines by.
schedules_plot_coloured <- function(dataset, colour_var){
  # Ensure the dataset has the necessary columns
  plot_dataset <- dataset |>
    # Fix 24h format for plotting
    mutate(
      shift_start_clean = if_else(
        shift_start_clean == as_hms("24:00:00"),
        as_hms("00:00:00"),
        shift_start_clean
      ),
      shift_end_clean = as_hms(
        if_else(
          shift_start_clean >= shift_end_clean |
            shift_duration_clean > 24,
          as.numeric(shift_end_clean) + 86400,  # Add 24 hours in seconds
          as.numeric(shift_end_clean)
        )
      )
    ) |> 
    mutate(
      shift_start_rounded = round_hms_to_nearest(
        shift_start_clean,
        interval_hours = 1.0
      ),
      shift_end_rounded = round_hms_to_nearest(
        shift_end_clean,
        interval_hours = 1.0
      )
    ) |>
    mutate(
      shift_start_rounded = as_hms(
        if_else(
          shift_start_rounded >= as_hms("24:00:00") &
            shift_end_rounded > as_hms("24:00:00"),
          as.numeric(shift_start_rounded) - 86400, # Substract 24 hours in seconds
          as.numeric(shift_start_rounded)
        )
      ),
      shift_end_rounded = as_hms(
        if_else(
          (shift_start_rounded == as_hms("00:00:00") |
             shift_start_rounded >= as_hms("24:00:00")) &
            shift_end_rounded > as_hms("24:00:00"),
          as.numeric(shift_end_rounded) - 86400,  # Substract 24 hours in seconds
          as.numeric(shift_end_rounded)
        )
      )
    ) |>
    arrange(shift_start_rounded, shift_end_rounded) |> 
    ungroup() |> 
    mutate(y_order = row_number())
  
  # Plot
  ggplot(
    plot_dataset,
    aes(
      y = y_order
    )
  ) +
    geom_segment(
      aes(
        x = shift_start_rounded,
        xend = shift_end_rounded,
        y = y_order,
        yend = y_order,
        color = {{colour_var}}
      ),
      linewidth = 0.5
    ) +
    scale_x_time(
      breaks = breaks_width("4 hours"),
      labels = time_format("%H:%M")
    ) +
    labs(
      x = "Time of Day",
      y = ""
    ) +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size=10),
      axis.ticks.y = element_blank(),
      axis.text.y =  element_blank(),
      panel.grid.major.x = element_line(
        colour = "darkgray",
        linewidth = 0.3
      ),
      panel.grid.minor.x = element_line(
        colour = "gray",
        linewidth = 0.2
      ),
      panel.background = element_blank()
    )+
    guides(
      color = guide_legend(
        override.aes = list(linewidth = 2)  # Thicker lines only in legend
      )
    )
}

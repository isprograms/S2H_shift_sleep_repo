################################################################################
###################SLEEP DURATION AND QUALITY BY SHIFT TYPE#####################
################################################################################

# by: Isabel Santonja
# Description: This script contains the models and code to create Figure 1 as  
# well as associated supplements, which show the relationship between shift type 
# sleep duration and quality, in addition to shift features. The models are 
# and linear mixed effects models, to account for repeated measures within 
# individuals.
# Last updated: 18/09/2026

# Load necessary libraries and functions
library(tidyverse)
library(hms)
library(brms)
library(tidybayes)
library(emmeans)
library(bayestestR)
library(patchwork)
library(gridExtra)

# Build the shift-specific sleep dataset ####
# Combines the shift-level cluster assignments 
# with the sleep duration/quality reported for each shift (1-4), pivoted to
# one row per shift per participant.
schedules <- readRDS(here::here("data", "schedules_cluster.RDS"))
dataset <- readRDS(here::here("data", "data_analysis"))

## Extract one sleep duration/quality value per reported shift (1-4)
make_sleep_vars <- function(data, n) {
  if (n == 1) {
    duration_vars <- c("oq_wd_sleep_duration",
                       "oq_es_sleep_duration",
                       "oq_ls_sleep_duration",
                       "oq_ns_sleep_duration")

    quality_vars  <- c("oq_wd_sleep_qual",
                       "oq_es_sleep_qual",
                       "oq_ls_sleep_qual",
                       "oq_ns_sleep_qual")
  } else {
    duration_vars <- paste0(c("oq_es_sleep_duration_",
                              "oq_ls_sleep_duration_",
                              "oq_ns_sleep_duration_"), n)

    quality_vars  <- paste0(c("oq_es_sleep_qual_",
                              "oq_ls_sleep_qual_",
                              "oq_ns_sleep_qual_"), n)
  }

  out_duration <- paste0("sleep_duration_", n)
  out_quality  <- paste0("sleep_quality_", n)

  data %>%
    mutate(
      !!out_duration := pmax(!!!syms(duration_vars), na.rm = TRUE),
      !!out_quality  := pmax(!!!syms(quality_vars),  na.rm = TRUE)
    )
}

## Pivot sleep variables to long format (one row per shift per participant)
long_sleep <- dataset |>
  select(id, matches(c("sleep_duration", "sleep_qual")), chronotype) %>%
  make_sleep_vars(1) |>
  make_sleep_vars(2) |>
  make_sleep_vars(3) |>
  make_sleep_vars(4) |>
  select(-matches("oq_")) |>
  pivot_longer(
    cols = matches("sleep_duration_|sleep_quality_"),
    names_to = c(".value", "shift_number"),
    names_pattern = "(sleep_duration|sleep_quality)_(\\d+)"
  ) |>
  mutate(shift_number = as.integer(shift_number),
         sleep_quality = factor(
           sleep_quality,
           levels = 1:5,
           labels = c("very poor", "poor", "moderate", "good", "very good"),
           ordered = TRUE
         ),
         chronotype_3cat = factor(
           case_when(
             chronotype == "Definitely morning" | chronotype == "Rather morning" ~ "Early",
             chronotype == "Definitely evening" | chronotype == "Rather evening" ~ "Late",
             .default = chronotype,
           )
         )
  ) |>
  drop_na(sleep_duration, sleep_quality)

## Merge with shift-cluster assignments
schedules_sleep <- schedules |>
  left_join(long_sleep, by = join_by(id, shift_number))

rm(dataset, long_sleep, schedules)

## Add a "very early" shift category, needed for the prespecified
## early/day/evening/night classification (distinct from the data-driven
## VMF clusters)
schedules_sleep <- schedules_sleep |>
  mutate(is_vearly_shift = if_else(
    shift_start_clean < as_hms("06:00:00") &
      (!is_night_shift | is.na(is_night_shift)),
    TRUE,
    FALSE
  )) |>
  mutate(shift_type_4 = as.factor(
    ifelse(
      is_vearly_shift,
      "Very Early",
      shift_type
    )
  ))

saveRDS(schedules_sleep, file = here::here("data", "schedules_sleep_long.RDS"))

# Load the data
data <- schedules_sleep
rm(schedules_sleep)


# SETTINGS ####
# Specify colours for plots
shift_type_colours <- c("Day shifts" = "#A6D854", 
                        "Evening shifts"   = "#FC8D62",
                        "Night shifts" = "#8DA0CB",
                        "Long day shifts" = "#66C2A5",
                        "Long evening shifts" = "#E78AC3"
                        )

sleep_qual_colours <- c("Very poor"= "#7B3294",
                        "Poor" = "#C2A5CF",
                        "Moderate" = "#F7F7F7",
                        "Good" = "#A6DBA0",
                        "Very good" = "#008837")

chronotype_colours <- c("Early" = "#d95f02",  
                        "Intermediate" = "#1b9e77",
                        "Late" = "#7570b3")

# Define custom theme for the tables
table1a_theme <- ttheme_default(
  colhead = list(
    fg_params = list(
      parse = FALSE, 
      lineheight = 1.1, # Adjusts space between the two lines
      fontface = "bold"
    ),
  base_size = 8, 
  core = list(
    fg_params = list(fontface = "plain", fontsize = 8) # Use fontsize, not size
  ),
  padding = unit(c(4, 8), "mm") 
  )
)

## Helper function to find cells in the table layout
find_cell <- function(
    table, row, col, name="core-bg"
    ){
  l <- table$layout
  which(l$t %in% row & l$l %in% col & l$name %in% name)
}

# Data wrangle
data <- data |>
  # Rename and reorder cluster levels
  mutate(
    cluster_k5 = case_when(
      cluster_k5 == "Afternoon shifts" ~ "Evening shifts",
      cluster_k5 == "Long early shifts" ~ "Long day shifts",
      cluster_k5 == "Long late shifts" ~ "Long evening shifts",
      .default = cluster_k5
    ),
    cluster_k5 = fct_relevel(
      cluster_k5,
      "Day shifts", "Evening shifts", "Night shifts"
    )
  )


# FIGURE 1A: Sleep duration by shift type ####
## Fit simple model
model_sleep_duration <- brm(
  sleep_duration ~ cluster_k5 + (1 | id),
  data = data,
  family = student(), # for coherence with main models + to model extreme values
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd"),
    set_prior("student_t(3, 0, 5)", class = "sigma")
  ),
  seed = 647,
  chains = 4, iter = 4000, warmup = 2000, cores = 4,
  backend = "cmdstanr",
)

summary(model_sleep_duration)
# Sleep duration different across shifts: evening> long day = long evening > 
# day > night

pp_check(model_sleep_duration)
# Looks fine

## Visualise shift type differences in sleep duration
### 1. Create a tiny reference grid (one row per shift)
reference_grid <- data |>
  distinct(cluster_k5) 

### 2. Get the predictive draws (expected total spread)
model_preds <- reference_grid |>
  add_predicted_draws(model_sleep_duration, re_formula = NA)

### 3. Plot
fig1a <- ggplot(
  data, 
  aes(x = cluster_k5, y = sleep_duration, fill = cluster_k5)
  ) +
  # Observed data violin plot
  geom_violin(width = 1.0, trim = FALSE) +
  # Model prediction represented as a boxplot-style interval
  stat_pointinterval(data = model_preds, 
                     aes(y = .prediction), 
                     .width = c(.5, .95), # 50% (box) and 95% (whiskers) intervals
                     color = "black", 
                     alpha = 0.7) +
  scale_y_continuous(limits = c(0, 12), n.breaks = 7) +
  scale_fill_manual(values = shift_type_colours) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(
    title = "A)", 
    x = "", 
    y = "Sleep duration (hours)",
    fill = "Shift cluster")

### 4. Associated table of differences
contrasts <- emmeans(model_sleep_duration, pairwise ~ cluster_k5)

#### Extract means
means <- contrasts$emmeans|>
  as_tibble()|>
  mutate(across(where(is.numeric), ~ round(.x, 2)))|>
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))|>
  mutate(
    "Mean (95% CI)" = paste0(
      emmean, 
      " [", 
      lower.HPD, 
      " to ", 
      upper.HPD, 
      "]"
    )
  )|>
  select(cluster_k5, "Mean (95% CI)")|>
  rename("Shift cluster" = cluster_k5)

#### Extract pds
pds <- describe_posterior(
  contrasts$contrasts,
  centrality = "mean", 
  ci = 0.95,
  test = c("p_direction")
  )|>
  as_tibble()|>
  select(contrast, pd)|>
  filter(grepl("Day shift", contrast))|>
  mutate(across(where(is.numeric), ~ round(.x, 2)))|>
  mutate(
    pd = case_when(
      pd == 1.00 ~ "> 0.99", 
      .default = as.character(pd)
    )
  )|>
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))|>
  mutate("Shift cluster" = gsub("Day shifts - ", "", contrast))|>
  select("Shift cluster", pd)


#### Assemble the table
table_fig1a_tb <- left_join(means, pds, by = "Shift cluster")|>
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))|>
  mutate("p-direction" = if_else(
    is.na(pd), 
    "REF", 
    pd
  ))|>
  select(-pd)|>
  rename(
    # Adding \n forces the text to a second line
    "Mean sleep duration\n(95% CI)" = "Mean (95% CI)"
  )

#### Convert tibble to a graphical table
table_fig1a <- tableGrob(table_fig1a_tb, rows = NULL, theme = table1a_theme)

#### Apply shift colours to table
for(i in 1:nrow(table_fig1a_tb)) {
  # Get the shift name 
  shift_name <- table_fig1a_tb$"Shift cluster"[i]
  # Get the matching color
  row_color <- shift_type_colours[shift_name]
  # Apply to the background
  # i+1 because row 1 is the header in graphical table
  cells <- find_cell(table_fig1a, row = i + 1, col = 1:ncol(table_fig1a_tb))
  for(cell_idx in cells){
    table_fig1a$grobs[[cell_idx]]$gp <- grid::gpar(fill = row_color, alpha = 1.0, col = "white")
  }
}

rm(contrasts, means, model_preds, pds, reference_grid, table_fig1a_tb)

# FIGURE 1B: Sleep quality by shift type ####
## Fit simple ordinal model
fit_quality <- brm(
  sleep_quality ~ cluster_k5 + (1 | id),
  data = data,
  family = cumulative("logit")
)
summary(fit_quality)
p_direction(fit_quality)
## pd evening and night 100%, pd long day 99.9% compared to day shifts. No credible difference between long day and day shifts.

## Get "Predictions" (probabilities) for each category
### Create the grid of clusters
grid <- data.frame(cluster_k5 = levels(data$cluster_k5))

### Get the expected probabilities (3D array)
probs_array <- posterior_epred(fit_quality, newdata = grid, re_formula = NA)

### Convert that array into a tidy data frame
fig1b_plot_data <- as.data.frame.table(
  apply(
    probs_array, 
    c(2, 3), 
    mean # Average across posterior draws to get 'Estimate'
    )) |>
  rename(
    cluster_idx = Var1, 
    quality_idx = Var2, 
    prop = Freq
    ) |>
  mutate(
    # Map back the cluster names
    cluster_k5 = levels(data$cluster_k5)[as.numeric(cluster_idx)],
    cluster_k5 = fct_relevel(
      factor(str_replace(cluster_k5, " shifts", "")),
      "Day", "Evening", "Night"),
    # Map back the quality levels (ensure these match your data labels exactly)
    sleep_quality = levels(data$sleep_quality)[as.numeric(quality_idx)],
    sleep_quality = fct_recode(
      sleep_quality,
      "Very poor" = "very poor",
      "Poor" = "poor",
      "Moderate" = "moderate",
      "Good" = "good",
      "Very good" = "very good"
      ),
    sleep_quality = fct_relevel(
      sleep_quality,
      "Very poor", "Poor", "Very good", "Good", "Moderate"
    ),
    prop = prop * 100, 
    prop_signed = ifelse(sleep_quality %in% c("Very poor", "Poor"), -prop, prop)
  )

## Plot
fig1b <- ggplot(
  fig1b_plot_data, 
  aes(x = cluster_k5, y = prop_signed, fill = sleep_quality)
  ) +
  geom_col(width = 0.7, ,
           colour = "darkgrey") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(
    name = "Sleep quality",
    values = sleep_qual_colours,
    breaks = c("Very good", "Good", "Moderate", "Poor", "Very poor")
    )+
  theme_minimal() +
  labs(title = "B)", 
       x = "Shift clusters", 
       y = "Adjusted percentage of shifts", 
       fill = "Sleep quality")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5))+ # wrap labels
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),       # Increase text size
    legend.title = element_text(size = 14, face = "bold"), # Increase title size
    legend.key.size = unit(1, "cm")#,           # Increase the size of the color boxes
    #legend.spacing.y = unit(0.5, "cm")           # Add space between legend items
  )

## Add annotations for credible differences
## Based on model estimates, if models are recalculated, these may 
## need to be updated. The stars indicate the direction and strength of the
## difference compared to day shifts.
annotations <- data.frame(
  cluster_k5 = c("Day", "Evening", "Night", "Long day"),
  label = c("Ref.", "*", "**", "*"),
  # Position them at the far end of the bars (e.g., at 105% or -105%)
  x_pos = c("Day", "Evening", "Night", "Long day"),
  y_pos = c(105, 105, 105, 105) 
)

# 2. Add to your existing plot
fig1b <- fig1b +
  geom_text(data = annotations, 
            aes(x = cluster_k5, y = y_pos, label = label), 
            inherit.aes = FALSE, size = 6, vjust = 0.8)+
  # Expand limits slightly to make room for stars
  scale_y_continuous(
    limits = c(-30, 105),
    breaks = seq(-30, 90, by = 20),
    labels = function(x) paste0(abs(x), "%")
  ) 


# FIGURE 1 ####
# 1. Prepare fig1b without its legend
fig1b_no_legend <- fig1b + theme(legend.position = "none")

# 2. Extract the legend as a separate object
# (Using the 'get_legend' function from cowplot or a patchwork hack)
legend_b <- cowplot::get_legend(fig1b)

# 3. Plot 
## Define layout
design <- "
  AT
  BL
"

# Combine using the wrap_plots function (cleaner for design layouts)
fig1 <- wrap_plots(
  A = fig1a, 
  T = table_fig1a, 
  B = fig1b_no_legend, 
  L = legend_b, 
  design = design
) +
  plot_layout(widths = c(1, 0.6, 0.8, 1)) # Forces the table area to be narrower than the plot

ggsave(here::here("figures", "fig1.pdf"), width = 12, height = 10)

# SUPPLEMENTARY FIGURE 4: figure 1 stratified for chronotype ####
## SUP. FIG. 4A ####
### Fit model
model_duration_chrono <- brm(
  sleep_duration ~ cluster_k5 * chronotype_3cat + (1 | id),
  data = data,
  family = student(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 2.5)", class = "sd"),
    set_prior("student_t(3, 0, 5)", class = "sigma")
  ),
  seed = 647,
  chains = 4, iter = 4000, warmup = 2000, cores = 4,
  backend = "cmdstanr",
)

summary(model_duration_chrono)
# sigma very similar (from 1.09 to 1.07), some credible interactions
pp_check(model_duration_chrono)
# looks fine

### Table of differences
contrast_duration_chrono <- emmeans(
  model_duration_chrono, 
  ~ chronotype_3cat | cluster_k5
  )

#### Extract means
dur_means_chrono_tb <- contrast_duration_chrono |>
  as_tibble()

#### Extract pds
contrasts_duration_within_shift <- contrast(
  contrast_duration_chrono,
  method = "trt.vs.ctrl", # pick first level (Day shifts) as baseline
  ref = 1)

pd_duration_contrasts <- describe_posterior(
  contrasts_duration_within_shift ,
  test = "p_direction"
) |> 
  as_tibble() |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
# early chronotypes slept longer than intermediate and late during day shifts 
# and long day shifts, also better than intermediate in long evening shifts
# late chronotypes slept longer than early in evening shifts
# no credible differences during night shifts

### Plot
### 1. Generate predictions including chronotype
reference_grid_chrono <- expand_grid(
  cluster_k5 = levels(data$cluster_k5),
  chronotype_3cat = levels(data$chronotype_3cat)
)
model_preds_chrono <- reference_grid_chrono |>
  add_predicted_draws(model_duration_chrono, re_formula = NA)

# 2. Plot with Facets
supfig_4a <- ggplot(
  data, 
  aes(
    x = chronotype_3cat, 
    y = sleep_duration, 
    fill = chronotype_3cat
    )
  ) +
  geom_violin(
    width = 1.0, 
    trim = FALSE
    ) +
  stat_pointinterval(
    data = model_preds_chrono, 
    aes(y = .prediction)
    ) + 
  theme_minimal() +
  theme(
    legend.position = "right", 
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 10) # facet title
      ) +
  scale_fill_manual(values = chronotype_colours) +
  labs(
    title = "A)",
    x = "Chronotype",
    y = "Sleep duration (hours)",
    fill = "Chronotype"
    )+
  facet_wrap(
    ~ factor(cluster_k5, 
             levels = c("Day shifts", "Evening shifts", "Night shifts", 
                        "Long day shifts", "Long evening shifts")),
    nrow = 1, ncol = 5
  )

### Annotate the plot
annotations_chrono_duration <- data.frame(
  cluster_k5 = factor(rep(c("Day shifts", "Evening shifts", "Night shifts", 
                            "Long day shifts", "Long evening shifts"), each = 3),
                      levels = c("Day shifts", "Evening shifts", "Night shifts", 
                                 "Long day shifts", "Long evening shifts")),
  chronotype_3cat = factor(rep(c("Early", "Intermediate", "Late"), 5),
                           levels = c("Early", "Intermediate", "Late")),
  label = c(
    "Ref.", "**", "**",  # Day: Early is ref, others slept less
    "Ref.", "", "*",   # Evening: Late slept longer
    "Ref.", "", "",  # Night: No differences
    "Ref.", "**", "**",  # Long Day: Early is ref, others slept less
    "Ref.", "**", ""   # Long Evening: Early > Interm
  )
)
supfig_4a <- supfig_4a +
  geom_text(data = annotations_chrono_duration, aes(x = chronotype_3cat, 
                                                    y = 14, 
                                                    label = label), 
            inherit.aes = FALSE, size = 3.5, vjust = 0.8)

## SUP. FIG. 4B ####
### Fit model
model_quality_chrono <- brm(
  sleep_quality ~ cluster_k5 * chronotype_3cat + (1 | id),
  data = data,
  family = cumulative("logit")
)

summary(model_quality_chrono)
# Some credible interactions
pp_check(model_quality_chrono)
# looks fine

### Get pds within each cluster compared to early chronotypes
### 1. Get means
eb_chrono <- emmeans(model_quality_chrono, ~ chronotype_3cat | cluster_k5)

### 2. Compare each chronotype to early within each shift cluster
within_shift_comparisons <- contrast(
  eb_chrono, 
  method = "trt.vs.ctrl", # pick first level (Day shifts) as baseline
  ref = 1
  )

# 3. Get pd and CIs
within_shift_results <- describe_posterior(
  within_shift_comparisons,
  test = "p_direction",
  centrality = "median"
)
within_shift_results
# Early chronotypes reported better sleep quality on day, long day and long 
# evening shifts than intermediate and late chronotypes (pd > 0.98), 
# they also reported better sleep quality on evening shifts than intermediate 
# chronotypes (pd > 0.98) 
# Late chronotypes reported better sleep quality on evening and night shifts 
# than early chronotypes (pd > 0.99)

### Plot
### 1. Create a grid with all combinations
grid_chrono <- expand.grid(
  cluster_k5 = levels(data$cluster_k5),
  chronotype_3cat = levels(data$chronotype_3cat)
)

### 2. Get the expected probabilities
probs_array <- posterior_epred(
  model_quality_chrono, 
  newdata = grid_chrono, 
  re_formula = NA
  )

### 3. Average across draws to get  "Predicted Probability"
probs_mean <- apply(probs_array, c(2, 3), mean)

### 4. Convert to a long data frame for ggplot
supfig4b_plot_data <- as.data.frame(probs_mean) |>
  # setNames(c("1", "2", "3", "4", "5")) |>
  bind_cols(grid_chrono)|>
  pivot_longer(cols = "very poor":"very good", 
               names_to = "sleep_quality", 
               values_to = "estimate")|>
  # Direction for the diverging bar
  mutate(
    sleep_quality = as.factor(sleep_quality),
    plot_estimate = case_when(
      sleep_quality %in% c("very poor", "poor") ~ -estimate*100,
      sleep_quality %in% c("moderate", "good", "very good") ~ estimate*100
    ),
    sleep_quality = fct_relevel(
    sleep_quality,
    "very poor", "poor", "very good", "good", "moderate"
  )
  )
#### Capitalize only the first letter of each level
levels(supfig4b_plot_data$sleep_quality) <- sub(
  "(^.)", "\\U\\1", 
  levels(supfig4b_plot_data$sleep_quality), 
  perl = TRUE
  )

### 5. Plot with facets for shifts clusters
supfig_4b <- ggplot(
  supfig4b_plot_data, 
  aes(x = chronotype_3cat, 
      y = plot_estimate, 
      fill = sleep_quality)
  ) +
  geom_col(width = 0.7,
           colour = "darkgrey") +
  geom_hline(yintercept = 0, color = "black") +
  scale_fill_manual(
    name = "Sleep quality",
    values = sleep_qual_colours,
    breaks = c("Very good", "Good", "Moderate", "Poor", "Very poor")
  )+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    strip.text = element_text(size = 10) # facet title
  ) +
  labs(title = "B)", 
       x = "Chronotype", 
       y = "Adjusted percentage of shifts", 
       fill = "Sleep quality")+
  facet_wrap(~ cluster_k5, nrow = 1, ncol = 5, drop = FALSE) +
  scale_y_continuous(
    limits = c(-30, 105),
    breaks = seq(-30, 90, by = 20),
    labels = function(x) paste0(abs(x), "%")
  ) 

### Annotate
annotations_chrono_quality <- data.frame(
  cluster_k5 = factor(rep(c("Day shifts", "Evening shifts", "Night shifts", 
                            "Long day shifts", "Long evening shifts"), each = 3),
                      levels = c("Day shifts", "Evening shifts", "Night shifts", 
                                 "Long day shifts", "Long evening shifts")),
  chronotype_3cat = factor(rep(c("Early", "Intermediate", "Late"), 5),
                           levels = c("Early", "Intermediate", "Late")),
  label = c(
    "Ref.", "**", "**",  # Day: Early is best
    "Ref.", "**", "*",   # Evening: Early > Interm, but Late > Early
    "Ref.", "", "*",   # Night: Late is best, no diff for Interm
    "Ref.", "**", "**",  # Long Day: Early is best
    "Ref.", "**", "**"   # Long Evening: Early is best
  )
)

supfig_4b <- supfig_4b +
  geom_text(data = annotations_chrono_quality, aes(x = chronotype_3cat, 
                                                    y = 100, 
                                                    label = label), 
            inherit.aes = FALSE, size = 3.5, vjust = 0.8)

## SUP. FIG. 4 ####
supfig_4 <- supfig_4a / supfig_4b
ggsave(here::here("figures", "supfig4.pdf"), plot = supfig_4, width = 12, height = 8)

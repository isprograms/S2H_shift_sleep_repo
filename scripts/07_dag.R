################################################################################
###################### DAG: SUPPLEMENTARY FIGURE S3 #############################
################################################################################

# by: Isabel Santonja
# description: Directed acyclic graph (DAG) of the modelled relationships
# between work schedules and sleep, used to select confounders and mediators
# for the regression models (Sup. Figure S3). Confounders (age, gender,
# country, education, marital status) are modelled as common causes of
# schedule and sleep. Income satisfaction and chronotype are modelled as
# variables that could be either confounders or mediators (bidirected edge
# to schedule). Smoking status, BMI and alcohol consumption are modelled as
# mediators on the causal pathway from schedule to sleep.

# last update: 18/09/2026

library(dagitty)

# Define DAG
dag <- dagitty("dag{
  sleep <- shift_scheduling

  sleep <- age_scaled -> shift_scheduling
  sleep <- gender_grouped -> shift_scheduling
  sleep <- country_grouped -> shift_scheduling
  sleep <- education_grouped -> shift_scheduling
  sleep <- marital_status_grouped -> shift_scheduling

  sleep <- income_satisfaction <-> shift_scheduling
  sleep <- chronotype <-> shift_scheduling

  sleep <- bmi_scaled <- shift_scheduling
  sleep <- smoking_status <- shift_scheduling
  sleep <- alcohol <- shift_scheduling
}")

coordinates(dag) <- list(
  x = c(
    shift_scheduling = 0, sleep = 2.2,

    age_scaled = 0.1, gender_grouped = 0.5,
    country_grouped = 0.9, education_grouped = 1.3,
    marital_status_grouped = 1.7,

    income_satisfaction = 0.6, chronotype = 1.6,

    bmi_scaled = 0.4, smoking_status = 1.1, alcohol = 1.8
  ),
  y = c(
    sleep = 0, shift_scheduling = 0,

    age_scaled = -1, gender_grouped = -1,
    country_grouped = -1, education_grouped = -1,
    marital_status_grouped = -1,

    income_satisfaction = -0.5, chronotype = -0.5,

    bmi_scaled = 0.6, smoking_status = 0.6, alcohol = 0.6
  )
)

pdf(here::here("figures", "sfig3.pdf"), width = 12, height = 3)
plot(dag)
dev.off()

# Programmatically derive the adjustment sets implied by the DAG (used to
# document/cross-check the confounders and mediators listed in the Methods
# and used as covariates)
predictors <- unlist(
  adjustmentSets(
    dag,
    exposure = "shift_scheduling",
    outcome = "sleep",
    effect = "direct"
  )
)

confounders <- unlist(
  adjustmentSets(
    dag,
    exposure = "shift_scheduling",
    outcome = "sleep",
    effect = "total"
  )
)

mediators <- unlist(
  setdiff(
    predictors,
    confounders
  )
)


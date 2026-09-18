################################################################################
######################## SAMPLE DESCRIPTION: TABLE S2 ###########################
################################################################################

# by: Isabel Santonja
# description: Sociodemographic and health-related characteristics of the
# study sample (Sup. Table S2), computed on the participant-level model
# dataset built by 05_prepare_model_data.R.
# last updated: 18/09/2026

library(table1)

data <- readRDS(here::here("data", "brm_data.RDS"))

table1::table1(~ oq_age + gender_grouped + country_grouped +
                 work_sector_grouped + education_grouped +
                 income_satisfaction + marital_status_grouped +
                 chronotype + alcohol + bmi_imputed + smoking_status +
                 night_shifts_updated + years_worked_ns,
               data = data)
# NOTE: the published Table S2 shows only age, gender, country, work sector,
# education, income satisfaction, marital status, chronotype, BMI, alcohol
# consumption and smoking history. 

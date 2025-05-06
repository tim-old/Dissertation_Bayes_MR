# Load packages
library(tidyverse)
library(TwoSampleMR)
library(rstan)
library(here)

# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))
source(here("MSc_Thesis_Split", "Script", "Hevo", "functions.mrhevo.R"))

# Set number of datasets used
n <- 1000

# Load data
causal_20k_point1_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point1_scen3_models.rds"))
causal_20k_point2_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point2_scen3_models.rds"))
causal_20k_point3_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point3_scen3_models.rds"))

# Generate rows
causal_20k_point1_scen3 <- get_summary_MR_tib_row(causal_20k_point1_scen3_models[1:n])
saveRDS(causal_20k_point1_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point1_scen3_row.rds"))

causal_20k_point2_scen3 <- get_summary_MR_tib_row(causal_20k_point2_scen3_models[1:n])
saveRDS(causal_20k_point2_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point2_scen3_row.rds"))

causal_20k_point3_scen3 <- get_summary_MR_tib_row(causal_20k_point3_scen3_models[1:n])
saveRDS(causal_20k_point3_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point3_scen3_row.rds"))

# Combine rows
causal_20k_scen3_sim_summ_tib <- bind_rows(causal_20k_point1_scen3,
                                           causal_20k_point2_scen3,
                                           causal_20k_point3_scen3)

# Save
saveRDS(causal_20k_scen3_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_scen3_sim_summ_tib.rds"))



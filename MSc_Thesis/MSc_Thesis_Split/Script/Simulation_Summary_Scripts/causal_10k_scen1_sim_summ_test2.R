# Load packages
library(tidyverse)
library(TwoSampleMR)
library(rstan)
library(here)

# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions_test.R"))
source(here("MSc_Thesis_Split", "Script", "Hevo", "functions.mrhevo.R"))

# Set number of datasets used
n <- 100

# Load data
#causal_10k_point0_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets_Test", "causal_10k_point0_scen3_models.rds"))
causal_10k_point1_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets_Test", "causal_10k_point1_scen3_models.rds"))
causal_10k_point2_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets_Test", "causal_10k_point2_scen3_models.rds"))
causal_10k_point3_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets_Test", "causal_10k_point3_scen3_models.rds"))

# Generate rows
#causal_10k_point0_scen3 <- get_summary_MR_tib_row(causal_10k_point0_scen3_models[1:n])
#saveRDS(causal_10k_point0_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point0_scen3_row_test.rds"))

causal_10k_point1_scen3 <- get_summary_MR_tib_row(causal_10k_point1_scen3_models[1:n])
#saveRDS(causal_10k_point1_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point1_scen3_row_test.rds"))

causal_10k_point2_scen3 <- get_summary_MR_tib_row(causal_10k_point2_scen3_models[1:n])
#saveRDS(causal_10k_point2_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point2_scen3_row_test.rds"))

causal_10k_point3_scen3 <- get_summary_MR_tib_row(causal_10k_point3_scen3_models[1:n])
#saveRDS(causal_10k_point3_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point3_scen3_row_test.rds"))

# Combine rows
causal_10k_scen3_sim_summ_tib <- bind_rows(#causal_10k_point0_scen3,
                                           causal_10k_point1_scen3,
                                           causal_10k_point2_scen3,
                                           causal_10k_point3_scen3)

# Save
#saveRDS(causal_10k_scen3_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_scen3_sim_summ_tib_test.rds"))



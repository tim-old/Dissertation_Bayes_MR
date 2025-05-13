# Load packages
library(tidyverse)
library(TwoSampleMR)
library(rstan)
library(here)

# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))

# Set number of datasets used
n <- 3

# Load data
causal_20k_point0_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point0_scen1_models.rds"))
causal_20k_point1_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point1_scen1_models.rds"))
causal_20k_point2_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point2_scen1_models.rds"))
causal_20k_point3_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point3_scen1_models.rds"))

# Generate rows
causal_20k_point0_scen1 <- get_summary_MR_tib_row(causal_20k_point0_scen1_models[1:n])
saveRDS(causal_20k_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point0_scen1_row.rds"))

causal_20k_point1_scen1 <- get_summary_MR_tib_row(causal_20k_point1_scen1_models[1:n])
saveRDS(causal_20k_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point1_scen1_row.rds"))

causal_20k_point2_scen1 <- get_summary_MR_tib_row(causal_20k_point2_scen1_models[1:n])
saveRDS(causal_20k_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point2_scen1_row.rds"))

causal_20k_point3_scen1 <- get_summary_MR_tib_row(causal_20k_point3_scen1_models[1:n])
saveRDS(causal_20k_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point3_scen1_row.rds"))

# Combine rows
causal_20k_scen1_sim_summ_tib <- bind_rows(causal_20k_point0_scen1,
                                           causal_20k_point1_scen1,
                                           causal_20k_point2_scen1,
                                           causal_20k_point3_scen1)

# Save
saveRDS(causal_20k_scen1_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_scen1_sim_summ_tib.rds"))



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
no_causal_20k_point0_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_point0_scen2_models.rds"))
no_causal_20k_point1_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_point1_scen2_models.rds"))
no_causal_20k_point2_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_point2_scen2_models.rds"))
no_causal_20k_point3_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_point3_scen2_models.rds"))

# Generate rows
no_causal_20k_point0_scen2 <- get_summary_MR_tib_row(no_causal_20k_point0_scen2_models[1:n])
saveRDS(no_causal_20k_point0_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "no_causal_20k_point0_scen2_row.rds"))

no_causal_20k_point1_scen2 <- get_summary_MR_tib_row(no_causal_20k_point1_scen2_models[1:n])
saveRDS(no_causal_20k_point1_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "no_causal_20k_point1_scen2_row.rds"))

no_causal_20k_point2_scen2 <- get_summary_MR_tib_row(no_causal_20k_point2_scen2_models[1:n])
saveRDS(no_causal_20k_point2_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "no_causal_20k_point2_scen2_row.rds"))

no_causal_20k_point3_scen2 <- get_summary_MR_tib_row(no_causal_20k_point3_scen2_models[1:n])
saveRDS(no_causal_20k_point3_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "no_causal_20k_point3_scen2_row.rds"))

# Combine rows
no_causal_20k_scen2_sim_summ_tib <- bind_rows(no_causal_20k_point0_scen2,
                                              no_causal_20k_point1_scen2,
                                              no_causal_20k_point2_scen2,
                                              no_causal_20k_point3_scen2)

# Save
saveRDS(no_causal_20k_scen2_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "no_causal_20k_scen2_sim_summ_tib.rds"))



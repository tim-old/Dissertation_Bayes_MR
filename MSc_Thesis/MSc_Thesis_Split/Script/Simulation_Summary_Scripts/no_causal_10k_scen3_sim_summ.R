# Load packages
library(tidyverse)
library(TwoSampleMR)
library(rstan)
library(here)

# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))
source(here("MSc_Thesis_Split", "Script", "Hevo", "functions.mrhevo.R"))

# Set number of datasets used
n <- 3

# Load data
no_causal_10k_point1_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_10k_point1_scen3_models.rds"))
no_causal_10k_point2_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_10k_point2_scen3_models.rds"))
no_causal_10k_point3_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_10k_point3_scen3_models.rds"))

# Generate rows
no_causal_10k_point1_scen3 <- get_summary_MR_tib_row(no_causal_10k_point1_scen3_models[1:n])
no_causal_10k_point2_scen3 <- get_summary_MR_tib_row(no_causal_10k_point2_scen3_models[1:n])
no_causal_10k_point3_scen3 <- get_summary_MR_tib_row(no_causal_10k_point3_scen3_models[1:n])

# Combine rows
no_causal_10k_scen3_sim_summ_tib <- bind_rows(no_causal_10k_point1_scen3,
                                              no_causal_10k_point2_scen3,
                                              no_causal_10k_point3_scen3)

# Save
saveRDS(no_causal_10k_scen3_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "no_causal_10k_scen3_sim_summ_tib.rds"))



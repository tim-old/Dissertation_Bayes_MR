# Load packages
library(tidyverse)
library(TwoSampleMR)
library(rstan)
library(here)

# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))

# Set number of datasets used
n <- 1000

# Load data
causal_10k_point0_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point0_scen3_models.rds"))
causal_10k_point1_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point1_scen3_models.rds"))
causal_10k_point2_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point2_scen3_models.rds"))
causal_10k_point3_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point3_scen3_models.rds"))

# Generate rows
causal_10k_point0_scen3 <- get_summary_MR_tib_row(causal_10k_point0_scen3_models[1:n])
saveRDS(causal_10k_point0_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point0_scen3_row.rds"))

causal_10k_point1_scen3 <- get_summary_MR_tib_row(causal_10k_point1_scen3_models[1:n])
saveRDS(causal_10k_point1_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point1_scen3_row.rds"))

causal_10k_point2_scen3 <- get_summary_MR_tib_row(causal_10k_point2_scen3_models[1:n])
saveRDS(causal_10k_point2_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point2_scen3_row.rds"))

causal_10k_point3_scen3 <- get_summary_MR_tib_row(causal_10k_point3_scen3_models[1:n])
saveRDS(causal_10k_point3_scen3, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_point3_scen3_row.rds"))

# Combine rows
causal_10k_scen3_sim_summ_tib <- bind_rows(causal_10k_point0_scen3,
                                           causal_10k_point1_scen3,
                                           causal_10k_point2_scen3,
                                           causal_10k_point3_scen3)


# Save
saveRDS(causal_10k_scen3_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_10k_scen3_sim_summ_tib.rds"))



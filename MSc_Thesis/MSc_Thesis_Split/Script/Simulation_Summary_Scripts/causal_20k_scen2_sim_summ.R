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

# Generate rows
#causal_20k_point0_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point0_scen2_models.rds"))
#causal_20k_point0_scen2 <- get_summary_MR_tib_row(causal_20k_point0_scen2_models[1:n])
#saveRDS(causal_20k_point0_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point0_scen2_row.rds"))
#rm(causal_20k_point0_scen2_models)

#causal_20k_point1_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point1_scen2_models.rds"))
#causal_20k_point1_scen2 <- get_summary_MR_tib_row(causal_20k_point1_scen2_models[1:n])
#saveRDS(causal_20k_point1_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point1_scen2_row.rds"))
#rm(causal_20k_point1_scen2_models)

#causal_20k_point2_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point2_scen2_models.rds"))
#causal_20k_point2_scen2 <- get_summary_MR_tib_row(causal_20k_point2_scen2_models[1:n])
#saveRDS(causal_20k_point2_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point2_scen2_row.rds"))
#rm(causal_20k_point2_scen2_models)

causal_20k_point3_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_point3_scen2_models.rds"))
causal_20k_point3_scen2 <- get_summary_MR_tib_row(causal_20k_point3_scen2_models[1:n])
saveRDS(causal_20k_point3_scen2, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_point3_scen2_row.rds"))
rm(causal_20k_point3_scen2)

# Combine rows
causal_20k_scen2_sim_summ_tib <- bind_rows(causal_20k_point0_scen2,
                                              causal_20k_point1_scen2,
                                              causal_20k_point2_scen2,
                                              causal_20k_point3_scen2)

# Save
saveRDS(causal_20k_scen2_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "causal_20k_scen2_sim_summ_tib.rds"))



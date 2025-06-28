# Load packages
library(tidyverse)
library(TwoSampleMR)
library(rstan)
library(here)

# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))

# Set number of datasets used
n <- 100

### --- Scenario 1 --- ###

# CAUSAL

# Generate rows
causal_20k_100snp_point0_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point0_scen1_models.rds"))
causal_20k_100snp_point0_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point0_scen1_models[1:n])
saveRDS(causal_20k_100snp_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point0_scen1_row.rds"))
rm(causal_20k_100snp_point0_scen1_models)

causal_20k_100snp_point1_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point1_scen1_models.rds"))
causal_20k_100snp_point1_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point1_scen1_models[1:n])
saveRDS(causal_20k_100snp_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point1_scen1_row.rds"))
rm(causal_20k_100snp_point1_scen1_models)

causal_20k_100snp_point2_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point2_scen1_models.rds"))
causal_20k_100snp_point2_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point2_scen1_models[1:n])
saveRDS(causal_20k_100snp_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point2_scen1_row.rds"))
rm(causal_20k_100snp_point2_scen1_models)

causal_20k_100snp_point3_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point3_scen1_models.rds"))
causal_20k_100snp_point3_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point3_scen1_models[1:n])
saveRDS(causal_20k_100snp_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point3_scen1_row.rds"))
rm(causal_20k_100snp_point3_scen1_models)

# Combine rows
causal_20k_100snp_scen1_sim_summ_tib <- bind_rows(causal_20k_100snp_point0_scen1,
                                            causal_20k_100snp_point1_scen1,
                                            causal_20k_100snp_point2_scen1,
                                            causal_20k_100snp_point3_scen1)

# Save
saveRDS(causal_20k_100snp_scen1_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_scen1_sim_summ_tib.rds"))
rm(causal_20k_100snp_point0_scen1,
   causal_20k_100snp_point1_scen1,
   causal_20k_100snp_point2_scen1,
   causal_20k_100snp_point3_scen1,
   causal_20k_100snp_scen1_sim_summ_tib)


# NO CAUSAL

# Generate rows
no_causal_20k_100snp_point0_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point0_scen1_models.rds"))
no_causal_20k_100snp_point0_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point0_scen1_models[1:n])
saveRDS(no_causal_20k_100snp_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point0_scen1_row.rds"))
rm(no_causal_20k_100snp_point0_scen1_models)

no_causal_20k_100snp_point1_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point1_scen1_models.rds"))
no_causal_20k_100snp_point1_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point1_scen1_models[1:n])
saveRDS(no_causal_20k_100snp_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point1_scen1_row.rds"))
rm(no_causal_20k_100snp_point1_scen1_models)

no_causal_20k_100snp_point2_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point2_scen1_models.rds"))
no_causal_20k_100snp_point2_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point2_scen1_models[1:n])
saveRDS(no_causal_20k_100snp_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point2_scen1_row.rds"))
rm(no_causal_20k_100snp_point2_scen1_models)

no_causal_20k_100snp_point3_scen1_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point3_scen1_models.rds"))
no_causal_20k_100snp_point3_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point3_scen1_models[1:n])
saveRDS(no_causal_20k_100snp_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point3_scen1_row.rds"))
rm(no_causal_20k_100snp_point3_scen1_models)

# Combine rows
no_causal_20k_100snp_scen1_sim_summ_tib <- bind_rows(no_causal_20k_100snp_point0_scen1,
                                               no_causal_20k_100snp_point1_scen1,
                                               no_causal_20k_100snp_point2_scen1,
                                               no_causal_20k_100snp_point3_scen1)

# Save
saveRDS(no_causal_20k_100snp_scen1_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_scen1_sim_summ_tib.rds"))
rm(no_causal_20k_100snp_point0_scen1,
   no_causal_20k_100snp_point1_scen1,
   no_causal_20k_100snp_point2_scen1,
   no_causal_20k_100snp_point3_scen1,
   no_causal_20k_100snp_scen1_sim_summ_tib)


### --- Scenario 2 --- ###

# CAUSAL

# Generate rows
causal_20k_100snp_point0_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point0_scen2_models.rds"))
causal_20k_100snp_point0_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point0_scen2_models[1:n])
saveRDS(causal_20k_100snp_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point0_scen2_row.rds"))
rm(causal_20k_100snp_point0_scen2_models)

causal_20k_100snp_point1_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point1_scen2_models.rds"))
causal_20k_100snp_point1_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point1_scen2_models[1:n])
saveRDS(causal_20k_100snp_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point1_scen2_row.rds"))
rm(causal_20k_100snp_point1_scen2_models)

causal_20k_100snp_point2_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point2_scen2_models.rds"))
causal_20k_100snp_point2_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point2_scen2_models[1:n])
saveRDS(causal_20k_100snp_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point2_scen2_row.rds"))
rm(causal_20k_100snp_point2_scen2_models)

causal_20k_100snp_point3_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point3_scen2_models.rds"))
causal_20k_100snp_point3_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point3_scen2_models[1:n])
saveRDS(causal_20k_100snp_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point3_scen2_row.rds"))
rm(causal_20k_100snp_point3_scen2_models)

# Combine rows
causal_20k_100snp_scen2_sim_summ_tib <- bind_rows(causal_20k_100snp_point0_scen1,
                                            causal_20k_100snp_point1_scen1,
                                            causal_20k_100snp_point2_scen1,
                                            causal_20k_100snp_point3_scen1)

# Save
saveRDS(causal_20k_100snp_scen2_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_scen2_sim_summ_tib.rds"))
rm(causal_20k_100snp_point0_scen1,
   causal_20k_100snp_point1_scen1,
   causal_20k_100snp_point2_scen1,
   causal_20k_100snp_point3_scen1,
   causal_20k_100snp_scen2_sim_summ_tib)


# NO CAUSAL

# Generate rows
no_causal_20k_100snp_point0_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point0_scen2_models.rds"))
no_causal_20k_100snp_point0_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point0_scen2_models[1:n])
saveRDS(no_causal_20k_100snp_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point0_scen2_row.rds"))
rm(no_causal_20k_100snp_point0_scen2_models)

no_causal_20k_100snp_point1_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point1_scen2_models.rds"))
no_causal_20k_100snp_point1_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point1_scen2_models[1:n])
saveRDS(no_causal_20k_100snp_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point1_scen2_row.rds"))
rm(no_causal_20k_100snp_point1_scen2_models)

no_causal_20k_100snp_point2_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point2_scen2_models.rds"))
no_causal_20k_100snp_point2_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point2_scen2_models[1:n])
saveRDS(no_causal_20k_100snp_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point2_scen2_row.rds"))
rm(no_causal_20k_100snp_point2_scen2_models)

no_causal_20k_100snp_point3_scen2_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point3_scen2_models.rds"))
no_causal_20k_100snp_point3_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point3_scen2_models[1:n])
saveRDS(no_causal_20k_100snp_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point3_scen2_row.rds"))
rm(no_causal_20k_100snp_point3_scen2_models)

# Combine rows
no_causal_20k_100snp_scen2_sim_summ_tib <- bind_rows(no_causal_20k_100snp_point0_scen1,
                                               no_causal_20k_100snp_point1_scen1,
                                               no_causal_20k_100snp_point2_scen1,
                                               no_causal_20k_100snp_point3_scen1)

# Save
saveRDS(no_causal_20k_100snp_scen2_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_scen2_sim_summ_tib.rds"))
rm(no_causal_20k_100snp_point0_scen1,
   no_causal_20k_100snp_point1_scen1,
   no_causal_20k_100snp_point2_scen1,
   no_causal_20k_100snp_point3_scen1,
   no_causal_20k_100snp_scen2_sim_summ_tib)


### --- Scenario 3 --- ###

# CAUSAL

# Generate rows
causal_20k_100snp_point0_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point0_scen3_models.rds"))
causal_20k_100snp_point0_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point0_scen3_models[1:n])
saveRDS(causal_20k_100snp_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point0_scen3_row.rds"))
rm(causal_20k_100snp_point0_scen3_models)

causal_20k_100snp_point1_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point1_scen3_models.rds"))
causal_20k_100snp_point1_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point1_scen3_models[1:n])
saveRDS(causal_20k_100snp_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point1_scen3_row.rds"))
rm(causal_20k_100snp_point1_scen3_models)

causal_20k_100snp_point2_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point2_scen3_models.rds"))
causal_20k_100snp_point2_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point2_scen3_models[1:n])
saveRDS(causal_20k_100snp_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point2_scen3_row.rds"))
rm(causal_20k_100snp_point2_scen3_models)

causal_20k_100snp_point3_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_20k_100snp_point3_scen3_models.rds"))
causal_20k_100snp_point3_scen1 <- get_summary_MR_tib_row(causal_20k_100snp_point3_scen3_models[1:n])
saveRDS(causal_20k_100snp_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_point3_scen3_row.rds"))
rm(causal_20k_100snp_point3_scen3_models)

# Combine rows
causal_20k_100snp_scen3_sim_summ_tib <- bind_rows(causal_20k_100snp_point0_scen1,
                                            causal_20k_100snp_point1_scen1,
                                            causal_20k_100snp_point2_scen1,
                                            causal_20k_100snp_point3_scen1)

# Save
saveRDS(causal_20k_100snp_scen3_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_causal_20k_100snp_scen3_sim_summ_tib.rds"))
rm(causal_20k_100snp_point0_scen1,
   causal_20k_100snp_point1_scen1,
   causal_20k_100snp_point2_scen1,
   causal_20k_100snp_point3_scen1,
   causal_20k_100snp_scen3_sim_summ_tib)


# NO CAUSAL

# Generate rows
no_causal_20k_100snp_point0_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point0_scen3_models.rds"))
no_causal_20k_100snp_point0_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point0_scen3_models[1:n])
saveRDS(no_causal_20k_100snp_point0_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point0_scen3_row.rds"))
rm(no_causal_20k_100snp_point0_scen3_models)

no_causal_20k_100snp_point1_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point1_scen3_models.rds"))
no_causal_20k_100snp_point1_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point1_scen3_models[1:n])
saveRDS(no_causal_20k_100snp_point1_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point1_scen3_row.rds"))
rm(no_causal_20k_100snp_point1_scen3_models)

no_causal_20k_100snp_point2_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point2_scen3_models.rds"))
no_causal_20k_100snp_point2_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point2_scen3_models[1:n])
saveRDS(no_causal_20k_100snp_point2_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point2_scen3_row.rds"))
rm(no_causal_20k_100snp_point2_scen3_models)

no_causal_20k_100snp_point3_scen3_models <- readRDS(file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "no_causal_20k_100snp_point3_scen3_models.rds"))
no_causal_20k_100snp_point3_scen1 <- get_summary_MR_tib_row(no_causal_20k_100snp_point3_scen3_models[1:n])
saveRDS(no_causal_20k_100snp_point3_scen1, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_point3_scen3_row.rds"))
rm(no_causal_20k_100snp_point3_scen3_models)

# Combine rows
no_causal_20k_100snp_scen3_sim_summ_tib <- bind_rows(no_causal_20k_100snp_point0_scen1,
                                               no_causal_20k_100snp_point1_scen1,
                                               no_causal_20k_100snp_point2_scen1,
                                               no_causal_20k_100snp_point3_scen1)

# Save
saveRDS(no_causal_20k_100snp_scen3_sim_summ_tib, file = here("MSc_Thesis_Split", "Data", "Summary_Tables", "sens_no_causal_20k_100snp_scen3_sim_summ_tib.rds"))
rm(no_causal_20k_100snp_point0_scen1,
   no_causal_20k_100snp_point1_scen1,
   no_causal_20k_100snp_point2_scen1,
   no_causal_20k_100snp_point3_scen1,
   no_causal_20k_100snp_scen3_sim_summ_tib)



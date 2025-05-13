# Load packages
library(tidyverse)
library(here)



# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))

# Set min datasets
n_datasets <- 1000

### ---------------------- Positive Causal Effect -------------------------- ###

# --- Scenario 1: Balanced Pleiotropy, InSIDE Assumption Satisfied --- #


## 10,000 Participants ##

# 0% Invalid

causal_10k_point0_scen1_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0,
                                                       causal_effect = TRUE, 
                                                       #two_sample = FALSE,     # for testing
                                                       #rand_error = FALSE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = TRUE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point0_scen1_models <- get_models(causal_10k_point0_scen1_data)


saveRDS(causal_10k_point0_scen1_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point0_scen1_models.rds"))

# 10% Invalid

causal_10k_point1_scen1_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.1,
                                                       causal_effect = TRUE,
                                                       #two_sample = FALSE,
                                                       #rand_error = FALSE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = TRUE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point1_scen1_models <- get_models(causal_10k_point1_scen1_data)


saveRDS(causal_10k_point1_scen1_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point1_scen1_models.rds"))

# 20% Invalid

causal_10k_point2_scen1_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.2,
                                                       causal_effect = TRUE,
                                                       #two_sample = FALSE,
                                                       #rand_error = FALSE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = TRUE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point2_scen1_models <- get_models(causal_10k_point2_scen1_data)


saveRDS(causal_10k_point2_scen1_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point2_scen1_models.rds"))

# 30% Invalid

causal_10k_point3_scen1_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.3,
                                                       causal_effect = TRUE,
                                                       #two_sample = FALSE,
                                                       #rand_error = FALSE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = TRUE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point3_scen1_models <- get_models(causal_10k_point3_scen1_data)


saveRDS(causal_10k_point3_scen1_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point3_scen1_models.rds"))

# --- Scenario 2: Directional Pleiotropy, InSIDE Assumption Satisfied --- #


## 10,000 Participants ##

# 0% Invalid

causal_10k_point0_scen2_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point0_scen2_models <- get_models(causal_10k_point0_scen2_data)


saveRDS(causal_10k_point0_scen2_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point0_scen2_models.rds"))

# 10% Invalid

causal_10k_point1_scen2_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.1,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point1_scen2_models <- get_models(causal_10k_point1_scen2_data)


saveRDS(causal_10k_point1_scen2_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point1_scen2_models.rds"))

# 20% Invalid

causal_10k_point2_scen2_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.2,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point2_scen2_models <- get_models(causal_10k_point2_scen2_data)


saveRDS(causal_10k_point2_scen2_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point2_scen2_models.rds"))

# 30% Invalid

causal_10k_point3_scen2_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.3,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = TRUE)

causal_10k_point3_scen2_models <- get_models(causal_10k_point3_scen2_data)


saveRDS(causal_10k_point3_scen2_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point3_scen2_models.rds"))

# --- Scenario 3: Directional Pleiotropy, InSIDE Assumption Not Satisfied --- #


## 10,000 Participants ##

# 0% Invalid

causal_10k_point0_scen3_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = FALSE)

causal_10k_point0_scen3_models <- get_models(causal_10k_point0_scen3_data)


saveRDS(causal_10k_point0_scen3_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point0_scen3_models.rds"))

# 10% Invalid

causal_10k_point1_scen3_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.1,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = FALSE)

causal_10k_point1_scen3_models <- get_models(causal_10k_point1_scen3_data)


saveRDS(causal_10k_point1_scen3_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point1_scen3_models.rds"))

# 20% Invalid

causal_10k_point2_scen3_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.2,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = FALSE)

causal_10k_point2_scen3_models <- get_models(causal_10k_point2_scen3_data)


saveRDS(causal_10k_point2_scen3_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point2_scen3_models.rds"))

# 30% Invalid

causal_10k_point3_scen3_data <-  get_simulated_MR_data(n_participants = 10000,
                                                       n_instruments = 25,
                                                       n_datasets = n_datasets,
                                                       prop_invalid = 0.3,
                                                       causal_effect = TRUE,
                                                       beta_val = 0.1,
                                                       balanced_pleio = FALSE,
                                                       InSIDE_satisfied = FALSE)

causal_10k_point3_scen3_models <- get_models(causal_10k_point3_scen3_data)


saveRDS(causal_10k_point3_scen3_models, file = here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point3_scen3_models.rds"))


# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(here)



# Load function scripts
source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))



### -----------------------No Causal Effect--------------------------------- ###

# --- Scenario 1: Balanced Pleiotropy, InSIDE Assumption Satisfied --- #


## 10,000 Participants ##

# 10% Invalid
set.seed(14101583)
no_causal_10k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 100,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point1_scen1_models <- extract_models(no_causal_10k_point1_scen1_data)

# 20% Invalid
set.seed(14101583)
no_causal_10k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 100,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point2_scen1_models <- extract_models(no_causal_10k_point2_scen1_data)

# # 30% Invalid
# set.seed(14101583)
# no_causal_10k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point3_scen1_models <- extract_models(no_causal_10k_point3_scen1_data)
# 
# 
# 
# ## 20,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_20k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point1_scen1_models <- extract_models(no_causal_20k_point1_scen1_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_20k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point2_scen1_models <- extract_models(no_causal_20k_point2_scen1_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_20k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point3_scen1_models <- extract_models(no_causal_20k_point3_scen1_data)
# 
# 
# 
# # --- Scenario 2: Directional Pleiotropy, InSIDE Assumption Satisfied --- #
# 
# 
# ## 10,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_10k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point1_scen2_models <- extract_models(no_causal_10k_point1_scen2_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_10k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point2_scen2_models <- extract_models(no_causal_10k_point2_scen2_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_10k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point3_scen2_models <- extract_models(no_causal_10k_point3_scen2_data)
# 
# 
# 
# ## 20,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_20k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point1_scen2_models <- extract_models(no_causal_20k_point1_scen2_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_20k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point2_scen2_models <- extract_models(no_causal_20k_point2_scen2_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_20k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point3_scen2_models <- extract_models(no_causal_20k_point3_scen2_data)
# 
# 
# 
# # --- Scenario 3: Directional Pleiotropy, InSIDE Assumption Not Satisfied --- #
# 
# 
# ## 10,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_10k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_10k_point1_scen3_models <- extract_models(no_causal_10k_point1_scen3_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_10k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_10k_point2_scen3_models <- extract_models(no_causal_10k_point2_scen3_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_10k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_10k_point3_scen3_models <- extract_models(no_causal_10k_point3_scen3_data)
# 
# 
# 
# ## 20,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_20k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_20k_point1_scen3_models <- extract_models(no_causal_20k_point1_scen3_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_20k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_20k_point2_scen3_models <- extract_models(no_causal_20k_point2_scen3_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_20k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = FALSE,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_20k_point3_scen3_models <- extract_models(no_causal_20k_point3_scen3_data)
# 
# 
# ### ---------------------- Positive Causal Effect -------------------------- ###
# 
# # --- Scenario 1: Balanced Pleiotropy, InSIDE Assumption Satisfied --- #
# 
# 
# ## 10,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_10k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point1_scen1_models <- extract_models(no_causal_10k_point1_scen1_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_10k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point2_scen1_models <- extract_models(no_causal_10k_point2_scen1_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_10k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point3_scen1_models <- extract_models(no_causal_10k_point3_scen1_data)
# 
# 
# 
# ## 20,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_20k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point1_scen1_models <- extract_models(no_causal_20k_point1_scen1_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_20k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point2_scen1_models <- extract_models(no_causal_20k_point2_scen1_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_20k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = TRUE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point3_scen1_models <- extract_models(no_causal_20k_point3_scen1_data)
# 
# 
# 
# # --- Scenario 2: Directional Pleiotropy, InSIDE Assumption Satisfied --- #
# 
# 
# ## 10,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_10k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point1_scen2_models <- extract_models(no_causal_10k_point1_scen2_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_10k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point2_scen2_models <- extract_models(no_causal_10k_point2_scen2_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_10k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_10k_point3_scen2_models <- extract_models(no_causal_10k_point3_scen2_data)
# 
# 
# 
# ## 20,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_20k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point1_scen2_models <- extract_models(no_causal_20k_point1_scen2_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_20k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point2_scen2_models <- extract_models(no_causal_20k_point2_scen2_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_20k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = TRUE)
# 
# no_causal_20k_point3_scen2_models <- extract_models(no_causal_20k_point3_scen2_data)
# 
# 
# 
# # --- Scenario 3: Directional Pleiotropy, InSIDE Assumption Not Satisfied --- #
# 
# 
# ## 10,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_10k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_10k_point1_scen3_models <- extract_models(no_causal_10k_point1_scen3_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_10k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_10k_point2_scen3_models <- extract_models(no_causal_10k_point2_scen3_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_10k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 10000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_10k_point3_scen3_models <- extract_models(no_causal_10k_point3_scen3_data)
# 
# 
# 
# ## 20,000 Participants ##
# 
# # 10% Invalid
# set.seed(14101583)
# no_causal_20k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.1,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_20k_point1_scen3_models <- extract_models(no_causal_20k_point1_scen3_data)
# 
# # 20% Invalid
# set.seed(14101583)
# no_causal_20k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.2,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_20k_point2_scen3_models <- extract_models(no_causal_20k_point2_scen3_data)
# 
# # 30% Invalid
# set.seed(14101583)
# no_causal_20k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
#                                                      n_instruments = 25,
#                                                      n_datasets = 20000,
#                                                      prop_invalid = 0.3,
#                                                      causal_effect = TRUE, 
#                                                      beta_val = 0.1,
#                                                      balanced_pleio = FALSE,
#                                                      InSIDE_satisfied = FALSE)
# 
# no_causal_20k_point3_scen3_models <- extract_models(no_causal_20k_point3_scen3_data)


save(as.list(environment()), file = here("MSc_Thesis_Split", "Data", "simulated_datasets_models.RData"))
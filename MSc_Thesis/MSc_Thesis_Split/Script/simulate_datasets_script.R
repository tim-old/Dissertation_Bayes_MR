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
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point1_scen1_models <- extract_models(no_causal_10k_point1_scen1_data)

saveRDS(no_causal_10k_point1_scen1_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point1_scen1.rds"))
 
# 20% Invalid
set.seed(14101583)
no_causal_10k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point2_scen1_models <- extract_models(no_causal_10k_point2_scen1_data)

saveRDS(no_causal_10k_point2_scen1_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point2_scen1.rds"))
 
# 30% Invalid
set.seed(14101583)
no_causal_10k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.3,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point3_scen1_models <- extract_models(no_causal_10k_point3_scen1_data)

saveRDS(no_causal_10k_point3_scen1_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point3_scen1.rds"))
 


## 20,000 Participants ##

# 10% Invalid
set.seed(14101583)
no_causal_20k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_20k_point1_scen1_models <- extract_models(no_causal_20k_point1_scen1_data)

saveRDS(no_causal_20k_point1_scen1_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point1_scen1.rds"))
 
# 20% Invalid
set.seed(14101583)
no_causal_20k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_20k_point2_scen1_models <- extract_models(no_causal_20k_point2_scen1_data)

saveRDS(no_causal_20k_point2_scen1_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point2_scen1.rds"))
 
# 30% Invalid
set.seed(14101583)
no_causal_20k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.3,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = TRUE,
                                                     InSIDE_satisfied = TRUE)

no_causal_20k_point3_scen1_models <- extract_models(no_causal_20k_point3_scen1_data)

saveRDS(no_causal_20k_point3_scen1_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point3_scen1.rds"))
 


# --- Scenario 2: Directional Pleiotropy, InSIDE Assumption Satisfied --- #


## 10,000 Participants ##

# 10% Invalid
set.seed(14101583)
no_causal_10k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point1_scen2_models <- extract_models(no_causal_10k_point1_scen2_data)

saveRDS(no_causal_10k_point1_scen2_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point1_scen2.rds"))
 
# 20% Invalid
set.seed(14101583)
no_causal_10k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point2_scen2_models <- extract_models(no_causal_10k_point2_scen2_data)

saveRDS(no_causal_10k_point2_scen2_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point2_scen2.rds"))
 
# 30% Invalid
set.seed(14101583)
no_causal_10k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.3,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = TRUE)

no_causal_10k_point3_scen2_models <- extract_models(no_causal_10k_point3_scen2_data)

saveRDS(no_causal_10k_point3_scen2_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point3_scen2.rds"))
 


## 20,000 Participants ##

# 10% Invalid
set.seed(14101583)
no_causal_20k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = TRUE)

no_causal_20k_point1_scen2_models <- extract_models(no_causal_20k_point1_scen2_data)

saveRDS(no_causal_20k_point1_scen2_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point1_scen2.rds"))
 
# 20% Invalid
set.seed(14101583)
no_causal_20k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = TRUE)

no_causal_20k_point2_scen2_models <- extract_models(no_causal_20k_point2_scen2_data)

saveRDS(no_causal_20k_point2_scen2_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point2_scen2.rds"))
 
# 30% Invalid
set.seed(14101583)
no_causal_20k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.3,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = TRUE)

no_causal_20k_point3_scen2_models <- extract_models(no_causal_20k_point3_scen2_data)

saveRDS(no_causal_20k_point3_scen2_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point3_scen2.rds"))
 


# --- Scenario 3: Directional Pleiotropy, InSIDE Assumption Not Satisfied --- #


## 10,000 Participants ##

# 10% Invalid
set.seed(14101583)
no_causal_10k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = FALSE)

no_causal_10k_point1_scen3_models <- extract_models(no_causal_10k_point1_scen3_data)

saveRDS(no_causal_10k_point1_scen3_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point1_scen3.rds"))
 
# 20% Invalid
set.seed(14101583)
no_causal_10k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = FALSE)

no_causal_10k_point2_scen3_models <- extract_models(no_causal_10k_point2_scen3_data)

saveRDS(no_causal_10k_point2_scen3_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point2_scen3.rds"))
 
# 30% Invalid
set.seed(14101583)
no_causal_10k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 10000,
                                                     prop_invalid = 0.3,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = FALSE)

no_causal_10k_point3_scen3_models <- extract_models(no_causal_10k_point3_scen3_data)

saveRDS(no_causal_10k_point3_scen3_models, file = here("MSc_Thesis_Split", "Data", "no_causal_10k_point3_scen3.rds"))
 


## 20,000 Participants ##

# 10% Invalid
set.seed(14101583)
no_causal_20k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.1,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = FALSE)

no_causal_20k_point1_scen3_models <- extract_models(no_causal_20k_point1_scen3_data)

saveRDS(no_causal_20k_point1_scen3_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point1_scen3.rds"))
 
# 20% Invalid
set.seed(14101583)
no_causal_20k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.2,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = FALSE)

no_causal_20k_point2_scen3_models <- extract_models(no_causal_20k_point2_scen3_data)

saveRDS(no_causal_20k_point2_scen3_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point2_scen3.rds"))
 
# 30% Invalid
set.seed(14101583)
no_causal_20k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                     n_instruments = 25,
                                                     n_datasets = 20000,
                                                     prop_invalid = 0.3,
                                                     causal_effect = FALSE,
                                                     balanced_pleio = FALSE,
                                                     InSIDE_satisfied = FALSE)

no_causal_20k_point3_scen3_models <- extract_models(no_causal_20k_point3_scen3_data)

saveRDS(no_causal_20k_point3_scen3_models, file = here("MSc_Thesis_Split", "Data", "no_causal_20k_point3_scen3.rds"))
 

### ---------------------- Positive Causal Effect -------------------------- ###

# --- Scenario 1: Balanced Pleiotropy, InSIDE Assumption Satisfied --- #


## 10,000 Participants ##

# 10% Invalid
set.seed(14101583)
causal_10k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.1,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = TRUE,
                                                  InSIDE_satisfied = TRUE)

causal_10k_point1_scen1_models <- extract_models(causal_10k_point1_scen1_data)

saveRDS(causal_10k_point1_scen1_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point1_scen1.rds"))

# 20% Invalid
set.seed(14101583)
causal_10k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.2,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = TRUE,
                                                  InSIDE_satisfied = TRUE)

causal_10k_point2_scen1_models <- extract_models(causal_10k_point2_scen1_data)

saveRDS(causal_10k_point2_scen1_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point2_scen1.rds"))

# 30% Invalid
set.seed(14101583)
causal_10k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.3,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = TRUE,
                                                  InSIDE_satisfied = TRUE)

causal_10k_point3_scen1_models <- extract_models(causal_10k_point3_scen1_data)

saveRDS(causal_10k_point3_scen1_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point3_scen1.rds"))



## 20,000 Participants ##

# 10% Invalid
set.seed(14101583)
causal_20k_point1_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.1,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = TRUE,
                                                  InSIDE_satisfied = TRUE)

causal_20k_point1_scen1_models <- extract_models(causal_20k_point1_scen1_data)

saveRDS(causal_20k_point1_scen1_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point1_scen1.rds"))

# 20% Invalid
set.seed(14101583)
causal_20k_point2_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.2,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = TRUE,
                                                  InSIDE_satisfied = TRUE)

causal_20k_point2_scen1_models <- extract_models(causal_20k_point2_scen1_data)

saveRDS(causal_20k_point2_scen1_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point2_scen1.rds"))

# 30% Invalid
set.seed(14101583)
causal_20k_point3_scen1_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.3,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = TRUE,
                                                  InSIDE_satisfied = TRUE)

causal_20k_point3_scen1_models <- extract_models(causal_20k_point3_scen1_data)

saveRDS(causal_20k_point3_scen1_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point3_scen1.rds"))



# --- Scenario 2: Directional Pleiotropy, InSIDE Assumption Satisfied --- #


## 10,000 Participants ##

# 10% Invalid
set.seed(14101583)
causal_10k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.1,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = TRUE)

causal_10k_point1_scen2_models <- extract_models(causal_10k_point1_scen2_data)

saveRDS(causal_10k_point1_scen2_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point1_scen2.rds"))

# 20% Invalid
set.seed(14101583)
causal_10k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.2,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = TRUE)

causal_10k_point2_scen2_models <- extract_models(causal_10k_point2_scen2_data)

saveRDS(causal_10k_point2_scen2_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point2_scen2.rds"))

# 30% Invalid
set.seed(14101583)
causal_10k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.3,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = TRUE)

causal_10k_point3_scen2_models <- extract_models(causal_10k_point3_scen2_data)

saveRDS(causal_10k_point3_scen2_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point3_scen2.rds"))



## 20,000 Participants ##

# 10% Invalid
set.seed(14101583)
causal_20k_point1_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.1,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = TRUE)

causal_20k_point1_scen2_models <- extract_models(causal_20k_point1_scen2_data)

saveRDS(causal_20k_point1_scen2_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point1_scen2.rds"))

# 20% Invalid
set.seed(14101583)
causal_20k_point2_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.2,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = TRUE)

causal_20k_point2_scen2_models <- extract_models(causal_20k_point2_scen2_data)

saveRDS(causal_20k_point2_scen2_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point2_scen2.rds"))

# 30% Invalid
set.seed(14101583)
causal_20k_point3_scen2_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.3,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = TRUE)

causal_20k_point3_scen2_models <- extract_models(causal_20k_point3_scen2_data)

saveRDS(causal_20k_point3_scen2_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point3_scen2.rds"))



# --- Scenario 3: Directional Pleiotropy, InSIDE Assumption Not Satisfied --- #


## 10,000 Participants ##

# 10% Invalid
set.seed(14101583)
causal_10k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.1,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = FALSE)

causal_10k_point1_scen3_models <- extract_models(causal_10k_point1_scen3_data)

saveRDS(causal_10k_point1_scen3_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point1_scen3.rds"))

# 20% Invalid
set.seed(14101583)
causal_10k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.2,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = FALSE)

causal_10k_point2_scen3_models <- extract_models(causal_10k_point2_scen3_data)

saveRDS(causal_10k_point2_scen3_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point2_scen3.rds"))

# 30% Invalid
set.seed(14101583)
causal_10k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 10000,
                                                  prop_invalid = 0.3,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = FALSE)

causal_10k_point3_scen3_models <- extract_models(causal_10k_point3_scen3_data)

saveRDS(causal_10k_point3_scen3_models, file = here("MSc_Thesis_Split", "Data", "causal_10k_point3_scen3.rds"))



## 20,000 Participants ##

# 10% Invalid
set.seed(14101583)
causal_20k_point1_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.1,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = FALSE)

causal_20k_point1_scen3_models <- extract_models(causal_20k_point1_scen3_data)

saveRDS(causal_20k_point1_scen3_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point1_scen3.rds"))

# 20% Invalid
set.seed(14101583)
causal_20k_point2_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.2,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = FALSE)

causal_20k_point2_scen3_models <- extract_models(causal_20k_point2_scen3_data)

saveRDS(causal_20k_point2_scen3_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point2_scen3.rds"))

# 30% Invalid
set.seed(14101583)
causal_20k_point3_scen3_data <-  simulate_MR_data(n_participants = 10000,
                                                  n_instruments = 25,
                                                  n_datasets = 20000,
                                                  prop_invalid = 0.3,
                                                  causal_effect = TRUE,
                                                  beta_val = 0.1,
                                                  balanced_pleio = FALSE,
                                                  InSIDE_satisfied = FALSE)

causal_20k_point3_scen3_models <- extract_models(causal_20k_point3_scen3_data)

saveRDS(causal_20k_point3_scen3_models, file = here("MSc_Thesis_Split", "Data", "causal_20k_point3_scen3.rds"))

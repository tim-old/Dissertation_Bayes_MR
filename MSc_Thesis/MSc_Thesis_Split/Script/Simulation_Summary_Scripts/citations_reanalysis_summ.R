library(tidyverse)
library(here)
library(TwoSampleMR)
library(rstan)


# Run local copy of MR-Hevo functions
# Not using full package due to conflicts with Windows
source(here("MSc_Thesis_Split", "Script", "Hevo", "functions.mrhevo.R"))

# Standard set-up for RStan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE, save_dso = TRUE)


# Compile model for MR-Hevo
mr.stanmodel <- stan_model(file= here("MSc_Thesis_Split", 
                                      "Script", 
                                      "Hevo", 
                                      "MRHevo_summarystats.stan"),
                           model_name="MRHevo.summarystats", 
                           verbose=FALSE,
                           save_dso = TRUE,
                           auto_write = TRUE)


# load data
Citations_Instrument_Data <- read.csv(here("MSc_Thesis_Split", 
                                           "Data", 
                                           "Citations_Datasets", 
                                           "Citations_Instrument_Data.csv")) %>% 
  as_tibble() 



# define function
get_reanalysis_tib <- function(models_in){
  
  output_tib_row <- tibble(N = as.integer(),
                           WME_Av = as.double(),
                           WME_SE = as.double(),
                           Hevo_Av = as.double(),
                           Hevo_SE = as.double()
  )
  
  # Create blank tibble to receive results of Weighted
  # Median Estimator function from MR-Base
  
  results_tib <-  tibble(
    Study_Author = as.character(),
    Study_Ref = as.character(),
    Study_DOI = as.character(),
    Study_Exposure = as.character(),
    Study_Outcome = as.character(),
    Study_Reported_Measure = as.character(),
    N = as.integer(),
    WME_est = as.double(),
    WME_se = as.double(),
    WME_pval = as.double(),
    WME_nsnp = as.integer(),
    Hevo_est = as.double(),
    Hevo_se = as.double(),
    Hevo_sd = as.double(),
    Hevo_est_lower_CI = as.double(),
    Hevo_est_upper_CI = as.double()
  )
  
  # Convert full tibble to list of 1 tibble per dataset
  model_list <- models_in %>% 
    group_by(Study_DOI) %>%
    group_split()
  
  
  # Should = 10
  n_datasets <- length(model_list)
  
  #Run WME and MR-Hevo for each dataset
  for(dataset in 1:n_datasets){
    
    # Extract study details
    results_tib[dataset, ]$Study_Author <- model_list[[dataset]]$Author[[1]]
    results_tib[dataset, ]$Study_Ref <- model_list[[dataset]]$Study_Ref[[1]]
    results_tib[dataset, ]$Study_DOI <- model_list[[dataset]]$Study_DOI[[1]]
    results_tib[dataset, ]$Study_Exposure <- model_list[[dataset]]$Exposure[[1]]
    results_tib[dataset, ]$Study_Outcome <- model_list[[dataset]]$Outcome[[1]]
    results_tib[dataset, ]$Study_Reported_Measure <- model_list[[dataset]]$Reported_Measure[[1]]
    
    # Stored as individual vectors for MR-Hevo/RStan - not
    # Tidyverse compatible
    coeff_G_X_vect <- model_list[[dataset]]$Coeff_G_X
    coeff_G_Y_vect <- model_list[[dataset]]$Coeff_G_Y
    coeff_G_X_SE_vect <- model_list[[dataset]]$Coeff_G_X_SE
    coeff_G_Y_SE_vect <- model_list[[dataset]]$Coeff_G_Y_SE
    
    
    # N.B. MR-Hevo terminology vs WME paper/other code:
    # alpha = effects of instruments on exposure, i.e. coeff_G_X
    # beta = pleiotropic effects of instruments on outcome, i.e. alpha in WME
    # gamma = effects of instruments on outcome, i.e. coeff_G_Y
    # theta = causal effect X on Y, i.e. b
    
    # Results from weighted median estimator method
    WME_results <- mr_weighted_median(b_exp = coeff_G_X_vect,
                                      b_out = coeff_G_Y_vect,
                                      se_exp = coeff_G_X_SE_vect,
                                      se_out = coeff_G_Y_SE_vect,
                                      parameters = list(nboot = 1000))
    
    # Results from MR-Hevo method
    Hevo_results <- run_mrhevo.sstats(alpha_hat = coeff_G_X_vect,
                                      se.alpha_hat = coeff_G_X_SE_vect,
                                      gamma_hat = coeff_G_Y_vect,
                                      se.gamma_hat = coeff_G_Y_SE_vect
    ) %>%
      summary()
    
    
    # Extract WME Results
    results_tib[dataset, ]$WME_est <- WME_results$b
    results_tib[dataset, ]$WME_se <- WME_results$se
    results_tib[dataset, ]$WME_pval <- WME_results$pval
    results_tib[dataset, ]$WME_nsnp <- WME_results$nsnp
    
    # Extract MR-Hevo Results
    results_tib[dataset, ]$Hevo_est <- Hevo_results$summary["theta","mean"]
    results_tib[dataset, ]$Hevo_se <- Hevo_results$summary["theta","se_mean"]
    results_tib[dataset, ]$Hevo_sd <- Hevo_results$summary["theta","sd"]
    results_tib[dataset, ]$Hevo_est_lower_CI <- Hevo_results$summary["theta","2.5%"]
    results_tib[dataset, ]$Hevo_est_upper_CI <- Hevo_results$summary["theta","97.5%"]
    
    results_tib[dataset, ]$N <- dataset
    
    
  }
  
  results_tib <- results_tib %>%
    mutate(WME_OR = exp(WME_est),
           WME_est_lower_CI = (WME_est - (1.96 * WME_se)),
           WME_est_upper_CI = (WME_est + (1.96 * WME_se)),
           WME_OR_lower_CI = exp(WME_est_lower_CI),
           WME_OR_upper_CI = exp(WME_est_upper_CI),
           WME_est_causal_detected = WME_pval < 0.05,
           #WME_est_causal_detected = (WME_est_lower_CI > 0  | WME_est_upper_CI < 0),
           #WME_OR_causal_detected = (WME_OR_lower_CI > 1  | WME_OR_upper_CI < 1),
           Hevo_OR = exp(Hevo_est),
           Hevo_OR_lower_CI = exp(Hevo_est_lower_CI),
           Hevo_OR_upper_CI = exp(Hevo_est_upper_CI),
           Hevo_OR_causal_detected = (Hevo_OR_lower_CI > 1  | Hevo_OR_upper_CI < 1),
           Hevo_est_causal_detected = (Hevo_est_lower_CI > 0  | Hevo_est_upper_CI < 0)
    )
  
  return(results_tib)
  
}

# process
citations_reanalysis_summ_tib <- Citations_Instrument_Data %>%
         # MR-Hevo unable to handle zero values
  mutate(across(Coeff_G_X:Coeff_G_Y_SE, ~replace(., . == 0, 10^-100)),
         # prevents accidental character type
         across(starts_with("Coeff_"), ~as.double(.))) %>% 
  get_reanalysis_tib() %>% 
  mutate(Discordant = (WME_est_causal_detected != Hevo_est_causal_detected))

saveRDS(citations_reanalysis_summ_tib, file = here("MSc_Thesis_Split", 
                                                   "Data", 
                                                   "Summary_Tables", 
                                                   "citations_reanalysis_summ_tib.rds"))
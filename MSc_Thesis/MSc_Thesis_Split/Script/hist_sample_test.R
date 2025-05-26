library(tidyverse)
library(rstan)
library(TwoSampleMR)

# Load function scripts
#source(here("MSc_Thesis_Split", "Script", "simulation_functions.R"))

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

# --- Need to comment out model.dir in Hevo_results --- #  


get_summary_MR_tib_row <- function(model_list){
  
  
  # Create output tibble in same format as Table 2/3 from
  # Bowden et al
  output_tib_row <- tibble(N = as.integer(),
                           Prop_Invalid = as.double(),
                           F_stat = as.double(),
                           R2_stat = as.double(),
                           WME_Av = as.double(),
                           WME_SE = as.double(),
                           WME_Pos_Rate = as.double(),
                           Hevo_Av = as.double(),
                           Hevo_SE = as.double(),
                           Hevo_Pos_Rate = as.double())
  
  n_datasets <- length(model_list)
  
  #output_tib_row$N <-  n_datasets
  
  # Create blank tibble to receive results of Weighted
  # Median Estimator function from MR-Base
  
  results_tib <-  tibble(WME_est = as.double(),
                         WME_se = as.double(),
                         WME_pval = as.double(),
                         WME_nsnp = as.integer(),
                         Hevo_est = as.double(),
                         Hevo_se = as.double(),
                         Hevo_sd = as.double(),
                         Hevo_est_lower_CI = as.double(),
                         Hevo_est_upper_CI = as.double(),
                         Hevo_causal_detected = as.logical()
  )
  
  
  # Run WME and MR-Hevo for each dataset 
  for(dataset in 1:n_datasets){
    
    # Stored as individual vectors for MR-Hevo/RStan - not
    # Tidyverse compatible
    coeff_G_X_vect <- model_list[[dataset]]$coeff_G_X
    coeff_G_Y_vect <- model_list[[dataset]]$coeff_G_Y
    coeff_G_X_SE_vect <- model_list[[dataset]]$coeff_G_X_SE
    coeff_G_Y_SE_vect <- model_list[[dataset]]$coeff_G_Y_SE
    prop_invalid <- min(model_list[[dataset]]$prop_invalid)
    F_stat <- min(model_list[[dataset]]$F_stat)
    R2_stat <- min(model_list[[dataset]]$R2_stat)
    n_instruments <- max(model_list[[dataset]]$Instrument)
    n_participants <- min(model_list[[dataset]]$n_participants)
    
    
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
    Hevo_results<- run_mrhevo.sstats(alpha_hat = coeff_G_X_vect,
                                     se.alpha_hat = coeff_G_X_SE_vect,
                                     gamma_hat = coeff_G_Y_vect,
                                     se.gamma_hat = coeff_G_Y_SE_vect) %>%
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
    
  }
  
  # Add confidence intervals/causality Boolean to MR-Hevo
  results_tib <- results_tib %>%
    mutate(
      #Hevo_est_lower_CI = (Hevo_est - (1.96 * Hevo_se)),
      #Hevo_est_upper_CI = (Hevo_est + (1.96 * Hevo_se)),
      Hevo_est_causal_detected = (Hevo_est_lower_CI > 0  | Hevo_est_upper_CI < 0)
    )
  
  
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC10616660/
  # https://mr-dictionary.mrcieu.ac.uk/term/r-squared/
  output_tib_row <- results_tib %>% 
    summarise(N = n_participants,
              Prop_Invalid = prop_invalid,
              F_stat = mean(F_stat),
              R2_stat = mean(R2_stat),
              WME_Av = mean(WME_est),
              WME_SE = mean(WME_se),
              WME_Pos_Rate = length(WME_pval[WME_pval < 0.05]) / n_datasets,
              Hevo_Av = mean(Hevo_est),
              Hevo_SE = mean(Hevo_se),
              Hevo_Lower_CI = mean(Hevo_est_lower_CI),
              Hevo_Upper_CI = mean(Hevo_est_upper_CI),
              Hevo_Pos_Rate = sum(Hevo_est_causal_detected) / n_datasets
    ) %>% 
    mutate(across(where(is.double), round, 3))
  
  #return(output_tib_row)
  return(results_tib)
  
}




set.seed(1701)
rand_samp <- sample(1:1000, 5, replace = FALSE)

causal_10k_point0_scen1_models <- readRDS(here("MSc_Thesis_Split", "Data", "Simulated_Datasets", "causal_10k_point0_scen1_models.rds"))

causal_10k_point0_scen1_sample <- causal_10k_point0_scen1_models[rand_samp] %>% 
  get_summary_MR_tib_row() #modify


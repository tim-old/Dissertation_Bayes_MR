library(tidyverse)
library(rstan)


###------------------------------Simulate MR Data-------------------------------------------###

# Define function to create data generating model
# Arguments/default values based on Bowden et al
simulate_MR_data <- function(n_participants = as.integer(), 
                             n_instruments = as.integer(),
                             n_datasets = as.integer(),
                             prop_invalid = 0.1,
                             causal_effect = TRUE,
                             balanced_pleio = TRUE,
                             InSIDE_satisfied = TRUE,
                             rand_error = TRUE,      # remove random errors, for testing
                             two_sample = TRUE,      # 1- or 2-sample MR toggle, for testing
                             beta_val = 0.1,         # size of causal effect
                             allele_freq_min = 0.01, # frequency of effect allele
                             allele_freq_max = 0.99,
                             gamma_min = 0.03,       # size of pleiotropic effects on exposure
                             gamma_max = 0.1,
                             alpha_min = -0.2,       # size of pleiotropic effects on outcome
                             alpha_max = 0.2,
                             phi_min = -0.2,         # size of additional pleiotropic effects
                             phi_max = 0.2){         # when InSIDE not satisfied
  
  
  
  # Initialise blank lists to receive datasets for
  # each of:
  #     U (vector: unmeasured confounding exposures per participant), 
  #     X (vector: exposure:outcome associations estimated per participant) 
  #     Y (vector: gene:outcome association estimated per participant), 
  #     G (Matrices: Genotype data)
  #
  #     gamma (vector: pleiotropic effects of each instrument on exposure)
  #     alpha (vector: pleiotropic effects of each instrument on outcome)
  #     phi (vector: additional pleiotropic effects of each instrument when InSIDE 
  #     assumption not satisfied)
  U_list <- list()
  X_list <- list()
  Y_list <- list()
  G_X_list <- list()
  G_Y_list <- list()
  
  gamma_list <- list()
  alpha_list <- list()
  phi_list <- list()
  beta_list <- list()
  #prop_invalid_list <- list()
  
  n_participants_list <- list()
  n_instruments_list <- list()
  n_datasets_list <- list()
  prop_invalid_list <- list()
  causal_effect_list <- list()
  balanced_pleio_list <- list()
  InSIDE_satisfied_list <- list()
  rand_error_list <- list()
  two_sample_list <- list()
  beta_val_list <- list()
  allele_freq_min_list <- list()
  allele_freq_max_list <- list()
  gamma_min_list <- list()
  gamma_max_list <- list()
  alpha_min_list <- list()
  alpha_max_list <- list()
  phi_min_list <- list()
  phi_max_list <- list()
  
  
  # --- Assign features common to all datasets --- #
  
  # size of causal effect
  beta <- if_else(causal_effect == TRUE, 
                  beta_val,
                  0)
  
  # create vector of participant indices for 1st n participants
  # i.e. participants used for estimating gene:exposure coefficient
  sample_1_ref <- 1:n_participants        
  
  
  # Default is to estimate gene:outcome coefficient from different sample
  # to gene:exposure coefficient (i.e. simulating 2-sample MR)
  # two_sample == FALSE toggles to single sample for testing simulation
  ifelse(two_sample == FALSE,
         sample_2_ref <- sample_1_ref, # 1 sample MR
         sample_2_ref <- (n_participants+1):(2*n_participants)) # 2 sample MR
  
  
  
  # --- Create separate datasets --- #
  
  # Create N datasets by simulating genotype matrices with
  # 1 row per participant, 1 column per genetic instrument
  # Use these to estimate U, X + Y
  
  for(n in 1:n_datasets){
    
    # Create error terms for U, X + Y per participant,
    # each drawn from standard normal distribution
    # unless random error turned off (for testing)
    
    ifelse(rand_error == TRUE,
           U_epsilon_vect <- rnorm(n = 2 * n_participants),
           U_epsilon_vect <- rep(0, 2 * n_participants))
    
    ifelse(rand_error == TRUE,
           X_epsilon_vect <- rnorm(n = n_participants),
           X_epsilon_vect <- rep(0, n_participants))
    
    ifelse(rand_error == TRUE,
           Y_epsilon_vect <- rnorm(n = n_participants),
           Y_epsilon_vect <- rep(0, n_participants))
    
    
    # --- Create matrix of genotypes --- #
    
    # 0 = reference, i.e. zero effect alleles, 
    # 1 = 1 effect allele, 2 = 2 effect alleles 
    
    
    # Probability of effect allele set per dataset  
    # for each instrument, default value set at  
    # random between 0.01-0.99 (i.e. both effect +
    # reference are common alleles)
    allele_freq_vect <- runif(n = n_instruments,
                              min = allele_freq_min,
                              max = allele_freq_max)
    
    
    
    # Assign genotypes by sampling from binomial distribution
    # twice (as two alleles) per participant with probability
    # equal to frequency of effect allele
    # Create twice as many genotypes as participants in sample
    # to simulate 2 sample MR, i.e. first half used to estimate
    # Gene:Exposure, second half used to estimate Gene:Outcome
    
    
    # Matrix where columns are instruments, rows are participants
    # Values 0, 1 or 2
    G_mat <- matrix(rbinom(n = 2 * n_participants * n_instruments,
                           size = 2,
                           prob = rep(allele_freq_vect, 2 * n_participants)),
                    nrow = 2 * n_participants,
                    ncol = n_instruments,
                    byrow = TRUE)
    
    
    # --- Set characteristics for each genetic instrument --- # 
    
    # Set which instruments invalid, 0 = valid, 1 = invalid
    invalid_instrument_vect <- rbinom(n = n_instruments,
                                      size = 1, 
                                      prob = prop_invalid)
    
    
    # Set genetic effects of each instrument on the exposure,
    # drawn from uniform distribution, min/max as per Bowden 
    # et al
    gamma_vect <- runif(n = n_instruments,
                        min = gamma_min,
                        max = gamma_max)
    
    
    # Set pleiotropic effects on outcome, Scenarios and 
    # min/max from Bowden et al
    alpha_vect <- double() # Pleiotropic effects of instruments on outcome
    phi_vect <- double() # Pleiotropic effects of confounders on outcome
    
    for(j in 1:n_instruments){
      ifelse(invalid_instrument_vect[j] == 0, # alpha = 0 if valid
             alpha_vect[j] <- 0,
             ifelse(balanced_pleio == TRUE,
                    alpha_vect[j] <- runif(n = n_instruments, # balanced
                                           min = alpha_min,
                                           max = alpha_max),
                    alpha_vect[j] <- runif(n = n_instruments, # directional
                                           min = 0,
                                           max = alpha_max)
             )
      )
      
      # Assign default phi = 0 unless directional pleiotropy & 
      # InSIDE assumption not satisfied & genetic instrument invalid
      if(balanced_pleio == FALSE & InSIDE_satisfied == FALSE){
        ifelse(invalid_instrument_vect[j] == 0,
               phi_vect[j] <- 0,
               phi_vect[j] <- runif(n = 1,
                                    min = phi_min,
                                    max = phi_max)
        )
        
      }
      else{
        phi_vect[j] <- 0
      }
    }
    
    
    # --- Combine Gene matrix/parameters to recreate model --- #
    
    # Create vectors of estimates for U, X and Y per individual,
    # i.e. Ui, Xi and Yi. Uses matrix inner product operator " %*%" 
    # https://stackoverflow.com/questions/22060515/the-r-operator 
    # http://matrixmultiplication.xyz/
    
    
    Ui_vect <-  G_mat %*% phi_vect + U_epsilon_vect
    
    Xi_vect <-  G_mat[sample_1_ref, ] %*% gamma_vect + 
      Ui_vect[sample_1_ref, ] + 
      X_epsilon_vect
    
    Yi_vect <-  G_mat[sample_2_ref, ] %*% alpha_vect + 
      beta * Xi_vect + 
      Ui_vect[sample_2_ref, ] + 
      Y_epsilon_vect
    
    
    # Add vectors of estimates from this dataset to lists of 
    # estimates from all datasets
    U_list[[n]] <- Ui_vect
    
    X_list[[n]] <- Xi_vect
    
    Y_list[[n]] <- Yi_vect
    
    G_X_list[[n]] <- G_mat[sample_1_ref, ]
    
    G_Y_list[[n]] <- G_mat[sample_2_ref, ]
    
    
    # Include actual parameter values generated for simulation 
    alpha_list[[n]] <- alpha_vect
    
    gamma_list[[n]] <- gamma_vect
    
    phi_list[[n]] <- phi_vect
    
    # Include inputs for reference/testing 
    n_participants_list[[n]] <- n_participants
    n_instruments_list[[n]] <- n_instruments
    #n_datasets_list[[n]] <- n_datasets
    prop_invalid_list[[n]] <- prop_invalid
    #causal_effect_list[[n]] <- causal_effect
    #balanced_pleio_list[[n]] <- balanced_pleio
    #InSIDE_satisfied_list[[n]] <- InSIDE_satisfied
    #rand_error_list[[n]] <- rand_error
    #two_sample_list[[n]] <- two_sample
    beta_val_list[[n]] <- beta_val
    #allele_freq_min_list[[n]] <- allele_freq_min
    #allele_freq_max_list[[n]] <- allele_freq_max
    #gamma_min_list[[n]] <- gamma_min
    #gamma_max_list[[n]] <- gamma_max
    #alpha_min_list[[n]] <- alpha_min
    #alpha_max_list[[n]] <- alpha_max
    #phi_min_list[[n]] <- phi_min
    #phi_max_list[[n]] <- phi_max
    
  }
  
  #     U (vector: unmeasured confounding exposures per participant), 
  #     X (vector: exposure:outcome associations estimated per participant) 
  #     Y (vector: gene:outcome association estimated per participant) 
  
  
  # --- Combine all outputs to return --- #
  
  combined_list <- list(U = U_list,         # Estimates 
                        X = X_list, 
                        Y = Y_list,
                        G_X = G_X_list,     # Genotypes of 1st sample
                        G_Y = G_Y_list,     # Genotypes of 2nd sample
                        
                        alpha = alpha_list, # Actual values for validating simulation
                        gamma = gamma_list,
                        phi = phi_list,
                        #beta = beta_list,
                        #prop_invalid = prop_invalid_list,
                        
                        n_participants = n_participants_list, # Inputs
                        n_instruments = n_instruments_list,
                        #n_datasets = n_datasets_list,
                        prop_invalid = prop_invalid_list,
                        #causal_effect = causal_effect_list,
                        #balanced_pleio = balanced_pleio_list,
                        #InSIDE_satisfied = InSIDE_satisfied_list,
                        #rand_error = rand_error_list,
                        #two_sample = two_sample_list,
                        beta_val = beta_val_list#,
                        #allele_freq_min = allele_freq_min_list,
                        #allele_freq_max = allele_freq_max_list,
                        #gamma_min = gamma_min_list,
                        #gamma_max = gamma_max_list,
                        #alpha_min = alpha_min_list,
                        #alpha_max = alpha_max_list,
                        #phi_min = phi_min_list,
                        #phi_max = phi_max_list
  )
  
  return(combined_list)
}

###------------------------------Extract Models-------------------------------------------###


# Create plotting tibble with Mean/SD X + Y grouped by
# Dataset + instrument
extract_models <- function(sim){
  
  output_list <- list()
  
  # Create linear models per dataset to get coefficients
  # for gene:exposure association (coeff_G_X) and gene:outcome
  # association (coeff_G_Y)
  for(dataset in 1:length(sim$X)){
    
    X <- sim$X[[dataset]]
    Y <- sim$Y[[dataset]]
    Instruments_X <- sim$G_X[[dataset]]
    Instruments_Y <- sim$G_Y[[dataset]]
    
    alpha <- sim$alpha[[dataset]]
    gamma <- sim$gamma[[dataset]]
    phi <- sim$phi[[dataset]]
    beta <- sim$beta_val[[dataset]]
    prop_invalid <- sim$prop_invalid[[dataset]]
    n_instruments <- sim$n_instruments[[dataset]]
    n_participants<- sim$n_participants[[dataset]]
    
    
    # Model for gene:exposure
    X_lm <- lm(X ~ 0 + Instruments_X)
    coeff_G_X_vect <- coef(summary(X_lm))[1:(ncol(Instruments_X)), 1]
    SE_coeff_G_X_vect <- coef(summary(X_lm))[1:(ncol(Instruments_X)), 2]
    
    R2_stat <- summary(lm(X ~ Instruments_X))$r.squared
    F_stat <- summary(lm(X ~ Instruments_X))$fstatistic[[1]]
    #R2_stat <- summary(X_lm)$r.squared
    #F_stat <- summary(X_lm)$fstatistic
    
    # Model for gene:outcome
    Y_lm <- lm(Y ~ 0 + Instruments_Y)
    coeff_G_Y_vect <- coef(summary(Y_lm))[1:(ncol(Instruments_Y)), 1] 
    SE_coeff_G_Y_vect <- coef(summary(Y_lm))[1:(ncol(Instruments_Y)), 2]
    
    output_list[[dataset]] <- as_tibble(list(dataset = dataset,
                                             Instrument = c(1:ncol(Instruments_X)),
                                             coeff_G_X = coeff_G_X_vect,
                                             coeff_G_X_SE = SE_coeff_G_X_vect,
                                             gamma = gamma,
                                             F_stat = F_stat,
                                             R2_stat = R2_stat,
                                             coeff_G_Y = coeff_G_Y_vect,
                                             coeff_G_Y_SE = SE_coeff_G_Y_vect,
                                             alpha = alpha,
                                             phi = phi,
                                             beta = beta,
                                             prop_invalid = prop_invalid,
                                             n_instruments = n_instruments,
                                             n_participants = n_participants),
                                        .name_repair = "unique")
  }
  
  return(output_list)
  
}  





###------------------------------Summary MR Row-------------------------------------------###

# Run local copy of MR-Hevo functions
# Not using full package due to conflicts with Windows
#source(here("MSc_Thesis_Split", "Script", "Hevo", "functions.mrhevo.R"))

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
                           auto_write = TRUE, )

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
                         WME_Q = as.double(),
                         WME_Q_df = as.double(),
                         WME_Q_pval = as.double(),
                         WME_nsnp = as.integer(),
                         Hevo_est = as.double(),
                         Hevo_se = as.double(),
                         Hevo_sd = as.double(),
                         Hevo_2.5 = as.double(),
                         Hevo_25 = as.double(),
                         Hevo_50 = as.double(),
                         Hevo_75 = as.double(),
                         Hevo_97.5 = as.double(),
                         Hevo_n_eff = as.double(),
                         Hevo_n_Rhat = as.double(),
                         Hevo_z_stat = as.double(),
                         Hevo_pval = as.double(),
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
                                     se.gamma_hat = coeff_G_Y_SE_vect
    ) %>%
      summary()
    
    
    # Extract WME Results
    results_tib[dataset, ]$WME_est <- WME_results$b
    results_tib[dataset, ]$WME_se <- WME_results$se
    results_tib[dataset, ]$WME_pval <- WME_results$pval
    results_tib[dataset, ]$WME_Q <- WME_results$Q
    results_tib[dataset, ]$WME_Q_df <- WME_results$Q_df
    results_tib[dataset, ]$WME_Q_pval <- WME_results$Q_pval
    results_tib[dataset, ]$WME_nsnp <- WME_results$nsnp
    
    # Extract MR-Hevo Results
    results_tib[dataset, ]$Hevo_est <- Hevo_results$summary["theta","mean"]
    results_tib[dataset, ]$Hevo_se <- Hevo_results$summary["theta","se_mean"]
    results_tib[dataset, ]$Hevo_sd <- Hevo_results$summary["theta","sd"]
    results_tib[dataset, ]$Hevo_2.5 <- Hevo_results$summary["theta","2.5%"]
    results_tib[dataset, ]$Hevo_25 <- Hevo_results$summary["theta","25%"]
    results_tib[dataset, ]$Hevo_50 <- Hevo_results$summary["theta","50%"]
    results_tib[dataset, ]$Hevo_75 <- Hevo_results$summary["theta","75%"]
    results_tib[dataset, ]$Hevo_97.5 <- Hevo_results$summary["theta","97.5%"]
    results_tib[dataset, ]$Hevo_n_eff <- Hevo_results$summary["theta","n_eff"]
    results_tib[dataset, ]$Hevo_n_Rhat <- Hevo_results$summary["theta","Rhat"]
    
    
  }
  
  
  results_tib <- results_tib %>%
    mutate(Hevo_causal_detected = !(Hevo_2.5 < 0  & Hevo_97.5 > 0))
  
  
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC10616660/
  # https://mr-dictionary.mrcieu.ac.uk/term/r-squared/
  output_tib_row <- results_tib %>% 
    summarise(N = n_participants,
              Prop_Invalid = prop_invalid,
              F_stat = F_stat,
              R2_stat = R2_stat,
              WME_Av = mean(WME_est),
              WME_SE = mean(WME_se),
              WME_Pos_Rate = length(WME_pval[WME_pval < 0.05]) / n_datasets,
              Hevo_Av = mean(Hevo_est),
              Hevo_SE = mean(Hevo_se),
              Hevo_Pos_Rate = sum(Hevo_causal_detected) / n_datasets
    ) %>% 
    mutate(across(where(is.double), round, 3))
  
  return(output_tib_row)
  #return(results_tib)
  
}



###------------------------------Extract Models-------------------------------------------###




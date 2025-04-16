
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
  prop_invalid_list <- list()
  
  
  # --- Assign features common to all datasets --- #
  
  beta <- if_else(causal_effect == TRUE, # size of causal effect
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
    
    
    # Include actual parameters used in simulation for testing
    alpha_list[[n]] <- alpha_vect
    
    gamma_list[[n]] <- gamma_vect
    
    phi_list[[n]] <- phi_vect
    
    beta_list[[n]] <- beta
    
    prop_invalid_list[[n]] <- prop_invalid
    
    
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
                        beta = beta_list,
                        prop_invalid = prop_invalid_list
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
    beta <- sim$beta[[dataset]]
    prop_invalid <- sim$prop_invalid[[dataset]]
    
    
    # Model for gene:exposure
    X_lm <- lm(X ~ 0 + Instruments_X)
    coeff_G_X_vect <- coef(summary(X_lm))[1:(ncol(Instruments_X)), 1]
    SE_coeff_G_X_vect <- coef(summary(X_lm))[1:(ncol(Instruments_X)), 2]
    
    # Model for gene:outcome
    Y_lm <- lm(Y ~ 0 + Instruments_Y)
    coeff_G_Y_vect <- coef(summary(Y_lm))[1:(ncol(Instruments_Y)), 1] 
    SE_coeff_G_Y_vect <- coef(summary(Y_lm))[1:(ncol(Instruments_Y)), 2]
    
    output_list[[dataset]] <- as_tibble(list(dataset = dataset,
                                             Instrument = c(1:ncol(Instruments_X)),
                                             coeff_G_X = coeff_G_X_vect,
                                             coeff_G_X_SE = SE_coeff_G_X_vect,
                                             gamma = gamma,
                                             coeff_G_Y = coeff_G_Y_vect,
                                             coeff_G_Y_SE = SE_coeff_G_Y_vect,
                                             alpha = alpha,
                                             phi = phi,
                                             beta = beta,
                                             prop_invalid = prop_invalid),
                                        .name_repair = "unique")
  }
  
  return(output_list)
  
}  




###------------------------------Summary MR Row-------------------------------------------###



###------------------------------Extract Models-------------------------------------------###




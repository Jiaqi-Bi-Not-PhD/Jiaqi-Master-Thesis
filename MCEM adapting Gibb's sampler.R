##########################################################################
################## MCEM adapting Gibb's Sampler ##########################
##########################################################################

################## Generating z as family level ##########################
check <- check |> mutate(famID_simpler = dense_rank(famID))
z_new <- c()
for (i in 1:n_fam) {
  z_new[i] = rnorm(1, mean = i, sd = 1)
  check$z[check$famID_simpler == i] <- z_new[i]
}
##########################################################################


mcem_step <- function(data, initial_theta, 
                      design, tol = 1e-6, 
                      max_iter = 1000, n_gibbs = 20,
                      frailty.dist) {
  
  current_theta <- initial_theta
  converged <- FALSE
  iter <- 0
  n_fams <- length(unique(data$famID))
  data$z <- 1 # Initialize z 
  
  theta_history <- list()
  
  ## Missing PRS data
  data$missing_prs <- is.na(data$PRS)
  
  last_imputed_state <- vector("list", m_imputations)
  for (i in 1:m_imputations) {
    last_imputed_state[[i]] <- data ## Initialize with the original data
  }
  
  while (!converged & iter < max_iter) {
    imputed_datasets <- mclapply(1:m_imputations, function(i) {
      imputed_data <- last_imputed_state[[i]] ## Start with the last state
      
      ## f(PRS)~N(...)
      lmer_model <- lmer(PRS ~ log(timeBC):BC + currentage + proband + mgeneI + (1|famID), data = imputed_data, REML = FALSE)
      sigma1 <- sigma(lmer_model)
      
      ## Iterate over families for imputation
      for (family in unique(imputed_data$famID)) { ## Family-wise
        family_data <- subset(imputed_data, famID == family)
        ## Only originally missing PRS for imputation
        missing_prs <- family_data$missing_prs
        
        if (any(missing_prs)) {
          for (j in which(missing_prs)) { ## Sampling if original PRS missing
            mu <- predict(lmer_model, newdata = family_data[j, ], re.form = NULL, type = "response")
            ## Directly update imputed_data for originally missing PRS
            imputed_data$PRS[imputed_data$famID == family & missing_prs][j] <- rnorm(1, mean = mu, sd = sigma1)
            
          } ## Sampling if original PRS missing
        } 
      } ## Family-wise
      
      return(imputed_data) # Return the updated dataset
    }, mc.cores = parallel::detectCores())
    
    ## Update the last_imputed_state for the next iteration
    last_imputed_state <- imputed_datasets
    
    ## Store the history thetas
    theta_history[[iter+1]] <- current_theta
    
    ## E-step Log-Normal
    objective_function <- function(theta) {
      mean(unlist(mclapply(imputed_datasets, function(dataset) {
        Y <- as.matrix(dataset[, c("timeBC", "BC")])
        X <- as.matrix(dataset[, c("PRS", "mgeneI")])
        lognormal_single(X = X, Y = Y, theta = theta, 
                         data = dataset, nbase = 2, 
                         base.dist = "Weibull", 
                         frailty.dist = "lognormal", agemin = 18,
                         design = design)
      }, mc.cores = parallel::detectCores())))
    }
    
    
    
    ## M-step: log-normal, Nelder-Mead
    theta_history[[iter + 1]] <- current_theta
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
    ## Stopping rule
    if (length(theta_history) >= 21) {
      squared_diffs <- sapply(1:length(current_theta), function(k) {
        (theta_history[[iter + 1]][k] - theta_history[[iter - 19]][k])^2
      })
      if (all(squared_diffs < tol)) {
        converged <- TRUE
      }
      else prev_theta <- current_theta ## Update theta
    }
    
    
    iter <- iter + 1
    print(iter)
    
  }
  
  return(list(final_theta = current_theta,
              convergence = converged,
              iterations = iter,
              Baseline = "Weibull", 
              Frailty_Dist = frailty.dist,
              theta_history = theta_history))
}
library(parallel)
##################### MCEM for Gamma & log-normal (Global Mean Imputation) #####################
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, max_iter = 1000, frailty.dist) {
  current_theta <- initial_theta
  prev_theta <- current_theta
  converged <- FALSE
  iter <- 0
  tol <- rep(tol, length(initial_theta))
  
  while (!converged & iter < max_iter) {
    ## MC samppling using X_{obs}, try Gibb's sampling
    imputed_datasets <- mclapply(1:m_imputations, function(i) {
      imputed_data <- data
      mu <- mean(data$PRS[!is.na(data$PRS)])
      sd <- sd(data$PRS[!is.na(data$PRS)])
      imputed_data$PRS[is.na(data$PRS)] <- rnorm(n = sum(is.na(data$PRS)), mean = mu, sd = sd)
      return(imputed_data)
    }, mc.cores = parallel::detectCores())
    

    ## E-step Gamma
    objective_function <- function(theta) {
      mean(unlist(mclapply(imputed_datasets, function(dataset) {
        Y <- as.matrix(dataset[, c("timeBC", "BC")])
        X <- as.matrix(dataset[, c("PRS", "mgeneI")])
        loglik_frailty_single_gamma(X = X, Y = Y, theta = theta, 
                                    data = dataset, nbase = 2, 
                                    base.dist = "Weibull", 
                                    frailty.dist = "gamma", agemin = 18,
                                    design = design)
      }, mc.cores = parallel::detectCores())))
    }
    
    ## M-step: Gamma and log-normal, Nelder-Mead
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
    ## Convergence rule
    if (sum(abs(current_theta - prev_theta) < tol) == length(initial_theta)) {
      converged <- TRUE
    } else {
      prev_theta <- current_theta
    }
    
    iter <- iter + 1
  }
  
  return(list(final_theta = current_theta,
       convergence = converged,
       iterations = iter,
       Baseline = "Weibull", 
       Frailty_Dist = frailty.dist))
}

initial_params <- c(1/41.41327,1,0, 0, 1)
results_gamma <- mcem_step(data = brca1_prs, initial_theta = initial_params,
                           design = "pop", frailty.dist = "gamma")
results_lognormal <- mcem_step(data = brca1_prs, initial_theta = initial_params, 
                     design = "pop", frailty.dist = "lognormal")

##################### MCEM for Gamma & log-normal (Gibb's Sampling, modelling each family) #####################
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, 
                      max_iter = 1000, n_gibbs = 10, burn_in = 5,
                      frailty.dist) {
  current_theta <- initial_theta
  prev_theta <- current_theta
  converged <- FALSE
  iter <- 0
  tol <- rep(tol, length(initial_theta))
  
  ## Store the history theta to plot the diagnostics
  theta_history <- list()
  
  while (!converged & iter < max_iter) {
    ## Gibb's sampling
    imputed_datasets <- mclapply(1:m_imputations, function(i) {
      imputed_data <- data
      for (g in 1:(n_gibbs + burn_in)) { ## burn-in iterations
        families <- unique(imputed_data$famID)
        for (family in families) { ## Family-wise
          family_data <- subset(imputed_data, famID == family)
          missing_prs <- is.na(family_data$PRS)
          if (any(missing_prs)) { ## Model the missing PRS
            fit <- lm(PRS ~ log(timeBC)*BC + currentage + proband + mgeneI, data = family_data[!missing_prs, ])
            sigma2 <- summary(fit)$sigma
            
            for (j in which(missing_prs)) { ## Within-family
              mu <- predict(fit, newdata = family_data[j, ], type = "response")
              ## Update PRS imputations only after burn-in period
              if (g > burn_in) {
                imputed_data$PRS[is.na(data$PRS) & data$famID == family & row.names(data) == row.names(family_data[j, ])] <- rnorm(1, mean = mu, sd = sigma2)
              }
            } ## Within-family imputation
          } ## Model the missing PRS
        } ## Family-wise
      } ## Burn-in iterations
      return(imputed_data)
    }, mc.cores = parallel::detectCores())
    
    ## E-step Log-Normal
    #objective_function <- function(theta) {
    #  mean(unlist(mclapply(imputed_datasets, function(dataset) {
    #    Y <- as.matrix(dataset[, c("timeBC", "BC")])
    #    X <- as.matrix(dataset[, c("PRS", "mgeneI")])
    #    lognormal_single(X = X, Y = Y, theta = theta, 
    #                                data = dataset, nbase = 2, 
    #                                base.dist = "Weibull", 
    #                                frailty.dist = "lognormal", agemin = 18,
    #                                design = design)
    #  }, mc.cores = parallel::detectCores())))
    #}
    
    ## E-step Gamma
    objective_function <- function(theta) {
      mean(unlist(mclapply(imputed_datasets, function(dataset) {
        Y <- as.matrix(dataset[, c("timeBC", "BC")])
        X <- as.matrix(dataset[, c("PRS", "mgeneI")])
        loglik_frailty_single_gamma(X = X, Y = Y, theta = theta, 
                                    data = dataset, nbase = 2, 
                                    base.dist = "Weibull", 
                                    frailty.dist = "gamma", agemin = 18,
                                    design = design)
      }, mc.cores = parallel::detectCores())))
    }
    
    ## M-step: Gamma and log-normal, Nelder-Mead
    theta_history[[iter + 1]] <- current_theta
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
    ## Convergence rule
    if (sum(abs(current_theta - prev_theta) < tol) == length(initial_theta)) {
      converged <- TRUE
    } else {
      prev_theta <- current_theta
    }
    
    iter <- iter + 1
  }
  
  return(list(final_theta = current_theta,
              convergence = converged,
              iterations = iter,
              Baseline = "Weibull", 
              Frailty_Dist = frailty.dist,
              theta_history = theta_history))
}

initial_params <- c(1/41.41327,1,0, 0, 1)
results_gamma <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params,
                           design = "pop", frailty.dist = "gamma", m_imputations = 5,
                           n_gibbs = 5, burn_in = 5)
results_lognormal <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params, 
                               design = "pop", frailty.dist = "lognormal")

## Diagnostics
theta_matrix <- do.call(rbind, results_gamma$theta_history)
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  plot(theta_matrix[, i], type = 'l', main = paste("Trace Plot for Parameter", i), xlab = "Iteration", ylab = "Parameter Estimate")
}
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  acf(theta_matrix[, i], main = paste("ACF for Parameter", i))
}

## E-step Gamma
#objective_function <- function(theta) {
#  mean(unlist(mclapply(imputed_datasets, function(dataset) {
#    Y <- as.matrix(dataset[, c("timeBC", "BC")])
#    X <- as.matrix(dataset[, c("PRS", "mgeneI")])
#    loglik_frailty_single_gamma(X = X, Y = Y, theta = theta, 
#                                data = dataset, nbase = 2, 
#                                base.dist = "Weibull", 
#                                frailty.dist = "gamma", agemin = 18,
#                                design = design)
#  }, mc.cores = parallel::detectCores())))
#}

##################### MCEM for Gamma & log-normal (Gibb's Sampling, modelling lmer) #####################
library(lme4)
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, 
                      max_iter = 1000, n_gibbs = 10, burn_in = 5,
                      frailty.dist) {
  current_theta <- initial_theta
  prev_theta <- current_theta
  converged <- FALSE
  iter <- 0
  tol <- rep(tol, length(initial_theta))
  
  ## Store the history theta to plot the diagnostics
  theta_history <- list()
  
  while (!converged & iter < max_iter) {
    ## Gibb's sampling
    imputed_datasets <- mclapply(1:m_imputations, function(i) {
      imputed_data <- data
      
      lmer_model <- lmer(PRS ~ log(timeBC)*BC + currentage + proband + mgeneI + (1|famID), data = imputed_data, REML = FALSE)
      sigma1 <- sigma(lmer_model)
      
      for (g in 1:(n_gibbs + burn_in)) { ## burn-in iterations
        families <- unique(imputed_data$famID)
        for (family in families) { ## Family-wise
          family_data <- subset(imputed_data, famID == family)
          missing_prs <- is.na(family_data$PRS)
          if (any(missing_prs)) { ## Model the missing PRS
            
            for (j in which(missing_prs)) { ## Within-family
              mu <- predict(lmer_model, newdata = family_data[j, ], re.form = NULL, type = "response")
              ## Update PRS imputations only after burn-in period
              if (g > burn_in) {
                imputed_data$PRS[is.na(data$PRS) & data$famID == family & row.names(data) == row.names(family_data[j, ])] <- rnorm(1, mean = mu, sd = sigma1)
              }
            } ## Within-family imputation
          } ## Model the missing PRS
        } ## Family-wise
      } ## Burn-in iterations
      return(imputed_data)
    }, mc.cores = parallel::detectCores())
    
    ## E-step Log-Normal
    #objective_function <- function(theta) {
    #  mean(unlist(mclapply(imputed_datasets, function(dataset) {
    #    Y <- as.matrix(dataset[, c("timeBC", "BC")])
    #    X <- as.matrix(dataset[, c("PRS", "mgeneI")])
    #    lognormal_single(X = X, Y = Y, theta = theta, 
    #                                data = dataset, nbase = 2, 
    #                                base.dist = "Weibull", 
    #                                frailty.dist = "lognormal", agemin = 18,
    #                                design = design)
    #  }, mc.cores = parallel::detectCores())))
    #}
    
    ## E-step Gamma
    objective_function <- function(theta) {
      mean(unlist(mclapply(imputed_datasets, function(dataset) {
        Y <- as.matrix(dataset[, c("timeBC", "BC")])
        X <- as.matrix(dataset[, c("PRS", "mgeneI")])
        loglik_frailty_single_gamma(X = X, Y = Y, theta = theta, 
                                    data = dataset, nbase = 2, 
                                    base.dist = "Weibull", 
                                    frailty.dist = "gamma", agemin = 18,
                                    design = design)
      }, mc.cores = parallel::detectCores())))
    }
    
    ## M-step: Gamma and log-normal, Nelder-Mead
    theta_history[[iter + 1]] <- current_theta
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
    ## Convergence rule
    if (sum(abs(current_theta - prev_theta) < tol) == length(initial_theta)) {
      converged <- TRUE
    } else {
      prev_theta <- current_theta
    }
    
    iter <- iter + 1
  }
  
  return(list(final_theta = current_theta,
              convergence = converged,
              iterations = iter,
              Baseline = "Weibull", 
              Frailty_Dist = frailty.dist,
              theta_history = theta_history))
}

initial_params <- c(1/41.41327,1,0, 0, 1)
results_gamma <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params,
                           design = "pop", frailty.dist = "gamma", m_imputations = 5,
                           n_gibbs = 5, burn_in = 5)
results_lognormal <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params, 
                               design = "pop", frailty.dist = "lognormal")

## Diagnostics
theta_matrix <- do.call(rbind, results_gamma$theta_history)
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  plot(theta_matrix[, i], type = 'l', main = paste("Trace Plot for Parameter", i), xlab = "Iteration", ylab = "Parameter Estimate")
}
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  acf(theta_matrix[, i], main = paste("ACF for Parameter", i))
}



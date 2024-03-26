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
      
      lmer_model <- lmer(PRS ~ log(timeBC):BC + currentage + proband + mgeneI + (1|famID), data = imputed_data, REML = FALSE)
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
    

    
    ## M-step: Gamma and log-normal, Nelder-Mead
    theta_history[[iter + 1]] <- current_theta
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
  
    prev_theta <- current_theta ## Give it no convergence rule
    
    iter <- iter + 1
  }
  
  return(list(final_theta = current_theta,
              convergence = FALSE,
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
                               design = "pop", frailty.dist = "lognormal", m_imputations = 5,
                               n_gibbs = 5, burn_in = 5)

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

## Convergence rule
#if (sum(abs(current_theta - prev_theta) < tol) == length(initial_theta)) {
#  converged <- TRUE
#} else {
#  prev_theta <- current_theta
#}

## Convergence rule (use loglikelihood value)
# if (sum(abs(current_theta - prev_theta) < tol) == length(initial_theta)) {
#  converged <- TRUE
#} else {
prev_theta <- current_theta ## Give it no convergence rule
#}

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

theta_matrix <- do.call(rbind, results_lognormal$theta_history)
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  plot(theta_matrix[, i], type = 'l', main = paste("Trace Plot for Parameter", i), xlab = "Iteration", ylab = "Parameter Estimate")
}
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  acf(theta_matrix[, i], main = paste("ACF for Parameter", i))
}

#######  MCEM for Gamma & log-normal (Gibb's Sampling, modelling lmer)  ############
library(foreach)
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, 
                      max_iter = 1000, n_gibbs = 20, burn_in = 5,
                      frailty.dist) {
  current_theta <- initial_theta
  converged <- FALSE
  iter <- 0
  
  ## Store the history theta to plot the diagnostics
  theta_history <- list()
  
  ## Initialize a structure to hold the last state of imputed data for each imputation
  last_imputed_state <- vector("list", m_imputations)
  
  imputed_data <- data
  
  while (!converged & iter < max_iter) {
    ## Gibb's sampling with conditioning on the last state
    imputed_datasets <- mclapply(1:m_imputations, function(i) {
      
      ## Fit the mixed-effects model on the current imputed data
      lmer_model <- lmer(PRS ~ log(timeBC):BC + currentage + proband + mgeneI + (1|famID), data = imputed_data, REML = FALSE)
      sigma1 <- sigma(lmer_model)
      
      ## Perform Gibbs sampling
      for (g in 1:(n_gibbs + burn_in)) { ## burn-in iterations
        foreach(family = unique(imputed_data$famID), .combine = 'rbind') %dopar% { ## Family-wise
          family_data <- subset(imputed_data, famID == family)
          missing_prs <- is.na(family_data$PRS)
          
          if (any(missing_prs)) { ## Model the missing PRS
            for (j in which(missing_prs)) { ## Within-family
              mu <- predict(lmer_model, newdata = family_data[j, ], re.form = NULL, type = "response")
              if (g > burn_in) {
                family_data$PRS[j] <- rnorm(1, mean = mu, sd = sigma1)
              }
            }
          } ## Model the missing PRS
        } ## Family-wise
      } ## burn-in iterations
      
      ## Update the last state for this imputation
      last_imputed_state[[i]] <- imputed_data
      return(imputed_data)
    }, mc.cores = parallel::detectCores())
    
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
    
    
    
    ## M-step: Gamma and log-normal, Nelder-Mead
    theta_history[[iter + 1]] <- current_theta
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
    ## Implement MC-based stopping rule here
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


initial_params <- c(1/41.41327,1,0, 0, 1)
results_gamma <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params,
                           design = "pop", frailty.dist = "gamma", m_imputations = 20,
                           n_gibbs = 5, burn_in = 5)
results_lognormal <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params, 
                               design = "pop", frailty.dist = "lognormal")
## Diagnostics
theta_matrix <- do.call(rbind, results_lognormal$theta_history)
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  plot(theta_matrix[, i], type = 'l', main = paste("Trace Plot for Parameter", i), xlab = "Iteration", ylab = "Parameter Estimate")
}
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  acf(theta_matrix[, i], main = paste("ACF for Parameter", i))
}

#######  MCEM for Gamma & log-normal (Gibb's Sampling, modelling lmer)  ############
library(foreach)
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, 
                      max_iter = 1000, n_gibbs = 20, burn_in = 5,
                      frailty.dist) {
  current_theta <- initial_theta
  converged <- FALSE
  iter <- 0
  
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
      
      ## Gibbs sampling with burn-in
      for (g in 1:(n_gibbs + burn_in)) {
        ## Iterate over families for imputation
        for (family in unique(imputed_data$famID)) { ## Family-wise
          family_data <- subset(imputed_data, famID == family)
          ## Only originally missing PRS for imputation
          missing_prs <- family_data$missing_prs
          
          if (any(missing_prs)) {
            for (j in which(missing_prs)) { ## Sampling if original PRS missing
              mu <- predict(lmer_model, newdata = family_data[j, ], re.form = NULL, type = "response")
              if (g > burn_in) {
                ## Directly update imputed_data for originally missing PRS
                imputed_data$PRS[imputed_data$famID == family & missing_prs][j] <- rnorm(1, mean = mu, sd = sigma1)
              }
            } ## Sampling if original PRS missing
          } 
        } ## Family-wise
      }
      
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
    
    
    
    ## M-step: Gamma and log-normal, Nelder-Mead
    theta_history[[iter + 1]] <- current_theta
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    
    ## Implement MC-based stopping rule here
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


initial_params <- c(1/41.41327,1,0, 0, 1)
log_normal_poss <- c(-19.0481278, 1.1189520, 0.4178628, 3.9732051, 2.8852498)
results_gamma <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params,
                           design = "pop", frailty.dist = "gamma", m_imputations = 20)
results_lognormal <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params, 
                               design = "pop", frailty.dist = "lognormal")
## Diagnostics
theta_matrix <- do.call(rbind, results_lognormal$theta_history)
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  plot(theta_matrix[, i], type = 'l', main = paste("Trace Plot for Parameter", i), xlab = "Iteration", ylab = "Parameter Estimate")
}
par(mfrow = c(3, 2)) 
for (i in 1:ncol(theta_matrix)) {
  acf(theta_matrix[, i], main = paste("ACF for Parameter", i))
}

#######  MCEM for Gamma & log-normal (Monte Carlo Imputation with EM Framework)  ############
library(foreach)
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, 
                      max_iter = 1000, n_gibbs = 20, burn_in = 5,
                      frailty.dist) {
  current_theta <- initial_theta
  converged <- FALSE
  iter <- 0
  
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
    
    
    
    ## M-step: Gamma and log-normal, Nelder-Mead
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


initial_params <- c(1/41.41327,1,0, 0, 1)
log_normal_poss <- c(-19.0481278, 1.1189520, 0.4178628, 3.9732051, 2.8852498)
results_gamma <- mcem_step(data = brca1_prs_MCEM, initial_theta = initial_params,
                           design = "pop", frailty.dist = "gamma", m_imputations = 20)
results_lognormal2 <- mcem_step(data = brca1_prs_MCEM, initial_theta = log_normal_poss, 
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

 
 ###########################################################################################################
 #######  MCEM for frailty and missing data jointly distributed  ###########################################
 #######  Adapting the missing PRS using the multivariate Normal considering the kinship matrix ############
 ###########################################################################################################
 library(foreach)
 mcem_step <- function(data, initial_theta, 
                       design, 
                       m_imputations = 20, tol = 1e-6, 
                       max_iter = 1000, n_gibbs = 20, burn_in = 5,
                       frailty.dist) {
   current_theta <- initial_theta
   converged <- FALSE
   iter <- 0
   
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

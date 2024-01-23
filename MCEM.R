##################### MCEM for Gamma (unconditional MC Sampling) #####################
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
    
    if (frailty.dist == "gamma") {
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
    }
    
    else if (frailty.dist == "lognormal") {
    ## E-step log-normal
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
    }
    
    ## M-step: Gamma and log-normal, Nelder-Mead
    optimized_result <- optim(par = current_theta, objective_function)
    current_theta <- optimized_result$par
    #optimized_result <- nloptr::nloptr(x0 = current_theta, eval_f = objective_function, 
    #                           opts = list("algorithm"="NLOPT_LN_NELDERMEAD"))
    #current_theta <- optimized_result$solution
    
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
GammaMCEM_results <- results

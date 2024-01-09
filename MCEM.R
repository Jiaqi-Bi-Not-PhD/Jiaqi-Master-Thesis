##################### MCEM #####################
mcem_step <- function(data, initial_theta, 
                      design, 
                      m_imputations = 20, tol = 1e-6, max_iter = 1000) {
  current_theta <- initial_theta
  prev_theta <- current_theta
  converged <- FALSE
  iter <- 0
  tol <- rep(tol, length(initial_theta))
  
  while (!converged & iter < max_iter) {
    ## MC samppling
    imputed_datasets <- lapply(1:m_imputations, function(i) {
      imputed_data <- data
      mu <- mean(data$PRS[!is.na(data$PRS)])
      sd <- sd(data$PRS[!is.na(data$PRS)])
      imputed_data$PRS[is.na(data$PRS)] <- rnorm(n = sum(is.na(data$PRS)), mean = mu, sd = sd)
      return(imputed_data)
    })
    
    ## E-step
    objective_function <- function(theta) {
      mean(unlist(lapply(imputed_datasets, function(dataset) {
        Y <- as.matrix(dataset[, c("timeBC", "BC")])
        X <- as.matrix(dataset[, c("PRS", "mgeneI")])
        loglik_frailty_single_gamma(X = X, Y = Y, theta = theta, 
                                    data = dataset, nbase = 2, 
                                    base.dist = "Weibull", 
                                    frailty.dist = "gamma", agemin = 18,
                                    design = design)
      })))
    }
    ## M-step
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
       iterations = iter))
}

initial_params <- c(1/41.41327,1,0.1, 0.1, 1)
results <- mcem_step(data = brca1_prs, initial_theta = initial_params, 
                     design = "pop")
GammaMCEM_results <- results

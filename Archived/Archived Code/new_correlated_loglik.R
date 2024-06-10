##### Correlated Frailty Weibull #####

# Ascertainment correction (design = pop, pop+)

#cagep <- data$currentage[ip]-agemin
#xbeta.p <- xbeta[ip]
#bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)

#H.p <- bcumhaz.p*exp(xbeta.p)
#logasc <- log(1-laplace(frailty.dist, H.p, sigma))
#logasc[logasc == -Inf] <- 0
#slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) 
#sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm = T) - slogasc

library(mvtnorm)
library(parallel)

corr_frailty_llhd <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE) {
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  nx <- dim(X)[2]
  beta <- theta[(nb+1):(nb+nx)]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[nb+nx+1])
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  eventsum <- sum(status1)
  print(eventsum)
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)
  logh <- log(bhaz) + xbeta
  loglik <- status1*logh 
  df1 <- data$df1[ip]
  Hsum <- sum(H, na.rm = TRUE)
  
  
  ##### Generate kinship matrix
  
  Sigma <- with(data, kinship2::kinship(id = indID, dadid = fatherID, momid = motherID,
                                        sex = rep(2,nrow(data))))
  
  ## MCEM
  max_iter <- 1000  
  tol <- 1e-6       
  iter <- 0         
  converged <- FALSE
  previous_llhd <- -Inf
  M <- 100  ## M samples
  
  while (!converged & iter < max_iter) {
    iter <- iter + 1
    
    if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
    
    nb <- nbase 
    
    if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
    else bparms1 <- exp(theta[1:nbase])
    beta <- theta[(nb+1):(nb+nx)]
    xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
    sigma <- exp(theta[nb+nx+1])
    
    bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
    bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
    
    H <- bcumhaz*exp(xbeta)
    logh <- log(bhaz) + xbeta
    loglik <- status1*logh 
    Hsum <- sum(H, na.rm = TRUE)
    
    L <- chol(sigma^2 * Sigma)
    
    ## E-step: Monte Carlo simulation
    mc_integrals <- replicate(M, {
      u_sample <- rnorm(nrow(data))
      z_sample <- L %*% u_sample
      z_power_d <- z_sample^eventsum

      exp_sum <- exp(sum(z_sample * Hsum))
      z_power_d * exp_sum
    })
 
    log_expected_integral <- log(mean(mc_integrals))
    
    ## M-step: Update theta using optim()
    optim_fn <- function(theta = theta) {
      beta <- theta[(nb+1):(nb+nx)]
      xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
      sigma <- exp(theta[nb+nx+1])
      
      bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
      bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
      
      H <- bcumhaz*exp(xbeta)
      logh <- log(bhaz) + xbeta
      loglik <- status1*logh 
      
      new_llhd <- sum(loglik, na.rm=TRUE) + log_expected_integral
      return(-new_llhd)  
    }
    optim_results <- optim(par = theta, fn = optim_fn, control = list(maxit = 2000))
    theta <- optim_results$par
    
    ## Current log-likelihood
    current_llhd <- -optim_results$value
    
    ## Convergence check
    if (abs(current_llhd - previous_llhd) < tol) {
      converged <- TRUE
    }
    previous_llhd <- current_llhd
  }
  
  return(list(theta = theta,
              loglikelihood = current_llhd))
}

corr_frailty_llhd(X = X, Y = Y, theta = initial_params, data = brca1_prs_cca,
                  nbase = 2, design = "pop", frailty.dist = "multivariate normal", base.dist = "Weibull",
                  agemin = 18)

# Example usage
# results <- corr_frailty_llhd(X, Y, initial_theta, ...)


################ Test Zone ################
initial_params <- c(1/41.41327,1,0,0, 1)
X <- as.matrix(brca1_prs_cca[,c("mgeneI", "PRS")], ncol = 2)
Y <- as.matrix(brca1_prs_cca[,c("timeBC", "BC")], ncol = 2)
results <- optim(par = initial_params, fn = corr_frailty_llhd,
                 data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
                 design = "pop", frailty.dist = "multivariate normal", base.dist = "Weibull",
                 agemin = 18, 
                 control = list(maxit = 2000))

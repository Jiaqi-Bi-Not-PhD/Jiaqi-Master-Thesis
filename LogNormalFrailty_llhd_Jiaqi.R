log_likelihood <- function(params, data, event, time, cluster, covariates, 
                           N_p = 15, sigma = 1) {
  
  alpha <- params[1]
  lambda <- params[2]
  beta <- params[3:(2 + length(covariates))]
  ip <- data$proband == 1
  
  ghq <- statmod::gauss.quad(N_p, kind="hermite")
  
  event <- data[[event]]
  time <- data[[time]]
  cluster <- data[[cluster]]
  covariates <- data[covariates]
  
  unique_clusters <- unique(cluster)
  
  sum_logL <- 0
  
  for (j in unique_clusters) {
    idx <- which(cluster == j)
    n_j <- length(idx)
    d_j <- sum(event[idx])
    
    product_term_i <- 1
    sum_p <- 0  # Sum over p
    
    for (i in idx) {
      xi_ij <- exp(sum(beta * covariates[i,]))
      t_ij <- time[i]
      delta_ij <- event[i]
      
      # Compute the term for the product over i
      product_term_i = product_term_i * ((alpha * lambda * (t_ij^(lambda - 1)) * xi_ij)^delta_ij) / (sqrt(2*pi)*sigma)
      
      inner_sum_p = 0
      for (p in 1:N_p) {
        s_jp <- ghq$nodes[p]
        omega_p <- ghq$weights[p]
        
        term2_exp <- sqrt(2) * s_jp * sigma / n_j
        term2 = exp(term2_exp)^(d_j - n_j)
        term3 = exp(-alpha * xi_ij * exp(term2_exp) * t_ij^lambda)
        term4 = exp(term2_exp) * (sqrt(2) * sigma / n_j)
        
        inner_sum_p = inner_sum_p + omega_p * term2 * term3 * term4
      }
      sum_p = sum_p + inner_sum_p
    }
    
    combined_term = product_term_i * sum_p
    
    sum_logL = sum_logL + log(combined_term)
  }
  
  return(-sum_logL)
}

log_likelihood(params = c(log(2), log(2), 1.117507, 2.418385), data = brca1_prs_cca,
               event = "BC", time = "timeBC", cluster = "famID", covariates = c("mgeneI", "PRS"))



######### Get rid of the for loop ##########
LogNormalFrailty_llhd2 <- function(params, data, event, time, cluster, covariates, N_p = 15, sigma = 1) {
  
  alpha <- params[1]
  lambda <- params[2]
  beta <- params[3:length(params)]
  
  ghq <- statmod::gauss.quad(15, kind="hermite")
  
  event <- data[[event]]
  time <- data[[time]]
  cluster <- data[[cluster]]
  covariates <- as.matrix(data[covariates])
  
  # Split data by cluster
  split_data <- split(1:length(cluster), cluster)
  
  compute_llhd_by_cluster <- function(idx) {
    n_j <- length(idx)
    d_j <- sum(event[idx])
    
    xi_ij <- exp(covariates[idx,] %*% beta)
    t_ij <- time[idx]
    delta_ij <- event[idx]
    
    # Compute the term for the product over i
    product_term_i <- prod((alpha * lambda * (t_ij^(lambda - 1)) * xi_ij)^delta_ij / (sqrt(2*pi)*sigma))
    
    s_jp <- ghq$nodes
    omega_p <- ghq$weights
    
    term2_exp <- matrix(sqrt(2) * s_jp * sigma / n_j, nrow=length(idx), ncol=N_p, byrow=TRUE)
    term2 <- exp(term2_exp * (d_j - n_j))
    term3_matrix <- sweep(exp(term2_exp), 1, xi_ij, "*")
    term3 <- exp(-alpha * term3_matrix * t_ij^lambda)
    term4 <- exp(term2_exp) * (sqrt(2) * sigma / n_j)
    
    sum_p <- sum(omega_p * rowSums(term2 * term3 * term4))
    log(product_term_i * sum_p)
  }
  
  # Apply function to each cluster
  sum_logL_vec <- sapply(split_data, compute_llhd_by_cluster)
  
  return(-sum(sum_logL_vec))
}
##################################################################################
########### Adjust the log-normal frailty log-likelihood from FamEvent ###########
##################################################################################
loglik_frailty_single_lognormal <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE, sigma = 1) {
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  sigma <- sigma
  
  if (base.dist=="lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx1 <- dim(X)[2]
  xbeta1 <- c(X%*%theta[(nb+1):length(theta)])
  #kappa <- exp(theta[(nb+nx1+1):length(theta)])
  
  time01 <- Y[,1] - agemin
  #time01 <- Y[,1] - agemin
  #cuts0 <- cuts - agemin
  status1 <- Y[,2]
  d <- aggregate(Y[,2], list(data$famID), sum)[,2]
  n <- data |>
    dplyr::group_by(famID) |>
    dplyr::summarise(n = n())
  n <- n$n
  
  
  ip <- data$proband == 1
  #wt <- data$weight
  #wt.p <- wt[ip]
  
  #bhaz1 <- hazards(base.dist, time01, bparms1, cuts=cuts0)
  #bcumhaz1 <- cumhaz(base.dist, time01, bparms1, cuts=cuts0)
  #bhaz1 <- hazards(base.dist, time01, bparms1)
  bhaz1 <- hazards(base.dist, time01, bparms1)
  bcumhaz1 <- cumhaz(base.dist, time01, bparms1)
  
  H1 <- bcumhaz1*exp(xbeta1)
  logh1 <- log(bhaz1) + xbeta1
  first_term <-  status1*logh1 - H1
  sum_first_term <- sum(first_term, na.rm = TRUE)
  
  df1 <- data$df1[ip]
  Hfam1 <- aggregate(H1, list(data$famID), sum)[,2]
  
  #### numerical integral univariate normal
  #pts <- gauss.hermite(10)
  #G_first <- matrix(pts$Points) %*% (d-n+1)*sqrt(2) * sigma
  #G_second <- exp((matrix(pts$Points) * sqrt(2) * sigma) %*% n^(-1))
  #G_inside <- sweep(G_first, 1, G_second, "+")
  #G_inside <- exp(G_inside)
  #print(G_first)

  
  
  
  #gfun <- function(x, d1, H1) exp(x*d1-exp(x)*H1)
  #sigma <- kappa
  #pts <- gauss.hermite(10) #6(n=32); 5(n=21)
  #logdL <- wt.p*log(apply(pts$points, 1, gfun, d1=df1, H1=Hfam1) %*% pts$weights)
  #logdL <- wt.p*log(apply(as.matrix(pts$Points), 1, 
  #                        gfun, d1=df1, H1=Hfam1) %*% as.matrix(pts$Weights))
  #logdL <- log(apply(as.matrix(pts$Points), 1, 
  #                        gfun, d1=df1, H1=Hfam1) %*% as.matrix(pts$Weights))
  
  
  # Ascertainment correction (design = pop, pop+)
  #cagep <- data$currentage[ip] - agemin
  #cagep <- data$currentage[ip]-agemin
  #xbeta1.p <- xbeta1[ip]
  #bcumhaz1.p <- cumhaz(base.dist, cagep, bparms1)
  #H1.p <- bcumhaz1.p*exp(xbeta1.p)
  
  #sfun <- function(x, H1) exp(-exp(x)*H1)
  #int.gh <- c(t(apply(pts$points, 1, sfun, H1=H1.p) %*% pts$weights))
  #int.gh <- c(t(apply(as.matrix(pts$Points), 1, sfun, H1=H1.p) %*% as.matrix(pts$Weights)))
  #logasc <- log(1-int.gh)
  
  #print(int.gh)
  #print(logasc)
  
  #slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) 
  #sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL, na.rm=T) - slogasc
  #loglik[ip] <- loglik[ip] + logdL - logasc
  
  #if(vec) return(-loglik)
  #else  return(-sloglik)
}
loglik_frailty_single_lognormal(X=X, Y=Y, theta = c(log(2), log(2), 0, 0), nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = min(brca1_prs_cca$timeBC[brca1_prs_cca$BC == 1]), sigma = 1)

## To do list: 
## Ascertainment correction: Sampling is not totally random => conditional pdf => adjustments for ascertainment
## Change the for loop
## Change the coding to Uni's FamEvent log-likelihood with ascertainment correction
## Use MCMC to approximate log-normal log-likelihood (How many samples?) => MCEM using the same sample
## * Type everything to LaTeX *
## Changchang Xu thesis uoft 

lognormal_single <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE, sigma = 1)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- sum(nbase) # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)
  logh <- log(bhaz) + xbeta
  loglik <- status1*logh - H
  
  df1 <- data$df1[ip]
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  
  #print(Hfam)
  #### numerical integral ####
  gfun <- function(x, d, H) exp(x*d-exp(x)*H)
  pts <- gauss.hermite(20) 
  logdL <- log(apply(matrix(pts$Points, nrow = length(pts$Points), ncol = 1), 
                     1, gfun, d=df1, H=Hfam) %*% pts$Weights)
  #print(logdL)
  # Ascertainment correction (design = pop, pop+)
  
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)

  H.p <- bcumhaz.p*exp(xbeta.p)
  
  
  sfun <- function(x, H) exp(-exp(x)*H)
  int.gh <- c(t(apply(matrix(pts$Points, nrow = length(pts$Points), ncol = 1), 1, sfun, H=H.p) %*% pts$Weights))
  logasc <- log(1-int.gh)
  
  slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) 
  sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL, na.rm=T) - slogasc
  #print(sum(loglik[loglik!=-Inf], na.rm = T))
  #print(logdL)
  return(-sloglik)
}

lognormal_single(X = X, Y = Y, theta = c(-4.1151644, 1.0977276, 0.2378448, 1.3960671), 
                 nbase = 2, data = brca1_prs_cca,
                 design = "pop", base.dist = "Weibull", frailty.dist = "lognormal",
                 agemin = 18)

optim(par = initial_params, fn = loglik_frailty_single_gamma,
      data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
      agemin = 18, k = 1, 
      control = list(maxit = 2000))

library(plotly)
param1_values <- seq(from = -7, to = 0, by = 0.05)
param2_values <- seq(from = 0, to = 7, by = 0.05)
param_grid <- expand.grid(param1 = param1_values, param2 = param2_values)
fixed_beta1 <- 0.2
fixed_beta2 <- -17.13
likelihood_values <- apply(param_grid, 1, function(params) {
  lognormal_single(X, Y, theta = c(params['param1'], params['param2'], fixed_beta1, fixed_beta2), nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 18)
})
likelihood_values <- -log(likelihood_values)
likelihood_matrix <- matrix(likelihood_values, nrow = length(param1_values), ncol = length(param2_values))

filled.contour(x = param1_values, y = param2_values, z = likelihood_matrix,
               xlab = "param1", ylab = "param2", main = "Llhd Contour Plot")
plot_ly(x = ~param1_values, y = ~param1_values, z = ~likelihood_matrix) |>
  add_surface()

fixed_param1 <- -4.2653076
fixed_param2 <- 1.7911203
beta1_values <- seq(0, 2, by = 0.1)
beta2_values <- seq(-23, -21, by = 0.1)
beta_grid <- expand.grid(beta1 = beta1_values, beta2 = beta2_values)
likelihood_values <- apply(beta_grid, 1, function(betas) {
  lognormal_single(X, Y, theta = c(fixed_param1, fixed_param2, betas['beta1'], betas['beta2']), nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 18)
})
likelihood_values <- -likelihood_values
likelihood_matrix <- matrix(likelihood_values, nrow = length(beta1_values), ncol = length(beta2_values))

filled.contour(x = beta1_values, y = beta2_values, z = likelihood_matrix,
               xlab = "beta1", ylab = "beta2", main = "Log-Normal Frailty Weibull")
plot_ly(x = ~beta1_values, y = ~beta2_values, z = ~likelihood_matrix) |>
  add_surface()


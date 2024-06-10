## Parametric Weibull with Gamma Frailty log-likelihood - Single risk, with Ascertainment correction
######### Gamma Frailty #########
loglik_frailty_single_gamma <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE) {
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  if(base.dist=="lognormal") bparms <- c(theta[1], exp(theta[2]))
  else bparms <- exp(theta[1:(nbase)])
  
  nx <- dim(X)[2]
  xbeta <- c(X %*% theta[(nbase+1):(nbase+nx)])
  k <- exp(theta[length(theta)])
  
  time0 <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status <- Y[,2]
  
  bhaz <- hazards(base.dist, time0, bparms, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time0, bparms, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta) # Cumulative Hazard
  logh <- log(bhaz) + xbeta
  #loglik <- status*logh - H
  loglik <- status*logh 
  h_eachfam <- status * logh
  
  ip <- data$proband == 1
  I_Tp_j.ap_j <- data$I_Tp_j.ap_j[ip]
  df <- subset(data, data$proband == 1)
  d <- aggregate(Y[,2], list(data$famID), sum)[,2]
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  #print(Hfam)
  first_term <- aggregate(h_eachfam, list(data$famID), sum)[,2]
  
  
  ## Terms with Gamma params within the likelihood ##
  term_Gamma <- log(((factorial(k + d - 1))/(factorial(k)*k^(d-1))) * ((1 + Hfam/k)^(-k-d)))
  #term_Gamma <- lfactorial(k + d - 1) - lfactorial(k) - (d - 1) * log(k) - (k + d) * log(1 + Hfam/k)
  #term_Gamma <-  log((gamma(k + d) / (gamma(k + 1) * k^(d - 1))) * ((1 + Hfam / k)^(-k - d)))
  total <- first_term + term_Gamma 
  #print(total)
  total_loglik <- sum(total, na.rm = TRUE)
  
  ## Ascertainment correction ##
  #cagep <- data$currentage[ip]-agemin
  cagep <- data$currentage[ip]
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms, cuts=cuts0)
  H.p <- bcumhaz.p*exp(xbeta.p)
  
  logasc <- I_Tp_j.ap_j * log(1-(1+H.p/k)^(-k)) + (1-I_Tp_j.ap_j) * log((1+H.p/k)^(-k))
  slogasc <- sum(logasc[logasc!=-Inf], na.rm = T)
  
  cloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + total_loglik - slogasc
  
  return(-cloglik)
}

############## Test Zone ################
initial_params <- c(1/41.41327,1,0,0, 1)
X <- as.matrix(brca1_prs[,c("mgeneI", "PRS")], ncol = 2)
Y <- as.matrix(brca1_prs[,c("timeBC", "BC")], ncol = 2)
gamma_forgraph <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
      data = brca1_prs, X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
      agemin = 18, 
      control = list(maxit = 2000), hessian = TRUE)


nomiss_gamma_hessian <- gamma_forgraph$hessian
nomiss_gamma_hessian <- as.matrix(nearPD(nomiss_gamma_hessian)$mat)
vcov_mat <- solve(nomiss_gamma_hessian)
param_est <- gamma_forgraph$par
standard_errors <- sqrt(diag(vcov_mat))
z_scores <- param_est/standard_errors
p_values <- 2*(1-pnorm(abs(z_scores)))
final_results <- list( Estimates = param_est, 
                       Std.Error = standard_errors, 
                       Z = z_scores, 
                       p_val = p_values) # p-value

cox_model <- coxph(Surv(timeBC, BC) ~ mgene + PRS + frailty(famID, distribution = "gamma"), data = brca1_prs)
summary(cox_model)

loglik_frailty_single_gamma(X, Y, theta = c(log(2), log(2), 0, 0), nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 18, k = 1)

## Plot the log-likelihood, fix beta
library(plotly)
param1_values <- seq(from = -5, to = -3, by = 0.01)
param2_values <- seq(from = 0, to = 2, by = 0.01)
param_grid <- expand.grid(param1 = param1_values, param2 = param2_values)
fixed_beta1 <- 1.0635242
fixed_beta2 <- 0.2219278
fixed_k <- 4.0229869
likelihood_values <- apply(param_grid, 1, function(params) {
  loglik_frailty_single_gamma(X, Y, theta = c(params['param1'], params['param2'], fixed_beta1, fixed_beta2, fixed_k), 
                              nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", 
                              frailty.dist = "gamma", agemin = min(brca1_prs_cca$currentage))
})
likelihood_values <- -log(likelihood_values)
likelihood_matrix <- matrix(likelihood_values, nrow = length(param1_values), ncol = length(param2_values))

filled.contour(x = param1_values, y = param2_values, z = likelihood_matrix,
               xlab = "param1", ylab = "param2", main = "Llhd Contour Plot") # Contour plot
plot_ly(x = ~param1_values, y = ~param1_values, z = ~likelihood_matrix) |>
  add_surface() # 3D plot

## Plot the log-likelihood, fix baseline
fixed_param1 <- -4.101604
fixed_param2 <- 1.064875
fixed_k <- 4.354003
beta1_values <- seq(1, 1.4, by = 0.05)
beta2_values <- seq(0, 0.4, by = 0.05)
beta_grid <- expand.grid(beta1 = beta1_values, beta2 = beta2_values)
likelihood_values <- apply(beta_grid, 1, function(betas) {
  loglik_frailty_single_gamma(X, Y, theta = c(fixed_param1, fixed_param2, betas['beta1'], betas['beta2'], fixed_k), 
                              nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", 
                              frailty.dist = "gamma", agemin = 18)
})
likelihood_values <- -likelihood_values
likelihood_matrix <- matrix(likelihood_values, nrow = length(beta1_values), ncol = length(beta2_values))

filled.contour(x = beta1_values, y = beta2_values, z = likelihood_matrix,
               xlab = "beta1", ylab = "beta2", main = "Gamma Frailty Weibull") # Contour plot
plot_ly(x = ~beta1_values, y = ~beta2_values, z = ~likelihood_matrix) |>
  add_surface() # 3D contour plot


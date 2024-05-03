##### Log-Normal Frailty Weibull #####
lognormal_single <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
   
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[length(theta)])
  #sigma <- 1
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)
  logh <- log(bhaz) + xbeta
  loglik <- status1*logh 
  #print(loglik)
  df <- data$df1[ip]
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  hfam <- aggregate(loglik, list(data$famID), sum)[,2]
  #print(hfam)
  #### numerical integral #### 
  logdL <- log(gh(g = Hfam, d = df, p = sigma^2))
  #print(logdL)
  #print(logdL)
  #logdL <- log(gh(g = Hfam, d = df, p = sigma))
  #print(logdL)
  #logdL <- gh(g = Hfam, d = df, p = sigma)
  # Ascertainment correction (design = pop, pop+)
  
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  H.p <- bcumhaz.p*exp(xbeta.p)
  #print(gh(H.p, 0, sigma))
  logasc <- log1p(-laplace(frailty.dist, H.p, sigma))
  #print(logasc)
  #logasc <- 1-laplace(frailty.dist, H.p, sigma)
  #print(logasc)
  #print(sum(logdL[logdL!=-Inf], na.rm = T))
  logasc[logasc == -Inf] <- 0
  slogasc <- sum(logasc[logasc!=-Inf], na.rm=T)
  #print(slogasc)
  #sloglik <- sum(loglik * sum(logdL[logdL!=-Inf], na.rm = T)) 
  #sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + log(sum(integrand[integrand!=-Inf], na.rm = T))
  #- slogasc
  sloglik <- sum(hfam[hfam!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm = T) - slogasc
  #print(theta)
  #return(list(-sloglik, logdL))
  return(-sloglik)
}

library(kinship2)
kin_corr <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID, sex = rep(2, nrow(brca1_prs))))



X <- as.matrix(data.frame(brca1_prs_cca$mgeneI, brca1_prs_cca$PRS), 
               nrow=nrow(brca1_prs_cca), 
               ncol = 2)
Y <- as.matrix(data.frame(brca1_prs_cca$timeBC, brca1_prs_cca$BC), 
               nrow = nrow(brca1_prs_cca), 
               ncol = 2)
gausshermite <- lognormal_single(X = X, Y = Y, theta = c(1/41.41327,1,0), 
                 nbase = 2, data = brca1_prs_cca,
                 design = "pop", base.dist = "Weibull", frailty.dist = "lognormal",
                 agemin = 18)

initial_params <- c(-4.101604,  1.064875,  1.260023,  0.229735,  4.354003)
log_norm_forgraph <- optim(par = initial_params, fn = lognormal_single,
      data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
      agemin = 18, control = list(maxit = 10000))


X1 <- as.matrix(data.frame(brca1_prs_cca$mgeneI), 
                nrow=nrow(brca1_prs_cca), 
                ncol = 1)
initial_params <- c(1/41.41327,1, 0, 1)
optim(par = initial_params, fn = loglik_frailty,
      data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
      agemin = 18, control = list(maxit = 2000))


library(plotly)
param1_values <- seq(from = -7, to = -3, by = 0.1)
param2_values <- seq(from = 3, to = 7, by = 0.1)
param_grid <- expand.grid(param1 = param1_values, param2 = param2_values)
fixed_alpha <- -10.910323
fixed_lambda <- 1.408696
likelihood_values <- apply(param_grid, 1, function(params) {
  lognormal_single(X, Y, theta = c(fixed_alpha, fixed_lambda, params['param1'], params['param2']), nbase = 2, data = brca1_prs_cca1, design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 18)
})
likelihood_values <- -log(likelihood_values)
likelihood_matrix <- matrix(likelihood_values, nrow = length(param1_values), ncol = length(param2_values))

filled.contour(x = param1_values, y = param2_values, z = likelihood_matrix,
               xlab = "param1", ylab = "param2", main = "Llhd Contour Plot")
plot_ly(x = ~param1_values, y = ~param1_values, z = ~likelihood_matrix) |>
  add_surface()

fixed_param1 <- -10.910323
fixed_param2 <- 1.408696
beta1_values <- seq(3, 8, by = 0.01)
beta2_values <- seq(-7, -2, by = 0.01)
beta_grid <- expand.grid(beta1 = beta1_values, beta2 = beta2_values)
#beta_grid <- expand.grid(beta1 = beta1_values)
likelihood_values <- apply(beta_grid, 1, function(betas) {
  lognormal_single(X, Y, theta = c(fixed_param1, fixed_param2, betas['beta1'], betas['beta2']), nbase = 2, data = brca1_prs_cca, design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 18)
})
likelihood_values <- -likelihood_values
likelihood_matrix <- matrix(likelihood_values, nrow = length(beta1_values), ncol = length(beta2_values))



plot(x = beta1_values, y = likelihood_values, type = "l")
lines(x = 1.7)

filled.contour(x = beta1_values, y = beta2_values, z = likelihood_matrix,
               xlab = "beta1", ylab = "beta2", main = "Log-Normal Frailty Contour Plot")
plot_ly(x = ~beta1_values, y = ~beta2_values, z = ~likelihood_matrix) |>
  add_surface()



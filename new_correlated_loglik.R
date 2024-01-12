##### Correlated Frailty Weibull #####
corr_frailty_llhd <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  #
  print(bparms1)
  nx <- dim(X)[2]
  #
  beta <- theta[(nb+1):(nb+nx)]
  #
  print(beta)
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[nb+nx+1])
  #
  print(sigma)
  z_initial <- theta[(nb+nx+2):length(theta)]
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)
  logh <- log(bhaz) + xbeta
  loglik <- status1*logh 
  df1 <- data$df1[ip]
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  
  ##### Sort Hfam to run mgauss.hermite for each family
  
  Sigma <- with(data, kinship2::kinship(id = indID, dadid = fatherID, momid = motherID,
                          sex = rep(2,nrow(data))))
  
  #### numerical integral #### 
  integrand_part <- function(z_star, Hfam1 = Hfam, sigma1 = sigma, Sigma1 = Sigma, status11 = status1) {
    term_1 <- log(z_star) * status11
    term_2 <- z_star * sum(Hfam1[Hfam1 != -Inf], na.rm = TRUE)
    term_3 <- -0.5 * t(z_star) %*% solve(sigma1^2 * Sigma1) %*% z_star
    integral <- sum(term_1) + sum(term_2) + sum(term_3)
  }
  
  hess_matrix <- numDeriv::hessian(integrand_part, z_initial)
  hess_det <- abs(det(hess_matrix))
  
  llhd <- sum(status1 * logh, na.rm = TRUE) + integrand_part(z_star = z_initial) - 0.5*log(hess_det)
  
  
  
  # Ascertainment correction (design = pop, pop+)
  
  #cagep <- data$currentage[ip]-agemin
  #xbeta.p <- xbeta[ip]
  #bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  #H.p <- bcumhaz.p*exp(xbeta.p)
  #logasc <- log(1-laplace(frailty.dist, H.p, sigma))
  #logasc[logasc == -Inf] <- 0
  #slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) 
  #sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm = T) - slogasc
  
  
  return(-llhd)
}

################ Test Zone ################
initial_params <- c(1/41.41327,1,0,0, 1, rep(0, nrow(brca1_prs_cca)))
X <- as.matrix(brca1_prs_cca[,c("mgeneI", "PRS")], ncol = 2)
Y <- as.matrix(brca1_prs_cca[,c("timeBC", "BC")], ncol = 2)
results <- optim(par = initial_params, fn = corr_frailty_llhd,
                 data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
                 design = "pop", frailty.dist = "multivariate normal", base.dist = "Weibull",
                 agemin = 18, 
                 control = list(maxit = 2000))

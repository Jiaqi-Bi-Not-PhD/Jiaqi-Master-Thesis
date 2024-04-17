##########################################
######### Full Log-Likelihood ############
##########################################
log_likelihood <- function(X, Y, z, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  frailty_par <- exp(theta[length(theta)])
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)*z
  logh <- log(bhaz) + xbeta
  loglik <- status1*logh*z 
  df <- data$df1[ip]
  Hfam <- aggregate(H, list(data$famID), sum)[,2] # Summation over individuals within one family
  hfam <- aggregate(loglik, list(data$famID), sum)[,2] # Summation over individuals within one family
  ## Ascertainment correction (design = pop, pop+) ##
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  H.p <- bcumhaz.p*exp(xbeta.p)*z

  logasc[logasc == -Inf] <- 0
  slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) # Family-wise summation of ascertainment term

  sloglik <- sum(hfam[hfam!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm = T) - slogasc

  return(-sloglik)
}

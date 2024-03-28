## Complete likelihood
complete_likelihood <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  if (frailty.dist == "lognormal") sigma <- exp(theta[length(theta)])
  else if (frailty.dist == "gamma") k <- theta[length(theta)]
  #sigma <- exp(theta[length(theta)])
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

  sloglik <- sum(hfam[hfam!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm = T) - slogasc
  return(-sloglik)
}
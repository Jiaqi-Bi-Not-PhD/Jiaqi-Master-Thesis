## Complete likelihood
complete_likelihood <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  if (frailty.dist == "lognormal") sigma <- exp(theta[length(theta)])
  else if (frailty.dist == "gamma") k <- theta[length(theta)]
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  comp_df <- data.frame(famID = integer(), H = numeric(), loglik = numeric())
  
  for (i in seq_along(list_z)) {
    current_famID <- data$famID[i]
    logz <- list_z[[i]]
    H <- bcumhaz*exp(xbeta)*exp(logz)
    logh <- log(bhaz) + xbeta + logz
    loglik <- status1*logh 
    comp_df <- rbind(comp_df, data.frame(famID = rep(current_famID, length(H)), H = H, loglik = loglik))
  }
  Hfam <- aggregate(H ~ famID, data = comp_df, sum)
  hfam <- aggregate(loglik ~ famID, data = comp_df, sum)
  
  #H <- bcumhaz*exp(xbeta)*exp(logz)
  #logh <- log(bhaz) + xbeta + logz
  #loglik <- status1*logh 
  df <- data$df1[ip]
  #Hfam <- aggregate(H, list(data$famID), sum)[,2]
  #hfam <- aggregate(loglik, list(data$famID), sum)[,2]
  
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  for (i in seq_along(list_z)) 
  H.p <- bcumhaz.p*exp(xbeta.p)*exp(logz)
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

lognormal_single <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[length(theta)])
  
  
  time <- Y[,1] - agemin
  cuts0 <- cuts - agemin
  status1 <- Y[,2]
  
  ip <- data$proband == 1
  
  bhaz <- hazards(base.dist, time, bparms1, cuts=cuts0)
  bcumhaz <- cumhaz(base.dist, time, bparms1, cuts=cuts0)
  
  H <- bcumhaz*exp(xbeta)*exp(logz)
  logh <- log(bhaz) + xbeta + logz
  loglik <- status1*logh 

  df <- data$df1[ip]
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  hfam <- aggregate(loglik, list(data$famID), sum)[,2]

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
  
  H.p <- bcumhaz.p*exp(xbeta.p)*exp(logz)
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
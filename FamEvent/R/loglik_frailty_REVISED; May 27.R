loglik_frailty<- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{

if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")

    
if(base.dist=="lognormal") bparms <- c(theta[1], exp(theta[2]))
else bparms <- exp(theta[1:nbase])

nX <- dim(X)[2]
xbeta <- c(X%*%theta[(nbase+1):(nbase+nX)])
kappa <- exp(theta[length(theta)])

time0 <- Y[,1] - agemin
cuts0 <- cuts - agemin
status <- Y[,2]
ip <- data$proband == 1
ip_fam <- aggregate(data$proband, list(unlist(data$famID.byuser)), sum)[,2] # indicates if family has proband or not
wt <- 1
wt_fam <- 1
#wt <- data$weight
#wt_fam <- wt[!duplicated(data$famID.byuser)]
#data <- data |> ### This creates the indicator to test whether proband is affected ###
#  dplyr::mutate(I_Tp_j.ap_j = ifelse(proband == 1 & time < currentage, 1, 0)) ### This creates the indicator to test whether proband is affected ###
i_ap <- with(data, ifelse(proband == 1 & time < currentage, 1, 0))[ip] ### indicates if proband is affected ###

bhaz <- hazards(base.dist, time0, bparms, cuts=cuts0)
bcumhaz <- cumhaz(base.dist, time0, bparms, cuts=cuts0)

H <- bcumhaz*exp(xbeta)
logh <- log(bhaz) + xbeta
loglik <-  wt * (status*logh )

df <- data$df[!duplicated(data$famID.byuser)]
s <- aggregate(H, list(unlist(data$famID.byuser)), sum)[,2]
logdL <- wt_fam*log( dlaplace(frailty.dist, g=s, d=df, k=kappa) )


# Ascertainment correction (design = pop, pop+)
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms, cuts=cuts0)
  laplace.p <- laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa)
  logasc.p <- ifelse(i_ap==1, log(1-laplace.p), log(laplace.p))
  logasc <- wt_fam*ifelse(ip_fam==0, 0, logasc.p)
  ip_fam[ip_fam!=0] <- i_ap
  #logasc <- wt.p*log(1-laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa))*I_ap + 
  #  wt.p*log(laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa))*(1-I_ap) ### Not all probands are affected ###
  
  
  logasc[logasc == -Inf] <- 0
  sloglik <- sum(loglik[loglik!=-Inf & loglik != Inf], na.rm=T) + sum(logdL[logdL!=-Inf & logdL != Inf], na.rm=T) - sum(logasc[logasc!=-Inf & logasc != Inf], na.rm=T)
  loglik[ip_fam] <- loglik[ip_fam] + logdL[ip_fam==1] - logasc[ip_fam==1]
  
  #print(c(theta, -sloglik))
  #print(c(theta, -sloglik, sum(logdL), sum(logasc)))
  if(vec) return(-loglik)
  else  return(-sloglik)
  
}


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

wt <- 1
wt.p <- 1

bhaz <- hazards(base.dist, time0, bparms, cuts=cuts0)
bcumhaz <- cumhaz(base.dist, time0, bparms, cuts=cuts0)

H <- bcumhaz*exp(xbeta)
logh <- log(bhaz) + xbeta
loglik <-  wt * (status*logh )

df <- data$df1[ip]
s <- aggregate(H, list(data$famID), sum)[,2]
logdL <- wt.p*log( dlaplace(frailty.dist, g=s, d=df, k=kappa) )

# Ascertainment correction (design = pop, pop+)
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms, cuts=cuts0)
  logasc <- wt.p*log(1-laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa))
  logasc[logasc == -Inf] <- 0
  sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL, na.rm=T) - sum(logasc, na.rm=T)
  loglik[ip] <- loglik[ip] + logdL - logasc
  
  #print(c(theta, -sloglik))
  
  if(vec) return(-loglik)
  else  return(-sloglik)
  
}

X <- as.matrix(data.frame(brca1_prs$mgeneI, brca1_prs$PRS), 
               nrow=nrow(brca1_prs), 
               ncol = 2)
Y <- as.matrix(data.frame(brca1_prs$timeBC, brca1_prs$BC), 
               nrow = nrow(brca1_prs), 
               ncol = 2)
initial_params <- c(1/41.41327,1, 0, 0, 1)
log_norm_forgraph <- optim(par = initial_params, fn = loglik_frailty,
                           data = brca1_prs, X = X, Y = Y, nbase = 2,
                           design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
                           agemin = 18, control = list(maxit = 2000))

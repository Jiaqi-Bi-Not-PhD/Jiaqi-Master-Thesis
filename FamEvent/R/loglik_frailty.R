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
data <- data |> ### This creates the indicator to test whether proband is affected ###
  dplyr::mutate(I_Tp_j.ap_j = ifelse(proband == 1 & time < currentage, 1, 0)) ### This creates the indicator to test whether proband is affected ###
I_Tp_j.ap_j <- data$I_Tp_j.ap_j[ip] ### This creates the indicator to test whether proband is affected ###


bhaz <- hazards(base.dist, time0, bparms, cuts=cuts0)
bcumhaz <- cumhaz(base.dist, time0, bparms, cuts=cuts0)

H <- bcumhaz*exp(xbeta)
logh <- log(bhaz) + xbeta
loglik <-  wt * (status*logh )

data_df <- data |> dplyr::group_by(famID) |> dplyr::summarise(df = sum(status))
#data <- merge(data, data_df, by = "famID")
df <- data_df$df
s <- aggregate(H, list(data$famID), sum)[,2]
logdL <- wt.p*log( dlaplace(frailty.dist, g=s, d=df, k=kappa) )


# Ascertainment correction (design = pop, pop+)
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms, cuts=cuts0)
  logasc <- wt.p*log(1-laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa))*I_Tp_j.ap_j + wt.p*log(laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa))*(1-I_Tp_j.ap_j) ### Not all probands are affected ###
  logasc[logasc == -Inf] <- 0
  sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm=T) - sum(logasc[logasc!=-Inf & logasc != Inf], na.rm=T)
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
                           design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                           agemin = 18, control = list(maxit = 2000))

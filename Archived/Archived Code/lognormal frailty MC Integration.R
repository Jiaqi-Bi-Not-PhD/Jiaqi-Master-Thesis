###### Log-Normal MC Integration #####
##### Log-Normal Frailty Weibull #####
lognormal_single_MC <- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE, nmc = 1000)
{
  
  if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")
  
  nb <- nbase # number of baselines parameters
  
  if(base.dist == "lognormal") bparms1 <- c(theta[1], exp(theta[2]))
  else bparms1 <- exp(theta[1:nbase])
  
  nx <- dim(X)[2]
  xbeta <- c(X%*%theta[(nb+1):(nb+nx)])
  sigma <- exp(theta[length(theta)])
  #print(sigma)
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
  df <- data$df1[ip]
  #print(length(df))
  Hfam <- aggregate(H, list(data$famID), sum)[,2]
  #print(length(Hfam))
  #### numerical integral using MC Integration #### 
  sample_z <- rlnorm(n = nmc, meanlog = 0, sdlog = sigma)
  integrand <- mapply(function(Hfamj, dfj) {
    exp(-sample_z * Hfamj) * sample_z^dfj
  }, Hfamj = Hfam, dfj = df)
 # print(dim(integrand))
  #print(c(Hfam, df))
  #integrand <- colMeans(integrand[integrand!=-Inf], na.rm = T)
  integrand_mean <- colMeans(integrand)# compute the MC means for each family
  #print(integrand_mean)
  integrand_mean <- log(integrand_mean)
  #print(integrand_mean)
  #logdL <- log(gh(g = Hfam, d = df, p = sigma))
  # Ascertainment correction (design = pop, pop+)
  
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  H.p <- bcumhaz.p*exp(xbeta.p)
  ascintegrand <- sapply(H.p, function(H.pj) {
    exp(-sample_z * H.pj)
  })
  #print(ascintegrand)
  ascintegrand <- colMeans(ascintegrand)
  #print(length(ascintegrand))
  #print(ascintegrand)
  logasc <- log(1-ascintegrand)
  #print(logasc)
  #logasc <- log(1-laplace(frailty.dist, H.p, sigma))
  #print(log(sum(integrand_mean[integrand_mean!=-Inf], na.rm = T)))
  #print(sum(integrand_mean[integrand_mean!=-Inf], na.rm = T))
  
  logasc[logasc == -Inf] <- 0
  slogasc <- sum(logasc[logasc!=-Inf], na.rm=T) 
  #print(slogasc)
  #sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(integrand_mean[integrand_mean!=-Inf], na.rm = T) - slogasc
  sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(integrand_mean[integrand_mean!=-Inf], na.rm = T) - slogasc
  
  #return(list(-sloglik, integrand_mean))
  return(-sloglik)
}

######### Testing Zone #########
X <- as.matrix(data.frame(brca1_prs_cca$mgeneI, brca1_prs_cca$PRS), 
               nrow=nrow(brca1_prs_cca), 
               ncol = 2)
Y <- as.matrix(data.frame(brca1_prs_cca$timeBC, brca1_prs_cca$BC), 
               nrow = nrow(brca1_prs_cca), 
               ncol = 2)
set.seed(123)
MCint <- lognormal_single_MC(X = X1, Y = Y, theta = c(1/41.41327,1,0, 2), 
                 nbase = 2, data = brca1_prs_cca,
                 design = "pop", base.dist = "Weibull", frailty.dist = "lognormal",
                 agemin = 18, nmc = 10000)
initial_params <- c(1/41.41327,1,0)

initial_params <- c(1, 1, 0, 0, 1)
MCinte_lognormal_results <- optim(par = initial_params, fn = lognormal_single_MC,
                data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
                design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
                agemin = 18,  
                control = list(maxit = 2000), nmc = 5000)

cbind(MCint[[2]], gausshermite[[2]])

initial_params <- c(1/41.41327,1,0,2)
MCinte_lognormal_results1 <- optim(par = initial_params, 
                                   fn = lognormal_single_MC,
                                  data = brca1_prs_cca, X = X1, Y = Y, 
                                  nbase = 2,
                                  design = "pop", frailty.dist = "lognormal", 
                                  base.dist = "Weibull",
                                  agemin = 18,  
                                  control = list(maxit = 2000))



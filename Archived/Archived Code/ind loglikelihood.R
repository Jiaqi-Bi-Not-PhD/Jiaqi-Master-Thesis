ind_loglikelihood <- function(X, Y, z, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE) {
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
  H <- exp(-bcumhaz*exp(xbeta)*z)
  logH <- log(H)
  h <- bhaz*exp(xbeta)*z
  term1 <- h^status1
  # df <- data$df1[ip] # Events within family

  ## Ascertainment correction (design = pop, pop+) ##
  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms1, cuts=cuts0)
  
  H.p <- bcumhaz.p*exp(xbeta.p)*z
  
  asc <- H.p
  
  llhd_ind <- log(term1) + logH - asc
  return(llhd_ind)
}

################ Test Zone ##################
X <- as.matrix(data.frame(brca1_prs_cca$mgeneI, brca1_prs_cca$PRS), 
               nrow=nrow(brca1_prs_cca), 
               ncol = 2)
Y <- as.matrix(data.frame(brca1_prs_cca$timeBC, brca1_prs_cca$BC), 
               nrow = nrow(brca1_prs_cca), 
               ncol = 2)
initial_params <- c(-19.0481278, 1.1189520, 0.4178628, 3.9732051, 2.8852498)
ind_loglikelihood(X = X, Y = Y, z = 1, theta = initial_params, nbase = 2, data = brca1_prs_cca, 
                  design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 18)

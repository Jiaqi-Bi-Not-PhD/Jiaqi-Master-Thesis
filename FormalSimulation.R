cumhaz  <-  function(dist="Weibull", t, parms, cuts=NULL){
  
  if(dist=="Weibull")	 	chaz <- (parms[1]*t)^parms[2] 
  else if(dist=="loglogistic")	chaz <- log(1+(parms[1]*t)^parms[2])
  else if(dist=="Gompertz") chaz <- parms[1]*(exp(parms[2]*t)-1)/parms[2]
  else if(dist=="lognormal") chaz <- -log(1-pnorm((log(t)-parms[1])/parms[2]))
  else if(dist=="gamma") chaz <- -log(1-pgamma(t,shape=parms[2], scale=1/parms[1] ))
  else if(dist=="logBurr") chaz <- parms[3]*log(1+(parms[1]*t)^parms[2]/parms[3])
  else if(dist=="piecewise") chaz <- Hpch(t, cuts=cuts, levels=parms)
  else stop("Unrecognized baseline distribution")
  chaz
}

hazards  <-  function(dist="Weibull", t, parms, cuts=NULL){
  
  if(dist=="Weibull")	 	haz <- (parms[1]^parms[2])*parms[2]*t^(parms[2]-1) 
  else if(dist=="loglogistic")	haz <- (parms[1]^parms[2])*parms[2]*t^(parms[2]-1)/(1+(parms[1]*t)^parms[2])
  else if(dist=="Gompertz") haz <- parms[1]*exp(parms[2]*t)
  else if(dist=="lognormal") haz <- dnorm((log(t)-parms[1])/parms[2])/(parms[2]*t*(1-pnorm((log(t)-parms[1])/parms[2])))
  else if(dist=="gamma") haz <- dgamma(t,shape=parms[2], scale=1/parms[1])/(1-pgamma(t,shape=parms[2], scale=1/parms[1]))
  else if(dist=="logBurr") haz <- (parms[1]^parms[2])*parms[2]*parms[3]*t^(parms[2]-1)/(parms[3]+(parms[1]*t)^parms[2])
  else if(dist=="piecewise") haz <- hpch(t, cuts=cuts, levels=parms)
  else stop("Unrecognized baseline distribution")
  haz
}

gh <- function(g,d,p){
  xk<- c(-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212, -1.2340762154,-0.737473728545, -0.245340708301,0.245340708301, 0.737473728545, 1.2340762154,1.73853771212,2.25497400209, 2.78880605843,  3.34785456738, 3.94476404012,4.60368244955, 5.38748089001)
  wk <- c(2.22939364553415129252E-13, 4.3993409922731805536E-10, 1.086069370769281694E-7, 7.80255647853206369415E-6, 2.28338636016353967257E-4, 0.00324377334223786183218, 0.0248105208874636108822, 0.1090172060200233200138, 0.2866755053628341297197, 0.46224366960061008965, 0.46224366960061008965, 0.28667550536283412972, 0.1090172060200233200138, 0.0248105208874636108822, 0.00324377334223786183218, 2.283386360163539672571E-4, 7.8025564785320636941E-6, 1.086069370769281694E-7, 4.39934099227318055363E-10, 2.22939364553415129252E-13)
  yk <- exp(sqrt(2/p)*xk) # p = 1/var = kappa  (previously p = var = 1/kappa)
  apply(cbind(d, g), 1, function(x) sum(wk*yk^x[1]*exp(-x[2]*yk))/sqrt(pi) )
  
  #wt <- c(0.898591961453,0.704332961176, 0.62227869619, 0.575262442852,0.544851742366, 0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667, 0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949,0.544851742366, 0.575262442852, 0.62227869619, 0.704332961176,0.898591961453)
  #apply(wt*exp(-t(as.matrix(g)%*%exp(as.matrix(t(sqrt(p)*xk)))) +d*sqrt(p)*xk - xk^2/2),2, sum)
  
}

laplace <- function(dist, g, k){
  if(dist=="gamma")  (1+g/k)^(-k)
  else if(dist=="positives") exp(-g^k)
  else if(dist=="powerseries") 1-(1-exp(-g))^k
  else if(dist=="logarithmic") -log(exp(-g)*(exp(-k)-1)+1)/k
  else if(dist=="lognormal") gh(g, 0, k)
  else stop("Unrecognized frailty distribution")
}

dlaplace <- function(dist, g, d, k){
  if(dist=="gamma")  factorial(k+d-1)*(1+g/k)^(-k-d)/factorial(k)/k^(d-1)
  #else if(dist=="positives") DD(expression(exp(-g^k)), name="g", order=d) # 0<k<=1
  #else if(dist=="power") DD(expression(1-(1-exp(-g))^k), name="g", order=d)
  #else if(dist=="logarithmic") DD(expression(-log(exp(-g)*(exp(-k)-1)+1)/k), name="g", order=d ) # 0<k<=1
  else if(dist=="lognormal") gh(g,d, k)
  else stop("Unrecognized frailty distribution")
}

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

# frailty.dist option 
penmodel <- function(formula, cluster="famID", gvar="mgene", parms, cuts=NULL, data, design="pop", base.dist="Weibull", frailty.dist="none", agemin=NULL, robust=FALSE){
  
  if(!frailty.dist%in%c("none", "gamma", "lognormal")) stop("frailty.dist should be one of \"none\", \"gamma\", or \"lognormal\".")
  nfp <- ifelse(frailty.dist=="none", 0, 1)
  
  if(any(is.na(data[, gvar]))) stop("data include missing genetic information, use penmodelEM function.")
  options(na.action='na.omit')
  
  if(is.null(agemin)){
    agemin <- attr(data, "agemin")
    if(is.null(agemin)) {
      agemin <- 0
      warning("agemin = 0 was used or assign agemin to attr(data, \"agemin\").")
    }
  }
  
  
  if(sum(data$time <=  agemin, na.rm = TRUE) > 0) cat("Individuals with time < agemin (", agemin,") were removed from the analysis.\n")
  
  data <- data[data$time >=  agemin, ]
  data$famID.byuser <- data[, cluster]
  
  x.names <- attr(terms(formula), "term.labels")
  i.missing <- apply(data[, x.names], 1, anyNA) # indicates if any x variable is missing
  newdata <- data[!i.missing,]
  
  
  m <- model.frame(formula, data) # missing data were removed
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  
  if (!inherits(Y, "Surv")) stop("Response must be a survival object.")
  type <- attr(Y, "type")
  if (type == "counting") stop("start-stop type Surv objects are not supported.")
  if (type == "mright" || type == "mcounting") stop("multi-state survival is not supported.")
  
  if(base.dist=="piecewise"){
    if(is.null(cuts)) stop("The cuts should be specified")
    if(any(cuts > max(Y[,1]) | cuts < min(Y[,1]))) stop("Some value(s) of the cuts are beyond the range.")
  }
  
  
  X <- model.matrix(Terms, m)
  
  nvar <- ncol(X)-1
  var.names <- colnames(X)[-1]
  if(nvar==1) X <- matrix(X[,-1])
  else X <- X[,-1]
  
  #number of parameters for baseline
  nbase <- ifelse(base.dist=="logBurr", 3, ifelse(base.dist=="piecewise", length(cuts)+1, 2)) 
  
  #nk <- ifelse(is.null(frailty.dist),0,1)
  colnames(X) <- var.names
  vbeta <- parms[(nbase+1):(nbase+nvar)]
  kappa <- parms[(nbase+nvar+1)]
  #  if(length(vbeta) != nvar) stop("The size of parms is incorrect.")
  if(length(parms) != (nvar+nbase+nfp) ) stop("The size of parms is incorrect.")
  
  if(is.null(data$weight)) data$weight <- 1
  
  if(base.dist=="lognormal"){
    if(parms[2] <= 0) stop("parms[2] has to be > 0")
    else parms[1] <- exp(parms[1])
  } 
  else if(any(parms[1:nbase]<=0)) stop("All baseline parameters should be > 0")
  
  if(frailty.dist=="none"){
    est1 <- optim(c(log(parms[1:nbase]), vbeta), loglik_ind, X=X, Y=Y, cuts=cuts, 
                  nbase=nbase, data=newdata, design=design, base.dist=base.dist, 
                  agemin=agemin, control = list(maxit = 50000), hessian=TRUE)
  }
  else{ # frailty model
    if(suppressWarnings(is.null(data$df))){
      df <- aggregate(Y[,2], list(newdata$famID), sum)[,2]
      fsize <- aggregate(Y[,2], list(newdata$famID), length)[,2]
      newdata$df <- rep(df, fsize)
    }
    est1 <- optim(c(log(parms[1:nbase]), vbeta, log(kappa)), loglik_frailty, X=X, Y=Y, 
                  cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, 
                  frailty.dist=frailty.dist, agemin=agemin, vec=FALSE,
                  control = list(maxit = 50000), hessian=TRUE)
  }
  
  logLik <- -est1$value
  EST <- est1$par
  H <- est1$hessian
  Var <- try(solve(H), TRUE)
  
  
  if(!is.null(attr(Var,"class"))) stop("Model did not converge.\n  Try again with different initial values.")
  else{ 
    bparms.name <- c("log.lambda","log.rho", "log.eta")
    if(base.dist=="lognormal") bparms.name[1] <- "lambda" 
    else if(base.dist=="piecewise") bparms.name <- paste0("log.q", 1:nbase,")")
    
    parms.cov <- Var
    parms.se <- sqrt(diag(parms.cov))
    parms.cov.robust <- parms.se.robust <- NULL
    if(frailty.dist=="none") parms.name <- c(bparms.name[1:nbase], colnames(X))
    else parms.name <- c(bparms.name[1:nbase], colnames(X), "log.kappa")
    
    names(EST) <- names(parms.se)  <- rownames(parms.cov) <- colnames(parms.cov) <- parms.name
    
    if(robust){
      if(frailty.dist=="none")  grad <- jacobian(loglik_ind, est1$par, X=X, Y=Y, cuts=cuts, nbase=nbase, data=data, design=design, base.dist=base.dist, agemin=agemin, vec=TRUE)
      else grad <- jacobian(loglik_frailty, est1$par, X=X, Y=Y, cuts=cuts, nbase=nbase, data=data, design=design, base.dist=base.dist, frailty.dist=frailty.dist, agemin=agemin, vec=TRUE)
      Jscore <- t(grad)%*%grad
      parms.cov.robust <- Var%*%(Jscore)%*%Var
      parms.se.robust <- sqrt(diag(parms.cov.robust))
      rownames(parms.cov.robust) <- colnames(parms.cov.robust) <- parms.name
    }
  }
  
  aic = 2*length(EST) - 2*logLik
  
  out <- list(estimates = EST, varcov = parms.cov, varcov.robust = parms.cov.robust, se = parms.se, se.robust = parms.se.robust,logLik = logLik, AIC=aic)  
  class(out) <- "penmodel"
  attr(out, "design") <- design
  attr(out, "base.dist") <- base.dist
  attr(out, "frailty.dist") <- frailty.dist
  attr(out, "agemin") <- agemin
  attr(out, "cuts") <- cuts
  attr(out, "nbase") <- nbase
  attr(out, "data") <- data
  attr(out, "robust") <- robust
  attr(out, "formula") <- formula
  attr(out, "X") <- X
  attr(out, "Y") <- Y
  invisible(out)
  
}#end




#################################################
################# Formal Simulation #############
################# 100 Simulations ###############
#################################################
## True thetas
true_beta1 <- 1
true_beta2 <- 3
true_beta3 <- 3
true_kappa <- log(2)
true_alpha <- log(0.035)
true_lambda <- log(2.3)
  
n_simulations <- 100
n_cores <- detectCores() - 1 

#start_time <- Sys.time() # Starting time 
#######################################################################################
################################ Generate 100 datasets ################################
#######################################################################################
## 500, 200, 50 families => Gamma or LogNormal
Nfam_values <- c(500, 200, 50)
frailty_dist_values <- c("gamma", "lognormal")

scenarios <- expand.grid(Nfam = Nfam_values, frailty_dist = frailty_dist_values)
scenarios_list <- split(scenarios, seq(nrow(scenarios)))

# Function to run a single simulation and analysis
run_simulation <- function(params, sim_index) {
  #library(survival)
  #library(survPen)
  
  Nfam <- params$Nfam
  frailty_dist <- params$frailty_dist
  
  repeat {
    ## Data Generation
    famx <- simfam(N.fam = Nfam, design = "pop", variation = "frailty", 
                   base.dist = "Weibull", frailty.dist = frailty_dist, interaction = FALSE,
                   add.x = TRUE, x.dist = "normal", x.parms = c(0, 1), depend = 2, 
                   base.parms = c(0.035,2.3), vbeta = c(1, 3, 3)) 
    
    ## If warning exists when the generated data runs the analysis => re-generate a new data and abandon the current
    result <- tryCatch({
      test_model <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                             design = "pop", base.dist = "Weibull", frailty.dist = frailty_dist, 
                             agemin = min(famx$currentage[famx$status == 1]), data = famx,
                             parms = c(1/41.41327,1,0,0,0, 1))
      list(model = test_model, warning = NULL)
    }, warning = function(w) {
      list(model = NULL, warning = w)
    }, error = function(e) {
      list(model = NULL, error = e)
    })
    
    if (is.null(result$warning) && is.null(result$error)) {
      return(famx) ## Only store the dataset
    } else {
      cat("Warning or error occurred in iteration", sim_index, "- rerunning\n")
    }
  }
}

## Generate 100 complete datasets for 2*3 = 6 types of data
n_simulations <- 100
n_cores <- detectCores() - 1  # 7 cores on my own computer

# Create a cluster
cl <- makeCluster(n_cores)
#clusterEvalQ(cl, {
#  library(survival)
#  library(survPen)
#})

clusterExport(cl, c("simfam", "penmodel", "run_simulation", "gh",
                    "loglik_frailty", "cumhaz", "hazards", "Surv", "dlaplace",
                    "laplace", "scenarios_list", "n_simulations"))  

## Apply the function in parallel for all scenarios
results_parallel <- parLapply(cl, seq_along(scenarios_list), function(scenario_idx) {
  scenario_params <- scenarios_list[[scenario_idx]]
  lapply(1:n_simulations, function(sim_index) {
    run_simulation(scenario_params, sim_index)
  })
})

stopCluster(cl)

for (i in seq_along(scenarios_list)) {
  scenario_params <- scenarios_list[[i]]
  scenario_name <- paste("simulated_dataset_100_Complete", scenario_params$frailty_dist, scenario_params$Nfam, "fam", sep = "_")
  
  assign(scenario_name, results_parallel[[i]])
}

# simulated_dataset_100_Complete_lognormal_50_fam
# simulated_dataset_100_Complete_lognormal_200_fam
# simulated_dataset_100_Complete_lognormal_500_fam
# simulated_dataset_100_Complete_gamma_50_fam
# simulated_dataset_100_Complete_gamma_200_fam
# simulated_dataset_100_Complete_gamma_500_fam
nomiss_datasets <- list(simulated_dataset_100_Complete_lognormal_50_fam,
                        simulated_dataset_100_Complete_lognormal_200_fam,
                        simulated_dataset_100_Complete_lognormal_500_fam,
                        simulated_dataset_100_Complete_gamma_50_fam,
                        simulated_dataset_100_Complete_gamma_200_fam,
                        simulated_dataset_100_Complete_gamma_500_fam)

#######################################################################################
################################ Generate 100 datasets ends ###########################
#######################################################################################

#######################################################################################
################################ Generate MAR #########################################
#######################################################################################
proportions <- c(0.20, 0.40, 0.60, 0.80)

miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)

## Missingness & proportion function
introduce_missingness <- function(index, dataset_list, prop, miss_pattern) {
  famx <- dataset_list[[index]]
  
  ampute_test <- ampute(data = famx, prop = prop, patterns = miss_pattern)
  reasonable_weights <- ampute_test$weights
  reasonable_weights[1,] <- c(0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0)
  ampute_test <- ampute(data = famx, prop = prop, patterns = miss_pattern, 
                        mech = "MAR", weights = reasonable_weights)
  miss_famx <- ampute_test$amp
  
  return(miss_famx)
}

## Function to run the missingness introduction for all proportions
run_for_all_proportions <- function(i, dataset_list, proportions, miss_pattern) {
  results <- lapply(proportions, function(prop) {
    introduce_missingness(i, dataset_list, prop, miss_pattern)
  })
  names(results) <- as.character(proportions)
  return(results)
}

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterEvalQ(cl, library(mice))  

## List of complete datasets
complete_datasets <- list(
  simulated_dataset_100_Complete_gamma_50_fam = simulated_dataset_100_Complete_gamma_50_fam,
  simulated_dataset_100_Complete_gamma_200_fam = simulated_dataset_100_Complete_gamma_200_fam,
  simulated_dataset_100_Complete_gamma_500_fam = simulated_dataset_100_Complete_gamma_500_fam,
  simulated_dataset_100_Complete_lognormal_50_fam = simulated_dataset_100_Complete_lognormal_50_fam,
  simulated_dataset_100_Complete_lognormal_200_fam = simulated_dataset_100_Complete_lognormal_200_fam,
  simulated_dataset_100_Complete_lognormal_500_fam = simulated_dataset_100_Complete_lognormal_500_fam
)

clusterExport(cl, c("introduce_missingness", "run_for_all_proportions", "miss_pattern", "proportions",
                    "complete_datasets"))

mar_datasets <- list()

## Iterate over each complete dataset and apply the missing data mechanism
for (dataset_name in names(complete_datasets)) {
  dataset <- complete_datasets[[dataset_name]]
  clusterExport(cl, "dataset")
  ## Apply the function in parallel
  results_parallel <- parLapply(cl, 1:100, function(i) {
    run_for_all_proportions(i, dataset, proportions, miss_pattern)
  })
  
  ## Store the results in the mar_datasets list
  for (prop in proportions) {
    mar_datasets[[paste0("simulated_dataset_100_MAR", prop * 100, "_", dataset_name)]] <- lapply(1:100, function(i) results_parallel[[i]][[as.character(prop)]])
  }
}

stopCluster(cl)

### mar_datasets contains 24 lists of 100 datasets -> 2400 datasets
### nomiss_datasets contains 6 lists of 100 datasets -> 600 datasets
### mar_datasets
names(mar_datasets)
#saveRDS(complete_datasets, file = "100 complete datasets lists.RData")
#saveRDS(mar_datasets, file = "100 mar datasets lists.RData")
#test_list <- readRDS("100 mar datasets lists.RData")


#######################################################################################
################################ Generate MAR Ends ####################################
#######################################################################################





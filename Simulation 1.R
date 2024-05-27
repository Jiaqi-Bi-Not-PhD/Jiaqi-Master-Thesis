######################################################
################## Simulation 1 ######################
######################################################
####### Missingness generated from mice package ######
######################################################
brca1_prs <- brca1_prs |>
  mutate(weight = 1)
brca1_prs_cca <- brca1_prs |>
  filter(!is.na(PRS))
## Generate kinship data - Number of family = 10
install.packages("FamEvent_3.1.tar.gz", type = "source")
library(FamEvent)
library(mice)

set.seed(123123)
famx <- simfam(N.fam = 200, design = "pop", variation = "frailty", 
               base.dist = "Weibull", frailty.dist = "gamma", interaction = FALSE,
               add.x = TRUE, x.dist = "normal", x.parms = c(0, 1),  depend = 2, 
               base.parms = c(0.016,3), vbeta = c(1, 3, 3)) # Number of family = 10

## Test if Gamma frailty model gives the correct estimates ##
test_model <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(test_model)
## Test if Gamma frailty model gives the correct estimates ##

##############################################################################################
IBDmatrix <- diag(1, dim(famx)[1])

inputdata <- famx[, !names(famx) %in%c("time", "status", "ageonset", "fsize", "naff")]

fam2 <- simfam2(design = "pop", variation = c("kinship", "IBD"),
                depend = c(2, 2), base.dist="Weibull", base.parms =c(0.016, 3),
                var_names = c("gender", "mgene", "newx"), vbeta = c(1, 3, 3),
                agemin=20, inputdata=inputdata, IBD=IBDmatrix) # covariates are gender, mgene, newx

test_model2 <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = fam2,
                       parms = c(1/41.41327,1,0,0,0, 1))
summary(test_model2)

fam2 <- fam2 |>
  mutate(gender2 = ifelse(gender == 1, "male", "female"))
kinship_mat <- with(fam2, kinship(id = indID, dadid = fatherID, momid = motherID, sex = gender))
kinship_mat <- Matrix(kinship_mat, sparse = TRUE)

test_model2 <- coxme(Surv(time, status) ~ gender + mgene + newx + (1|indID), data = fam2, varlist = list(kinship_mat))
summary(test_model2)


test_model <- penmodel(Surv(time, status) ~ mgene + newx + gender, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1))
summary(test_model)

fam2 |>
  ggplot(aes(x = newx)) + 
  geom_density() # Check the normality

fam2 <- fam2 |>
  mutate(I_Tp_j.ap_j = ifelse(proband == 1 & time < currentage, 1, 0))
###############################################################################################

## Estimates when there are no missing data
nomiss_gamma <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(nomiss_gamma)

gammafr <- as.vector(summary(nomiss_gamma)$estimates[1:2,1])
logalpha <- gammafr[1]
loglambda <- gammafr[2]

famx <- famx |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)))

nomiss_lognorm <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = famx,
                         parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(nomiss_lognorm)


################# Suppose 50% missing - Missing at Random ####################
################# Suppose a reasonable weights ###############################
miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)
ampute_test <- ampute(data = famx, prop = 0.2, patterns = miss_pattern)
md.pattern(ampute_test$amp, rotate.names = TRUE) # Missing pattern
ampute_test$mech
ampute_test$weights
reasonable_weights <- ampute_test$weights
reasonable_weights[1,] <- c(0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0)
ampute_test <- ampute(data = famx, prop = 0.2, patterns = miss_pattern, 
                      mech = "MAR", weights = reasonable_weights) # Missing proportions
miss50_famx <- ampute_test$amp
## CCA Gamma
miss50_famx_CCA <- miss50_famx |>
  filter(!is.na(newx)) 
miss50_gamma <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx_CCA,
                         parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(miss50_gamma)
## CCA Log-Normal
miss50_lognorm <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = miss50_famx_CCA,
                         parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(miss50_lognorm)

########### Multiple Imputation using kinship and H0 #############
famx_newx_imp_fam <- list() 
## Step 1 - Kinship matrix
kinship_mat <- with(miss50_famx, kinship(id = indID, dadid = fatherID, momid = motherID,
                                         sex = gender))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(miss50_famx)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step X - Initial baseline cumulative hazard
miss50_gamma_cca <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx_CCA,
                             parms = c(1/41.41327,1,0,0,0, 1))
baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
miss50_famx <- miss50_famx |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)))

start_time <- Sys.time() # Starting time 
for (m in 1:20) {
  ## Step 2 - empirical estimates
  model_test <- relmatLmer(newx ~ H0 * (status + proband) + (1|indID), data = miss50_famx, relmat = list(indID = kinship_mat))
  summary(model_test)
  X <- model.matrix(~ H0 + status + proband + H0:status + H0:proband, data = miss50_famx) # imputation model design matrix
  V <- vcov(model_test)
  
  #betas <- coef(model_test)
  betas <- as.vector(summary(model_test)$coefficients[,1]) # beta coefficients
  sigma_g_2 <- attr(VarCorr(model_test)$indID, "stddev")^2 # genetic variance
  sigma_e_2 <- attr(VarCorr(model_test), "sc")^2 # residual variance
  
  Sigma <- sigma_g_2*kinship_mat_sparse + sigma_e_2*Iden_mat_sparse # Sparse matrix for Sigma
  Sigma_mat <- as.matrix(Sigma) # Non-sparse
  
  ## Step 3 - conditional variance
  num_cores <- detectCores() - 2 # 6 cores
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist = c("Sigma")) 
  clusterEvalQ(cl, library(Matrix))
  
  cond_var <- function(i) {
    conditional_var <- Sigma[i,i] - Sigma[i,-i] %*% solve(Sigma[-i,-i]) %*% Sigma[-i,i]
    return(conditional_var)
  }
  conditional_variances <- parSapply(cl, 1:nrow(miss50_famx), cond_var)
  conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))
  
  stopCluster(cl)
  
  ## Imputation steps
  miss50_famx <- miss50_famx |>
    mutate(cond_var = conditional_variances_temp)
  
  newx <- miss50_famx$newx
  

  ## Step 4
  w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  
  ## Step 5
  betastar <- betas + w_1 %*% chol(V)
  
  ## Step 6
  mu_star <- X %*% as.vector(betastar)
  
  ## Step 7
  num_cores <- detectCores() - 2 # 6 cores
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist = c("Sigma", "newx", "mu_star")) 
  clusterEvalQ(cl, library(Matrix))
  
  E_cond <- function(i) {
    y_minus_i <- newx[-i]
    mu_star_minus_i <- mu_star[-i]
    non_NA <- which(!is.na(y_minus_i))
    y_minus_i <- y_minus_i[non_NA]
    mu_star_minus_i <- mu_star_minus_i[non_NA]
    Sigma_nonNA <- Sigma[-i,-i][non_NA, non_NA]
    Sigmahat <- Sigma[i,-i][non_NA]
    
    conditional_Expect <- mu_star[i] + Sigmahat %*% solve(Sigma_nonNA) %*% (y_minus_i - mu_star_minus_i)
    return(conditional_Expect)
  }
  conditional_expectations <- parSapply(cl, 1:nrow(miss50_famx), E_cond)
  conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
  
  conditional_vars <- miss50_famx$cond_var
  stopCluster(cl)
  
  ## Step 8
  w2i <- rnorm(n = nrow(miss50_famx), mean = 0, sd = 1)
  
  ## Step 9
  newx_star <- conditional_expectations + w2i * sqrt(conditional_vars)
  miss50_famx <- miss50_famx |>
    mutate(newx_I = newx_star) |>
    mutate(newx_I = ifelse(!is.na(newx), newx, newx_I)) |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)) )
  famx_newx_imp_fam[[m]] <- miss50_famx
  
  
  ## Step X - update baseline cumulative hazard
  updates_impdata <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                              design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx,
                              parms = c(1/41.41327,1,0,0,0, 1)) 
  baseline_gammafr_updates <- as.vector(summary(updates_impdata)$estimates[1:2,1])
  logalpha <- baseline_gammafr_updates[1]
  loglambda <- baseline_gammafr_updates[2]
  
}
end_time <- Sys.time() # End time
computational_time <- end_time - start_time


## Analysis Gamma
gamma_results <- list()
for (i in 1:20) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                           design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx_newx_imp_fam[[i]],
                           parms = c(1/41.41327,1,0,0,0, 1)) 
}





###############################################################################################
########### Multiple Imputation without Considering Family Structure ##########################
###############################################################################################

## Step 1 - Empirical estimates
imp_model <- lm(newx ~ gender + proband + mgene + log(time) + status, data = miss50_famx) 
summary(imp_model)
betas <- coef(imp_model)
sigmahat <- sigma(imp_model)
V <- vcov(imp_model)
SSE <- sum((imp_model$residuals)^2)

fam2_imp <- list()
for (i in 1:20) {
  ## Step 2 - g ~ chi^2 nobs-p
  g <- rchisq(n = 1, df = 1651)
  
  ## Step 3 - sigma star
  sigmastar <- sigmahat/(SSE/g)
  
  ## Step 4 - u1
  u1 <- c()
  u1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  
  ## Step 5
  betastar <- betas + (sigmastar/sigmahat) * u1 %*% chol(V)
  
  ## Step 6
  #num_miss <- sum(is.na(brca1_prs$PRS))
  #u2i <- rnorm(n = 1, mean = 0, sd = 1)
  
  ## Step 7
  fam2_imp[[i]] <- miss20_fam2 |>
    mutate(newx_I = ifelse(is.na(newx), betastar[,1] + betastar[,2]*gender + betastar[,3]*proband + betastar[,4]*mgene + betastar[,5]*log(time) + betastar[,6]*status + rnorm(n = 1, mean = 0, sd = 1)*sigmastar, newx))
}

gamma_results_nofamstruc <- list()
for (i in 1:20) {
  gamma_results_nofamstruc[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                                 design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = fam2_imp[[i]],
                                 parms = c(1/41.41327,1,0,0,0, 1)) 
}

results2 <- colMeans(do.call(rbind, gamma_results_nofamstruc))
nomiss_gamma$par
miss20_gamma_CCA$par
results1

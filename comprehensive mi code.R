####################################################################################
############# Multiple Imputations WITHOUT kinship using log(T) - draft ############
####################################################################################
################################# Imputation Step ##################################
####################################################################################

## Step 1 - Empirical estimates
imp_model <- lm(newx ~  ageonset + log(time) + status  , data = miss50_famx) 
X <- model.matrix( ~ ageonset + log(time) + status  , data = miss50_famx)
summary(imp_model)
betas <- coef(imp_model)
sigmahat <- sigma(imp_model)
V <- vcov(imp_model)
SSE <- sum((imp_model$residuals)^2)
Y_obs <- miss50_famx$newx

miss50_famx_imp <- list()
for (i in 1:5) {
  ## Step 2 - g ~ chi^2 nobs-p
  g <- rchisq(n = 1, df = summary(imp_model)$df[2])
  
  ## Step 3 - sigma star
  sigmastar <- sigmahat/(SSE/g)
  
  ## Step 4 - u1
  #u1 <- c()
  #u1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  mean_vector <- rep(0, summary(imp_model)$df[1])
  Iden_mat <- diag(summary(imp_model)$df[1])
  
  ## Step 5
  betastar <- betas + (sigmastar/sigmahat) * MASS::mvrnorm(n = 1, mu = mean_vector, Sigma = Iden_mat) %*% chol(V)
  betastar <- as.vector(betastar)
  
  ## Step 6
  #num_miss <- sum(is.na(brca1_prs$PRS))%
  #u2i <- rnorm(n = 1, mean = 0, sd = 1)
  
  newx_Imp <- X %*% betastar 
  newx_Imp <- as.vector(newx_Imp)
  #newx_Imp <- sapply(newx_Imp, find_closest, Y_obs)  # Comment out for regular MI
  
  ## Step 7
  miss50_famx_imp[[i]] <- miss50_famx |>
    mutate(newx_Imp_temp = newx_Imp) |>
    mutate(newx_I = ifelse(is.na(newx), newx_Imp_temp, newx))
  
  #miss50_famx_imp[[i]] <- miss50_famx |>
  #  mutate(newx_Imp = newx_Imp) |>
  #  mutate(newx_I = ifelse(is.na(newx), newx_Imp, newx))
}

####################################################################################
################################# Analysis Step ####################################
####################################################################################
## Analysis Gamma
library(foreach)
library(doParallel)

for (i in 1:5) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                                 gvar = "mgene", design = "pop", base.dist = "Weibull", 
                                 frailty.dist = "gamma", agemin = min(miss50_famx_imp[[i]]$currentage[miss50_famx_imp[[i]]$status == 1]), 
                                 data = miss50_famx_imp[[i]], parms = c(0.016,3,0, 0, 0,2)) 
}

MI_noK_logT <- gamma_results 
est_list <- lapply(MI_noK_logT, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))

####################################################################################
################################# Pooling Step #####################################
####################################################################################
final_est <- Pooling(model = gamma_results, imputed_data = miss50_famx_imp)


####################################################################################
############# Multiple Imputations WITHOUT kinship using H0(T) - draft #############
####################################################################################
################################# Imputation Step ##################################
####################################################################################
## Step X - Initial H0(T)
miss50_gamma_cca <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx,
                             parms = c(1/41.41327,1,0,0,0, 1))
summary(miss50_gamma_cca)
baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
miss50_famx_H0 <- miss50_famx |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0

Y_obs <- miss50_famx_H0$newx
miss50_famx_imp <- list()
for (i in 1:5) {
  ## Step 1 - Empirical estimates
  imp_model <- lm(newx ~ gender + currentage + log(H0) + status + mgene, data = miss50_famx_H0) 
  X <- model.matrix( ~ gender + currentage + log(H0) + status + mgene, data = miss50_famx_H0)
  
  betas <- coef(imp_model)
  sigmahat <- sigma(imp_model)
  V <- vcov(imp_model)
  SSE <- sum((imp_model$residuals)^2)
  
  
  ## Step 2 - g ~ chi^2 nobs-p
  g <- rchisq(n = 1, df = summary(imp_model)$df[2])
  
  ## Step 3 - sigma star
  sigmastar <- sigmahat/(SSE/g)
  
  ## Step 4 - u1
  #u1 <- c()
  #u1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  mean_vector <- rep(0, summary(imp_model)$df[1])
  Iden_mat <- diag(summary(imp_model)$df[1])
  
  ## Step 5
  betastar <- betas + (sigmastar/sigmahat) * MASS::mvrnorm(n = 1, mu = mean_vector, Sigma = Iden_mat) %*% chol(V)
  betastar <- as.vector(betastar)
  
  ## Step 6
  #num_miss <- sum(is.na(brca1_prs$PRS))%
  #u2i <- rnorm(n = 1, mean = 0, sd = 1)
  newx_Imp <- X %*% betastar 
  newx_Imp <- as.vector(newx_Imp)
  #newx_Imp <- sapply(newx_Imp, find_closest, Y_obs) # Comment out for regular MI
  
  ## Step 7
  miss50_famx_H0 <- miss50_famx_H0 |>
    mutate(newx_I = newx_Imp) |>
    mutate(newx_I = ifelse(is.na(newx), newx_I, newx)) |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)))
  miss50_famx_imp[[i]] <- miss50_famx_H0
  ## Step X - update baseline cumulative hazard
  if (i == 5) break
  updates_impdata <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                              design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                              agemin = min(miss50_famx_H0$currentage[miss50_famx_H0$status == 1]), 
                              data = miss50_famx_H0,
                              parms = c(0.016,3,0, 0, 0,2)) 
  baseline_gammafr_updates <- as.vector(summary(updates_impdata)$estimates[1:2,1])
  logalpha <- baseline_gammafr_updates[1]
  loglambda <- baseline_gammafr_updates[2]
  
}

####################################################################################
################################# Analysis Step ####################################
####################################################################################
## Analysis Gamma
library(foreach)
library(doParallel)

for (i in 1:5) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                                 gvar = "mgene", design = "pop", base.dist = "Weibull", 
                                 frailty.dist = "gamma", agemin = min(miss50_famx_imp[[i]]$currentage[miss50_famx_imp[[i]]$status == 1]), 
                                 data = miss50_famx_imp[[i]], parms = c(0.016,3,0, 0, 0,2)) 
}
est_list <- lapply(gamma_results, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))
stopCluster(cl)

####################################################################################
############### Multiple Imputation with kinship using log(T) - draft ##############
####################################################################################
################################# Imputation Step ##################################
####################################################################################
## Step 1 - Kinship matrix
library(kinship2)
library(lme4qtl)
library(lme4)
kinship_mat <- with(miss50_famx, kinship(id = indID, dadid = fatherID, momid = motherID,
                                         sex = gender))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(miss50_famx)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step 2 - empirical estimates
model_test <- relmatLmer(newx ~ gender + ageonset + currentage + log(time) + status + mgene + (1|indID), data = miss50_famx, relmat = list(indID = kinship_mat))
summary(model_test)
X <- model.matrix(~ gender + ageonset + currentage + log(time) + status + mgene, data = miss50_famx) # imputation model design matrix

#betas <- coef(model_test)
betas <- as.vector(summary(model_test)$coefficients[,1]) # beta coefficients
sigma_g_2 <- attr(lme4::VarCorr(model_test)$indID, "stddev")^2 # genetic variance
sigma_e_2 <- attr(lme4::VarCorr(model_test), "sc")^2 # residual variance

Sigma <- sigma_g_2*kinship_mat_sparse + sigma_e_2*Iden_mat_sparse # Sparse matrix for Sigma
Sigma_mat <- as.matrix(Sigma) # Non-sparse

V <- vcov(model_test)

## Step 3 - conditional variance
#num_cores <- detectCores() - 2 # 6 cores
#cl <- makeCluster(num_cores)
#clusterExport(cl, varlist = c("Sigma")) 
#clusterEvalQ(cl, library(Matrix))

#cond_var <- function(i) {
#  conditional_var <- Sigma[i,i] - Sigma[i,-i] %*% solve(Sigma[-i,-i]) %*% Sigma[-i,i]
#  return(conditional_var)
#}
#conditional_variances <- parSapply(cl, 1:nrow(miss50_famx), cond_var)
#conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))

#stopCluster(cl)

## Imputation steps
miss50_famx <- miss50_famx |>
  mutate(cond_var = conditional_variances_temp)

newx <- miss50_famx$newx
Y_obs <- miss50_famx$newx
famx_newx_imp_fam <- list()

for (m in 1:5) {
  
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
    
    conditional_Expect <- mu_star[i] + Sigma[i,-i][non_NA] %*% solve(Sigma[-i,-i][non_NA, non_NA]) %*% (y_minus_i - mu_star_minus_i)
    return(conditional_Expect)
  }
  conditional_expectations <- parSapply(cl, 1:nrow(miss50_famx), E_cond)
  conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
  
  #conditional_vars <- miss50_famx$cond_var
  
  ## Step 8
  w2i <- rnorm(n = nrow(miss50_famx), mean = 0, sd = 1)
  
  ## Step 9
  newx_star <- conditional_expectations 
  #newx_star <- sapply(newx_star, find_closest, Y_obs) # Comment out for regular MI
  #newx_star <- conditional_expectations + conditional_vars
  famx_newx_imp_fam[[m]] <- miss50_famx |>
    mutate(newx_I = newx_star) |>
    mutate(newx_I = ifelse(!is.na(newx), newx, newx_I))
  
  stopCluster(cl)
}

####################################################################################
################################# Analysis Step ####################################
####################################################################################
gamma_results <- list()
for (i in 1:5) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                                 gvar = "mgene", design = "pop", base.dist = "Weibull", 
                                 frailty.dist = "gamma", agemin = min(miss50_famx_imp[[i]]$currentage[miss50_famx_imp[[i]]$status == 1]), 
                                 data = famx_newx_imp_fam[[i]], parms = c(0.016,3,0, 0, 0,2)) 
}

cl <- 6
registerDoParallel(cl)
gamma_results <- foreach(i = 1:5) %dopar% {
  penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
           gvar = "mgene", design = "pop", base.dist = "Weibull", 
           frailty.dist = "gamma", agemin = 20, 
           data = famx_newx_imp_fam[[i]], parms = c(0.016,3,0, 0, 0,2))
}

mi_K_logT <- gamma_results
est_list <- lapply(mi_K_logT, function(x) summary(x)$estimates[,1]) 
SE_list <- lapply(mi_K_logT, function(x) x$se)
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))

famx_newx_imp_fam[[1]] |>
  ggplot(aes(x = newx_I)) + geom_density()

####################################################################################
##### Multiple Imputations with kinship using H0(T) Baseline cumulative hazard #####
####################################################################################
################################# Imputation Step ##################################
####################################################################################
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
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                             agemin = min(miss50_famx$currentage[miss50_famx$status == 1]), data = miss50_famx,
                             parms = c(1/41.41327,1,0,0,0, 1))
baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
miss50_famx <- miss50_famx |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
Y_obs <- miss50_famx$newx

start_time <- Sys.time() # Starting time 
for (m in 1:5) {
  ## Step 2 - empirical estimates
  model_test <- relmatLmer(newx ~ gender + ageonset + currentage + log(H0) + status + mgene + (1|indID), data = miss50_famx, relmat = list(indID = kinship_mat))
  #summary(model_test)
  X <- model.matrix(~  gender + ageonset + currentage + log(H0) + status + mgene, data = miss50_famx) # imputation model design matrix
  V <- vcov(model_test)
  
  #betas <- coef(model_test)
  betas <- as.vector(summary(model_test)$coefficients[,1]) # beta coefficients
  sigma_g_2 <- attr(VarCorr(model_test)$indID, "stddev")^2 # genetic variance
  sigma_e_2 <- attr(VarCorr(model_test), "sc")^2 # residual variance
  
  Sigma <- sigma_g_2*kinship_mat_sparse + sigma_e_2*Iden_mat_sparse # Sparse matrix for Sigma
  Sigma_mat <- as.matrix(Sigma) # Non-sparse
  
  ## Step 3 - conditional variance
  #num_cores <- detectCores() - 2 # 6 cores
  #cl <- makeCluster(num_cores)
  #clusterExport(cl, varlist = c("Sigma")) 
  #clusterEvalQ(cl, library(Matrix))
  
  #cond_var <- function(i) {
  #  conditional_var <- Sigma[i,i] - Sigma[i,-i] %*% solve(Sigma[-i,-i]) %*% Sigma[-i,i]
  #  return(conditional_var)
  #}
  #conditional_variances <- parSapply(cl, 1:nrow(miss50_famx), cond_var)
  #conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))
  
  #stopCluster(cl)
  
  ## Imputation steps
  #miss50_famx <- miss50_famx |>
  #  mutate(cond_var = conditional_variances_temp)
  
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
  
  #conditional_vars <- miss50_famx$cond_var
  stopCluster(cl)
  
  ## Step 8
  #w2i <- rnorm(n = nrow(miss50_famx), mean = 0, sd = 1)
  
  ## Step 9
  newx_star <- conditional_expectations 
  #newx_star <- sapply(newx_star, find_closest, Y_obs) # Comment out for regular MI
  
  miss50_famx <- miss50_famx |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)) ) |>
    mutate(newx_I = newx_star) |>
    mutate(newx_I = ifelse(!is.na(newx), newx, newx_I)) |>
    group_by(famID) |>
    mutate(newx_I = ifelse(!is.na(newx), find_closest(newx_I, newx), newx_I)) |>
    ungroup() 
    
  famx_newx_imp_fam[[m]] <- miss50_famx
  
  if (m >= 5) break
  ## Step X - update baseline cumulative hazard
  updates_impdata <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                              gvar = "mgene", 
                              design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                              agemin = min(miss50_famx$currentage[miss50_famx$status == 1]), 
                              data = miss50_famx,
                              parms = c(1/41.41327,1,0,0,0, 1)) 
  baseline_gammafr_updates <- as.vector(summary(updates_impdata)$estimates[1:2,1])
  logalpha <- baseline_gammafr_updates[1]
  loglambda <- baseline_gammafr_updates[2]
  
}
end_time <- Sys.time() # End time
computational_time <- end_time - start_time

####################################################################################
################################# Analysis Step ####################################
####################################################################################
## Analysis Gamma
gamma_results <- list()
for (i in 1:5) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                                 design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                                 agemin = min(miss50_famx_imp[[i]]$currentage[miss50_famx_imp[[i]]$status == 1]), data = famx_newx_imp_fam[[i]],
                                 parms = c(1/41.41327,1,0,0,0, 1)) 
}
MI_K_H0_gamma_results <- gamma_results
est_list <- lapply(gamma_results, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))

famx_newx_imp_fam[[5]] |>
  ggplot(aes(x = newx_I)) + geom_density()
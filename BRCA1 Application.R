###############################################################
###################### BRCA1 Data Application #################
####################### Complete Case Analysis ################
###############################################################
CCA_model_gamma <- penmodel(Surv(timeBC, BC) ~ mgeneI + PRS, cluster = "famID", 
                      gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
                      frailty.dist = "gamma", agemin = min(brca1_prs$timeBC[status == 1]), 
                      data = brca1_prs, parms = c(0.016,3,0, 0,2))
summary(CCA_model_gamma)
CCA_model_lognormal <- penmodel(Surv(timeBC, BC) ~ mgeneI + PRS, cluster = "famID", 
                                gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
                                frailty.dist = "lognormal", agemin = min(brca1_prs$timeBC[status == 1]), 
                                data = brca1_prs, parms = c(0.016,3,0, 0,2))
summary(CCA_model_lognormal)

coxph_1 <- coxph(Surv(timeBC, BC) ~ mgeneI + PRS, data = brca1_prs[brca1_prs$proband == 0,])
summary(coxph_1)





###############################################################
###################### BRCA1 Data Application #################
####################### MI without K log(T) ###################
###############################################################
## Step 1 - Empirical estimates
imp_model <- lm(PRS ~ log(timeBC) + BC + mgeneI + proband , data = brca1_prs) 
X <- model.matrix( ~ log(timeBC) + BC + mgeneI + proband , data = brca1_prs)
summary(imp_model)
betas <- coef(imp_model)
sigmahat <- sigma(imp_model)
V <- vcov(imp_model)
SSE <- sum((imp_model$residuals)^2)
Y_obs <- brca1_prs$PRS

brca1_prs_imp <- list()
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
  
  PRS_Imp <- X %*% betastar 
  PRS_Imp <- as.vector(PRS_Imp)
  #PRS_Imp <- sapply(PRS_Imp, find_closest, Y_obs)
  
  ## Step 7
  brca1_prs_imp[[i]] <- brca1_prs |>
    mutate(PRS_Imp_temp = PRS_Imp) |>
    mutate(PRS_I = ifelse(is.na(PRS), PRS_Imp_temp, PRS))
  
  #brca1_prs_imp[[i]] <- brca1_prs |>
  #  mutate(PRS_Imp = PRS_Imp) |>
  #  mutate(PRS_I = ifelse(is.na(PRS), PRS_Imp, PRS))
}

####################################################################################
################################# Analysis Step ####################################
####################################################################################
## Analysis Gamma
library(foreach)
library(doParallel)

num_cores <- detectCores() - 2 # 6 cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

gamma_results <- foreach(i = 1:5, .packages = c("survival")) %dopar% {
  penmodel(Surv(timeBC, BC) ~ mgeneI + PRS_I, cluster = "famID", 
           gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
           frailty.dist = "gamma", agemin = 20, 
           data = brca1_prs_imp[[i]], parms = c(0.016,3,0, 0,2)) 
}
MI_noK_logT <- gamma_results 
est_list <- lapply(MI_noK_logT, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))
stopCluster(cl)

lognormal_results <- foreach(i = 1:5, .packages = c("survival")) %dopar% {
  penmodel(Surv(timeBC, BC) ~ mgeneI + PRS_I, cluster = "famID", 
           gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
           frailty.dist = "lognormal", agemin = 20, 
           data = brca1_prs_imp[[i]], parms = c(0.016,3,0, 0,2)) 
}
MI_noK_logT <- lognormal_results 
est_list <- lapply(MI_noK_logT, function(x) summary(x)$estimates[,1]) 
est_lognormal_aftermi <- colMeans(do.call(rbind, est_list))
stopCluster(cl)


#################### Check Zone ####################
check_model <- penmodel(Surv(time, status) ~ gender + mgene + PRS_I, cluster = "famID", 
                        gvar = "mgene", design = "pop", base.dist = "Weibull", 
                        frailty.dist = "gamma", agemin = 20, 
                        data = brca1_prs_imp[[3]], parms = c(0.45,1,1,1,1, 1)) 
summary(check_model)

brca1_prs_imp[[5]] |> ggplot(aes(x = PRS_I)) + geom_density()

###############################################################
#################### BRCA1 Data Application ###################
#################### MI without K H0(T) #######################
###############################################################
## Step X - Initial H0(T)
brca1_gamma_cca <- penmodel(Surv(time, status) ~ mgeneI + PRS, cluster = "famID", gvar = "mgeneI", 
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                             agemin = min(brca1_prs$timeBC[status == 1]), data = brca1_prs,
                             parms = c(0.016,3,0, 0,2))
summary(brca1_gamma_cca)
baseline_gammafr <- as.vector(summary(brca1_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
brca1_prs_H0 <- brca1_prs |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0

Y_obs <- brca1_prs_H0$PRS
brca1_imp <- list()
for (i in 1:5) {
  ## Step 1 - Empirical estimates
  imp_model <- lm(PRS ~ log(H0) + mgeneI + status + proband, data = brca1_prs_H0) 
  X <- model.matrix( ~ log(H0) + mgeneI + status + proband, data = brca1_prs_H0)
  
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
  PRS_Imp <- X %*% betastar 
  PRS_Imp <- as.vector(PRS_Imp)
  #PRS_Imp <- sapply(PRS_Imp, find_closest, Y_obs)
  
  ## Step 7
  brca1_prs_H0 <- brca1_prs_H0 |>
    mutate(PRS_I = PRS_Imp) |>
    mutate(PRS_I = ifelse(is.na(PRS), PRS_I, PRS)) |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)))
  brca1_imp[[i]] <- brca1_prs_H0
  ## Step X - update baseline cumulative hazard
  if (i == 5) break
  updates_impdata <- penmodel(Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", gvar = "mgeneI", 
                              design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                              agemin = min(brca1_prs$timeBC[status == 1]), data = brca1_imp[[i]],
                              parms = c(0.016,3,0, 0,2)) 
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
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", 
                                 gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
                                 frailty.dist = "gamma", agemin = min(brca1_prs$timeBC[status == 1]), 
                                 data = brca1_imp[[i]], parms = c(0.016,1,0, 0, 2)) 
}
est_list <- lapply(gamma_results, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))
stopCluster(cl)

###############################################################
###################### BRCA1 Data Application #################
####################### MI with K log(T) ######################
###############################################################
## Step 1 - Kinship matrix
library(kinship2)
library(lme4qtl)
library(lme4)
kinship_mat <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID,
                                         sex = 0))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(brca1_prs)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step 2 - empirical estimates
model_test <- relmatLmer(PRS ~ log(time) + status + proband + mgeneI + (1|indID), data = brca1_prs, relmat = list(indID = kinship_mat))
summary(model_test)
X <- model.matrix(~ log(time) + status + proband + mgeneI, data = brca1_prs) # imputation model design matrix

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
#conditional_variances <- parSapply(cl, 1:nrow(brca1_prs), cond_var)
#conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))

#stopCluster(cl)

## Imputation steps
#brca1_prs <- brca1_prs |>
#  mutate(cond_var = conditional_variances_temp)

PRS <- brca1_prs$PRS
Y_obs <- brca1_prs$PRS
famx_PRS_imp_fam <- list()

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
  clusterExport(cl, varlist = c("Sigma", "PRS", "mu_star")) 
  clusterEvalQ(cl, library(Matrix))
  
  E_cond <- function(i) {
    y_minus_i <- PRS[-i]
    mu_star_minus_i <- mu_star[-i]
    non_NA <- which(!is.na(y_minus_i))
    y_minus_i <- y_minus_i[non_NA]
    mu_star_minus_i <- mu_star_minus_i[non_NA]
    
    conditional_Expect <- mu_star[i] + Sigma[i,-i][non_NA] %*% solve(Sigma[-i,-i][non_NA, non_NA]) %*% (y_minus_i - mu_star_minus_i)
    return(conditional_Expect)
  }
  conditional_expectations <- parSapply(cl, 1:nrow(brca1_prs), E_cond)
  conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
  
  #conditional_vars <- brca1_prs$cond_var
  
  ## Step 8
  w2i <- rnorm(n = nrow(brca1_prs), mean = 0, sd = 1)
  
  ## Step 9
  PRS_star <- conditional_expectations 
  #PRS_star <- sapply(PRS_star, find_closest, Y_obs)
  #PRS_star <- conditional_expectations + conditional_vars
  famx_PRS_imp_fam[[m]] <- brca1_prs |>
    mutate(PRS_I = PRS_star) |>
    mutate(PRS_I = ifelse(!is.na(PRS), PRS, PRS_I))
  
  stopCluster(cl)
}

####################################################################################
################################# Analysis Step ####################################
####################################################################################
gamma_results <- list()
for (i in 1:5) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", 
                                 gvar = "mgeneI", design = "pop", base.dist = "Weibull", 
                                 frailty.dist = "gamma", agemin = 20, 
                                 data = famx_PRS_imp_fam[[i]], parms = c(0.016,3,0, 0,2)) 
}
mi_K_logT <- gamma_results
est_list <- lapply(mi_K_logT, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))

###############################################################
###################### BRCA1 Data Application #################
####################### MI with K H0(T) #######################
###############################################################
famx_PRS_imp_fam <- list() 
## Step 1 - Kinship matrix
kinship_mat <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID,
                                         sex = gender))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(brca1_prs)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step X - Initial baseline cumulative hazard
brca1_gamma_cca <- penmodel(Surv(time, status) ~ mgeneI + PRS, cluster = "famID", gvar = "mgeneI", 
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                             agemin = min(brca1_prs$timeBC[status == 1]), data = brca1_prs,
                             parms = c(1/41.41327,1,0,0, 1))
baseline_gammafr <- as.vector(summary(brca1_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
brca1_prs <- brca1_prs |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
Y_obs <- brca1_prs$PRS

start_time <- Sys.time() # Starting time 
for (m in 1:5) {
  ## Step 2 - empirical estimates
  model_test <- relmatLmer(PRS ~ log(H0) + status + mgeneI + proband + (1|indID), data = brca1_prs, relmat = list(indID = kinship_mat))
  #summary(model_test)
  X <- model.matrix(~  log(H0) + status + mgeneI + proband, data = brca1_prs) # imputation model design matrix
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
  #conditional_variances <- parSapply(cl, 1:nrow(brca1_prs), cond_var)
  #conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))
  
  #stopCluster(cl)
  
  ## Imputation steps
  #brca1_prs <- brca1_prs |>
  #  mutate(cond_var = conditional_variances_temp)
  
  PRS <- brca1_prs$PRS
  
  
  ## Step 4
  w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  
  ## Step 5
  betastar <- betas + w_1 %*% chol(V)
  
  ## Step 6
  mu_star <- X %*% as.vector(betastar)
  
  ## Step 7
  num_cores <- detectCores() - 2 # 6 cores
  cl <- makeCluster(num_cores)
  clusterExport(cl, varlist = c("Sigma", "PRS", "mu_star")) 
  clusterEvalQ(cl, library(Matrix))
  
  E_cond <- function(i) {
    y_minus_i <- PRS[-i]
    mu_star_minus_i <- mu_star[-i]
    non_NA <- which(!is.na(y_minus_i))
    y_minus_i <- y_minus_i[non_NA]
    mu_star_minus_i <- mu_star_minus_i[non_NA]
    Sigma_nonNA <- Sigma[-i,-i][non_NA, non_NA]
    Sigmahat <- Sigma[i,-i][non_NA]
    
    conditional_Expect <- mu_star[i] + Sigmahat %*% solve(Sigma_nonNA) %*% (y_minus_i - mu_star_minus_i)
    return(conditional_Expect)
  }
  conditional_expectations <- parSapply(cl, 1:nrow(brca1_prs), E_cond)
  conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
  
  #conditional_vars <- brca1_prs$cond_var
  stopCluster(cl)
  
  ## Step 8
  #w2i <- rnorm(n = nrow(brca1_prs), mean = 0, sd = 1)
  
  ## Step 9
  PRS_star <- conditional_expectations 
  #PRS_star <- sapply(PRS_star, find_closest, Y_obs)
  
  brca1_prs <- brca1_prs |>
    mutate(PRS_I = PRS_star) |>
    mutate(PRS_I = ifelse(!is.na(PRS), PRS, PRS_I)) |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)) )
  famx_PRS_imp_fam[[m]] <- brca1_prs
  
  if (m >= 5) break
  ## Step X - update baseline cumulative hazard
  updates_impdata <- penmodel(Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", gvar = "mgeneI", 
                              design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                              agemin = min(brca1_prs$timeBC[status == 1]), data = brca1_prs,
                              parms = c(1/41.41327,1,0,0, 1)) 
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
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ mgeneI + PRS_I, cluster = "famID", gvar = "mgeneI", 
                                 design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                                 agemin = min(brca1_prs$timeBC[status == 1]), data = famx_PRS_imp_fam[[i]],
                                 parms = c(1/41.41327,1,0,0, 1)) 
}
MI_K_H0_gamma_results <- gamma_results
est_list <- lapply(gamma_results, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))

famx_PRS_imp_fam[[5]] |>
  ggplot(aes(x = PRS_I)) + geom_density()


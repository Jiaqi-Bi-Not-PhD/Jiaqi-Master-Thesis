##########################################################
########### Multiple Imputation ##########################
##########################################################
## Kinship matrix
devtools::install_github("variani/lme4qtl")
library(kinship2)
library(Matrix)
library("lme4qtl")
library(lmerTest)
library(parallel)
library(MASS)

kinship_mat <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID,
                                       sex = rep(2,nrow(brca1_prs))))
kinship_mat_sparse <- as(kinship_mat, "sparseMatrix")
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

## Linear mixed effect model using kinship matrix
model_test <- relmatLmer(PRS ~ proband + mgeneI + currentage + timeBC * BC + (1|indID), 
                         data = brca1_prs, relmat = list(indID = kinship_mat))
summary(model_test)
brca1_prs_MI <- brca1_prs |>
  mutate(mu_PRS = predict(model_test, newdata = brca1_prs, re.form = NA))

## Constructing kinship and covariance
sigma_g <- attr(VarCorr(model_test)$indID, "stddev") # genetic variance
sigma_e <- attr(VarCorr(model_test), "sc") # residual
Sigma <- sigma_g^2 * as(kinship_mat, "sparseMatrix") + sigma_e^2 # Is this formula correct?
Sigma <- as.matrix(Sigma)
kinship_ids <- rownames(Sigma)

## Imputation step
#numCores <- detectCores() - 2
imputed_data_list <- list()

missing_indices <- brca1_prs[which(brca1_prs$miss_index == 1), "indID"]
missing_indices <- as.vector(missing_indices)$indID

for (m in 1:20) {
  brca1_prs_MI_completed <- brca1_prs_MI
  for (i in missing_indices) {
    indID <- brca1_prs_MI_completed$indID[brca1_prs_MI_completed$indID == i]
    mu_PRS <- brca1_prs_MI_completed$mu_PRS[brca1_prs_MI_completed$indID == i]
    Sigma_index <- match(indID, kinship_ids)
    Sigma_i <- Sigma[Sigma_index, Sigma_index]
  
    imputed_value <- mvrnorm(n = 1, mu = mu_PRS, Sigma = matrix(Sigma_i, nrow = 1)) # change to family level
    brca1_prs_MI_completed$PRS_I[brca1_prs_MI_completed$indID == i] <- imputed_value # Gibb's sampling
  }
  brca1_prs_MI_completed <- brca1_prs_MI_completed |>
    mutate(PRS_I = ifelse(is.na(PRS), PRS_I, PRS))
  imputed_data_list[[m]] <- brca1_prs_MI_completed
}

## Analysis step (log-normal)
initial_params <- c(1/41.41327,1, 0, 0, 1)
results_list <- list()
for (m in 1:20) {
X <- as.matrix(data.frame(imputed_data_list[[m]]$mgeneI, imputed_data_list[[m]]$PRS_I), 
               nrow=nrow(imputed_data_list[[m]]), 
               ncol = 2)
Y <- as.matrix(data.frame(imputed_data_list[[m]]$timeBC, imputed_data_list[[m]]$BC), 
               nrow = nrow(imputed_data_list[[m]]), 
               ncol = 2)

results_list[[m]] <- optim(par = initial_params, fn = lognormal_single,
      data = imputed_data_list[[m]], X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
      agemin = 18, control = list(maxit = 2000))
}

results_list[[5]]

###############################################################################################
########### Multiple Imputation without Considering Family Structure ##########################
###############################################################################################

## Step 1 - Empirical estimates
imp_model <- lm(PRS ~ proband + proband:currentage + mgeneI + log(timeBC)*BC, data = brca1_prs) 
summary(imp_model)
betas <- coef(imp_model)
sigmahat <- sigma(imp_model)
V <- vcov(imp_model)
SSE <- sum((imp_model$residuals)^2)

brca1_prs_imp <- list()
for (i in 1:20) {
  ## Step 2 - g ~ chi^2 nobs-p
  g <- rchisq(n = 1, df = 497)
  
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
  brca1_prs_imp[[i]] <- brca1_prs |>
    mutate(PRS_I = ifelse(is.na(PRS), betastar[,1] + betastar[,2]*proband + betastar[,3]*mgeneI + betastar[,4]*log(timeBC) + betastar[,5]*BC + betastar[,6]*proband*currentage + betastar[,7]*log(timeBC)*BC + rnorm(n = 1, mean = 0, sd = 1)*sigmastar, PRS))
}

## Analysis log-normal
log_norm_results <- list()
for (i in 1:20) {
  X <- as.matrix(data.frame(brca1_prs_imp[[i]]$mgeneI, brca1_prs_imp[[i]]$PRS_I), 
                 nrow=nrow(brca1_prs_imp[[i]]), 
                 ncol = 2)
  Y <- as.matrix(data.frame(brca1_prs_imp[[i]]$timeBC, brca1_prs_imp[[i]]$BC), 
                 nrow = nrow(brca1_prs_imp[[i]]), 
                 ncol = 2)
  
  initial_params <- c(1/41.41327,1, 0, 0, 1)
  log_norm_forgraph <- optim(par = initial_params, fn = lognormal_single,
                             data = brca1_prs_imp[[i]], X = X, Y = Y, nbase = 2,
                             design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
                             agemin = 18, control = list(maxit = 2000))
  log_norm_results[[i]] <- log_norm_forgraph$par
}

mean(sapply(log_norm_results, function(x) x[1]))
mean(sapply(log_norm_results, function(x) x[2]))
mean(sapply(log_norm_results, function(x) x[3]))
mean(sapply(log_norm_results, function(x) x[4]))
mean(sapply(log_norm_results, function(x) x[5]))


## Analysis Gamma
gamma_results <- list()
for (i in 1:20) {
  initial_params <- c(1/41.41327,1,0,0, 1)
  X <- as.matrix(brca1_prs_imp[[i]][,c("mgeneI", "PRS_I")], ncol = 2)
  Y <- as.matrix(brca1_prs_imp[[i]][,c("timeBC", "BC")], ncol = 2)
  gamma_forgraph <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                          data = brca1_prs_imp[[i]], X = X, Y = Y, nbase = 2,
                          design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                          agemin = 18, 
                          control = list(maxit = 2000))
  gamma_results[[i]] <- gamma_forgraph$par
}

mean(sapply(gamma_results, function(x) x[1]))
mean(sapply(gamma_results, function(x) x[2]))
mean(sapply(gamma_results, function(x) x[3]))
mean(sapply(gamma_results, function(x) x[4]))
mean(sapply(gamma_results, function(x) x[5]))

colMeans(do.call(rbind, gamma_results))

## Analysis coxph
coxph_results <- list()
for (i in 1:20) {
  m <- coxph(Surv(timeBC, BC) ~ mgeneI + PRS_I + frailty(famID, distribution = "gamma"), data = brca1_prs_imp[[i]])
  coxph_results[[i]] <- m$coefficients
}

mean(sapply(coxph_results, function(x) x[1]))
mean(sapply(coxph_results, function(x) x[2]))

## Analysis coxme
coxme_results <- list()
for (i in 1:20) {
  m <- coxme(Surv(timeBC, BC) ~ mgeneI + PRS_I + (1|famID), data = brca1_prs_imp[[i]])
  coxme_results[[i]] <- m$coefficients
}

mean(sapply(coxme_results, function(x) x[1]))
mean(sapply(coxme_results, function(x) x[2]))

###############################################################################################
########### Multiple Imputation Considering Family Structure ##################################
###############################################################################################
library(lme4qtl)
library(kinship2)
library(Matrix)
library(MASS)
library(lmerTest)
library(parallel)
library(tidyverse)
library(survival)

brca1_prs <- brca1_prs |>
  mutate(index = seq(1, 2650, 1))

## Step 1 - Kinship matrix
kinship_mat <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID,
                                       sex = rep(2,nrow(brca1_prs))))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(brca1_prs)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step 2 - empirical estimates
model_test <- relmatLmer(PRS ~ proband + proband:currentage + mgeneI + log(timeBC)*BC + (1|indID), data = brca1_prs, relmat = list(indID = kinship_mat))
summary(model_test)
X <- model.matrix( ~ proband + proband:currentage + mgeneI + log(timeBC)*BC, data = brca1_prs) 

#betas <- coef(model_test)
betas <- c(0.309458, 0.190865, -0.238208, -0.115634, -0.835923, -0.002240, 0.233941)
sigma_g_2 <- attr(VarCorr(model_test)$indID, "stddev")^2 # genetic variance
sigma_e_2 <- attr(VarCorr(model_test), "sc")^2 # residual variance

Sigma <- sigma_g_2*kinship_mat_sparse + sigma_e_2*Iden_mat_sparse # Sparse matrix for Sigma
Sigma_mat <- as.matrix(Sigma) # Non-sparse

V <- vcov(model_test)

## Step 3 - conditional variance
num_cores <- detectCores() - 2 # 6 cores
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("Sigma")) 
clusterEvalQ(cl, library(Matrix))

cond_var <- function(i) {
  conditional_var <- Sigma[i,i] - Sigma[i,-i] %*% solve(Sigma[-i,-i]) %*% Sigma[-i,i]
  return(conditional_var)
}
conditional_variances <- parSapply(cl, 1:nrow(brca1_prs), cond_var)
conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))

stopCluster(cl) # Free up memory after everything

## Imputation steps
brca1_prs <- brca1_prs |>
  mutate(cond_var = conditional_variances_temp)
  
brca1_prs_imp_fam <- list() 
for (m in 1:20) {
  ## Step 4
  w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  
  ## Step 5
  betastar <- betas + w_1 %*% chol(V) # w_1 ~ N(0,1) always?
  
  ## Step 6
  mu_star <- X %*% as.vector(betastar)
  
  ## Step 7
  PRS <- brca1_prs$PRS
  
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
  
  conditional_vars <- brca1_prs$cond_var
  
  ## Step 8
  w2i <- rnorm(n = nrow(brca1_prs), mean = 0, sd = 1)
  
  ## Step 9
  PRS_star <- conditional_expectations + w2i * conditional_vars
  brca1_prs_imp_fam[[m]] <- brca1_prs |>
    mutate(PRS_I = PRS_star) |>
    mutate(PRS_I = ifelse(!is.na(PRS), PRS, PRS_I))
  
  stopCluster(cl)
}

## Analysis log-normal
log_norm_results <- list()
for (i in 1:20) {
  X <- as.matrix(data.frame(brca1_prs_imp_fam[[i]]$mgeneI, brca1_prs_imp_fam[[i]]$PRS_I), 
                 nrow=nrow(brca1_prs_imp_fam[[i]]), 
                 ncol = 2)
  Y <- as.matrix(data.frame(brca1_prs_imp_fam[[i]]$timeBC, brca1_prs_imp_fam[[i]]$BC), 
                 nrow = nrow(brca1_prs_imp_fam[[i]]), 
                 ncol = 2)
  
  initial_params <- initial_params <- c(1/41.41327,1,0,0, 1)
  log_norm_forgraph <- optim(par = initial_params, fn = lognormal_single,
                             data = brca1_prs_imp_fam[[i]], X = X, Y = Y, nbase = 2,
                             design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
                             agemin = 18, control = list(maxit = 2000))
  log_norm_results[[i]] <- log_norm_forgraph$par
}

mean(sapply(log_norm_results, function(x) x[1]))
mean(sapply(log_norm_results, function(x) x[2]))
mean(sapply(log_norm_results, function(x) x[3]))
mean(sapply(log_norm_results, function(x) x[4]))
mean(sapply(log_norm_results, function(x) x[5]))

colMeans(do.call(rbind, log_norm_results))


## Analysis Gamma
gamma_results <- list()
for (i in 1:20) {
  initial_params <- c(1/41.41327,1,0,0, 1)
  X <- as.matrix(brca1_prs_imp_fam[[i]][,c("mgeneI", "PRS_I")], ncol = 2)
  Y <- as.matrix(brca1_prs_imp_fam[[i]][,c("timeBC", "BC")], ncol = 2)
  gamma_forgraph <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                          data = brca1_prs_imp_fam[[i]], X = X, Y = Y, nbase = 2,
                          design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                          agemin = 18, 
                          control = list(maxit = 2000))
  gamma_results[[i]] <- gamma_forgraph$par
}

mean(sapply(gamma_results, function(x) x[1]))
mean(sapply(gamma_results, function(x) x[2]))
mean(sapply(gamma_results, function(x) x[3]))
mean(sapply(gamma_results, function(x) x[4]))
mean(sapply(gamma_results, function(x) x[5]))

gamma_mi_kinship <- colMeans(do.call(rbind, gamma_results))

## Analysis coxph
coxph_results <- list()
for (i in 1:20) {
  m <- coxph(Surv(timeBC, BC) ~ mgeneI + PRS_I + frailty(famID, distribution = "gamma"), data = brca1_prs_imp_fam[[i]])
  coxph_results[[i]] <- m$coefficients
}

mean(sapply(coxph_results, function(x) x[1]))
mean(sapply(coxph_results, function(x) x[2]))

colMeans(do.call(rbind, coxph_results))

cca <- coxph(Surv(timeBC, BC) ~ mgeneI + PRS + frailty(famID, distribution = "gamma"), data = brca1_prs_cca1)
summary(cca)
brca1_prs_cca1 <- brca1_prs_cca |>
  mutate(proband = 0)

## Analysis coxme
coxme_results <- list()
for (i in 1:20) {
  m <- coxme(Surv(timeBC, BC) ~ mgeneI + PRS_I + (1|famID), data = brca1_prs_imp_fam[[i]])
  coxme_results[[i]] <- m$coefficients
}

mean(sapply(coxme_results, function(x) x[1]))
mean(sapply(coxme_results, function(x) x[2]))

colMeans(do.call(rbind, coxme_results))


#######################
s <- Sigma[3,3] - Sigma[3, -3] %*% solve(Sigma[-3, -3]) %*% Sigma[-3, 3] # Test
Sigma_test <- Sigma[1:6, 1:6]
Sigma_test[3,3] - Sigma_test[3, -3] %*% solve(Sigma_test[-3, -3]) %*% Sigma_test[-3, 3]

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
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx_CCA,
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
  model_test <- relmatLmer(newx ~ gender + ageonset + currentage + H0 + status + mgene + (1|indID), data = miss50_famx, relmat = list(indID = kinship_mat))
  #summary(model_test)
  X <- model.matrix(~  gender + ageonset + currentage + H0 + status + mgene, data = miss50_famx) # imputation model design matrix
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
  newx_star <- sapply(newx_star, find_closest, Y_obs)
  
  miss50_famx <- miss50_famx |>
    mutate(newx_I = newx_star) |>
    mutate(newx_I = ifelse(!is.na(newx), newx, newx_I)) |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)) )
  famx_newx_imp_fam[[m]] <- miss50_famx
  
  if (m >= 5) break
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

####################################################################################
################################# Analysis Step ####################################
####################################################################################
## Analysis Gamma
gamma_results <- list()
for (i in 1:5) {
  gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                                 design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx_newx_imp_fam[[i]],
                                 parms = c(1/41.41327,1,0,0,0, 1)) 
}
MI_K_H0_gamma_results <- gamma_results
est_list <- lapply(gamma_results, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))

famx_newx_imp_fam[[5]] |>
  ggplot(aes(x = newx_I)) + geom_density()




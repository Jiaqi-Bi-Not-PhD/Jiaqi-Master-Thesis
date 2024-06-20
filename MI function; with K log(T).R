####################################################################################
####################################################################################
####################################################################################
################## Multiple Imputations WITH kinship using log(T) ##################
####################################################################################
####################################################################################
####################################################################################

####################################################################################
############### Multiple Imputation with kinship using log(T) - draft ##############
####################################################################################
################################# Imputation Step ##################################
####################################################################################
MI_FamEvent_K_logT <- function(data, M = 5, option = "General", frailty.dist = "gamma") {
  attempts <- 0
  max_attempts <- 3
  success <- FALSE
  while(attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
  ## Step 1 - Kinship matrix
  kinship_mat <- with(data, kinship(id = indID, dadid = fatherID, momid = motherID,
                                           sex = gender))
  kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
  kinship_mat_sparse <- 2 * kinship_mat_sparse
  kinship_mat <- as.matrix(kinship_mat_sparse)
  
  Iden_mat <- diag(nrow = nrow(data)) # Identity matrix n*n
  Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)
  
  ## Step 2 - empirical estimates
  model_test <- relmatLmer(newx ~ currentage + ageonset + log(time) + status + (1|indID), data = data, relmat = list(indID = kinship_mat))
  #summary(model_test)
  X <- model.matrix(~ currentage + ageonset + log(time) + status, data = data) # imputation model design matrix
  
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
  #conditional_variances <- parSapply(cl, 1:nrow(data), cond_var)
  #conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))
  
  #stopCluster(cl)
  
  ## Imputation steps
  #data <- data |>
  #  mutate(cond_var = conditional_variances_temp)
  
  newx <- data$newx
  Y_obs <- data$newx
  data_imp <- list()
  
  for (m in 1:M) {
    
    ## Step 4
    w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
    
    ## Step 5
    betastar <- betas + w_1 %*% chol(V)
    
    ## Step 6
    mu_star <- X %*% as.vector(betastar)
    
    ## Step 7
    #num_cores <- detectCores() - 2 # 6 cores
    #cl <- makeCluster(num_cores)
    #clusterExport(cl, varlist = c("Sigma", "newx", "mu_star")) 
    #clusterEvalQ(cl, library(Matrix))
    
    E_cond <- function(i) {
      y_minus_i <- newx[-i]
      mu_star_minus_i <- mu_star[-i]
      non_NA <- which(!is.na(y_minus_i))
      y_minus_i <- y_minus_i[non_NA]
      mu_star_minus_i <- mu_star_minus_i[non_NA]
      
      conditional_Expect <- mu_star[i] + Sigma[i,-i][non_NA] %*% solve(Sigma[-i,-i][non_NA, non_NA]) %*% (y_minus_i - mu_star_minus_i)
      return(conditional_Expect)
    }
    conditional_expectations <- sapply(1:nrow(data), E_cond)
    conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
    
    #conditional_vars <- data$cond_var
    
    ## Step 8
    w2i <- rnorm(n = nrow(data), mean = 0, sd = 1)
    
    ## Step 9
    newx_star <- conditional_expectations 
    if (option == "PMM") newx_star <- sapply(newx_star, find_closest, Y_obs) # Comment out for regular MI
    #newx_star <- conditional_expectations + conditional_vars
    data_imp[[m]] <- data |>
      mutate(newx_I = newx_star) |>
      mutate(newx_I = ifelse(!is.na(newx), newx, newx_I))
    
    #stopCluster(cl)
  }
  
  ####################################################################################
  ################################# Analysis Step ####################################
  ####################################################################################
  gamma_results <- list()
  for (i in 1:M) {
    gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                                   gvar = "mgene", design = "pop", base.dist = "Weibull", 
                                   frailty.dist = frailty.dist, 
                                   agemin = min(data_imp[[i]]$currentage[data_imp[[i]]$status == 1]), 
                                   data = data_imp[[i]], parms = c(0.016,3,0, 0, 0,2)) 
  }
  
  model_list <- gamma_results
  pooled_est <- Pooling(model = model_list, imputed_data = data_imp)
  return(pooled_est)
    }, error = function(e) {
      message("Error occurred: ", e$message, "on attempt ", attempts)
      if (attempts == max_attempts) {
        message("Max attempts reached, abandoning this dataset")
        return(NULL)
      }
    })
  }
}
#test_imp <- MI_FamEvent_K_logT(data = data)

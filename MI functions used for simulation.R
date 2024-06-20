############################################################
############# MI Functions used for simulation #############
############################################################

####################################################################################
####################################################################################
####################################################################################
################# Multiple Imputations WITHOUT kinship using log(T) ################
####################################################################################
####################################################################################
####################################################################################
MI_FamEvent_noK_logT <- function(data, M = 5, option = "General", frailty.dist = "gamma") {
  attempts <- 0
  max_attempts <- 3
  success <- FALSE
  while (attempts < max_attempts && !success) {
    attempts <- attempts + 1
    tryCatch({
  ## Step 1 - Empirical estimates
  imp_model <- lm(newx ~ ageonset + log(time) + status , data = data) 
  X <- model.matrix( ~ ageonset + log(time) + status, data = data) 
  betas <- coef(imp_model)
  sigmahat <- sigma(imp_model)
  V <- vcov(imp_model)
  SSE <- sum((imp_model$residuals)^2)
  Y_obs <- data$newx
  
  data_imp <- list()
  for (i in 1:M) {
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
    if (option == "PMM") newx_Imp <- sapply(newx_Imp, find_closest, Y_obs)  # Comment out for regular MI
    
    ## Step 7
    data_imp[[i]] <- data |>
      dplyr::mutate(newx_Imp_temp = newx_Imp) |>
      dplyr::mutate(newx_I = ifelse(is.na(newx), newx_Imp_temp, newx))
    
    #data_imp[[i]] <- data |>
    #  mutate(newx_Imp = newx_Imp) |>
    #  mutate(newx_I = ifelse(is.na(newx), newx_Imp, newx))
  }
  #return(data_imp)
  for (i in 1:M) {
    gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                                   gvar = "mgene", design = "pop", base.dist = "Weibull", 
                                   frailty.dist = frailty.dist, 
                                   agemin = min(data_imp[[i]]$currentage[data_imp[[i]]$status == 1]), 
                                   data = data_imp[[i]], parms = c(0.035,2.3,1, 3, 3,2)) 
  }
  model_list <- gamma_results
  #est_list <- lapply(model_list, function(x) summary(x)$estimates[,1]) 
  #est_aftermi <- colMeans(do.call(rbind, est_list))
  pooled_est <- Pooling(model = model_list, imputed_data = data_imp)
  return(pooled_est)
    }, error = function(e) {
      message("Error occurred: ", e$message, "on attempt ", attempts)
      if (attempts == max_attempts) {
        message("Max attempts reached, abandoning this dataset")
        return(NULL)
      }
    }, warning = function(w) {
      message("Warning occurred: ", w$message, "on attempt ", attempts)
      if (attempts == max_attempts) {
        message("Max attempts reached, abandoning this dataset")
        return(NULL)
      }
    })
  }
}
#MI_FamEvent_noK_logT(data = data)
#MI_function_test <- MI_FamEvent_noK_logT(data = data, option = "PMM")

####################################################################################
####################################################################################
####################################################################################
################# Multiple Imputations WITH kinship using log(T) ###################
####################################################################################
####################################################################################
####################################################################################

####################################################################################
####################################################################################
####################################################################################
################# Multiple Imputations WITH kinship using log(H_0(T)) ##############
####################################################################################
####################################################################################
####################################################################################
MI_FamEvent_K_H0T <- function(data, M = 5, option = "General") {
  attempts <- 0
  max_attempts <- 3
  success <- FALSE
  while(attempts < max_attempts && !success) {
    attempts <- attempts + 1
  tryCatch({
  data_imp <- list() 
  ## Step 1 - Kinship matrix
  kinship_mat <- with(data, kinship(id = indID, dadid = fatherID, momid = motherID,
                                           sex = gender))
  kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
  kinship_mat_sparse <- 2 * kinship_mat_sparse
  kinship_mat <- as.matrix(kinship_mat_sparse)
  
  Iden_mat <- diag(nrow = nrow(data)) # Identity matrix n*n
  Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)
  
  ## Step X - Initial baseline cumulative hazard
  miss50_gamma_cca <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                               design = "pop", base.dist = "Weibull", frailty.dist = frailty.dist, 
                               agemin = min(data$currentage[data$status == 1]), data = data,
                               parms = c(0.035,2.3,1, 3, 3,2))
  baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
  logalpha <- baseline_gammafr[1]
  loglambda <- baseline_gammafr[2]
  data <- data |>
    mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
  Y_obs <- data$newx
  
  #start_time <- Sys.time() # Starting time 
  for (m in 1:M) {
    ## Step 2 - empirical estimates
    model_test <- relmatLmer(newx ~ gender + ageonset + currentage + log(H0) + status + mgene + (1|indID), data = data, relmat = list(indID = kinship_mat))
    #summary(model_test)
    X <- model.matrix(~  gender + ageonset + currentage + log(H0) + status + mgene, data = data) # imputation model design matrix
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
    #conditional_variances <- parSapply(cl, 1:nrow(data), cond_var)
    #conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))
    
    #stopCluster(cl)
    
    ## Imputation steps
    #data <- data |>
    #  mutate(cond_var = conditional_variances_temp)
    
    newx <- data$newx
    
    
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
      Sigma_nonNA <- Sigma[-i,-i][non_NA, non_NA]
      Sigmahat <- Sigma[i,-i][non_NA]
      
      conditional_Expect <- mu_star[i] + Sigmahat %*% solve(Sigma_nonNA) %*% (y_minus_i - mu_star_minus_i)
      return(conditional_Expect)
    }
    conditional_expectations <- sapply(1:nrow(data), E_cond)
    conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
    
    #conditional_vars <- data$cond_var
    #stopCluster(cl)
    
    ## Step 8
    #w2i <- rnorm(n = nrow(data), mean = 0, sd = 1)
    
    ## Step 9
    newx_star <- conditional_expectations 
    if (option == "PMM") newx_star <- sapply(newx_star, find_closest, Y_obs) # Comment out for regular MI
    data$newx_star <- newx_star
    data <- data |>
      mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda)) ) |>
      mutate(newx_I = ifelse(!is.na(newx), newx, newx_star)) 
      #group_by(famID) |>
      #mutate(newx_I = ifelse(!is.na(newx), find_closest(newx_I, newx), newx_I)) |>
      #ungroup() 
    
    data_imp[[m]] <- data
    
    if (m <= 2) {
    ## Step X - update baseline cumulative hazard
    updates_impdata <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
                                gvar = "mgene", 
                                design = "pop", base.dist = "Weibull", frailty.dist = frailty.dist, 
                                agemin = min(data$currentage[data$status == 1]), 
                                data = data,
                                parms = c(0.035,2.3,1, 3, 3,2)) 
    baseline_gammafr_updates <- as.vector(summary(updates_impdata)$estimates[1:2,1])
    logalpha <- baseline_gammafr_updates[1]
    loglambda <- baseline_gammafr_updates[2]
    }
    
  }
  #end_time <- Sys.time() # End time
  #computational_time <- end_time - start_time
  
  ####################################################################################
  ################################# Analysis Step ####################################
  ####################################################################################
  ## Analysis Gamma
  
  gamma_results <- list()
  for (i in 1:M) {
    gamma_results[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                                   design = "pop", base.dist = "Weibull", frailty.dist = frailty.dist, 
                                   agemin = min(data_imp[[i]]$currentage[data_imp[[i]]$status == 1]), data = data_imp[[i]],
                                   parms = c(0.035,2.3,1, 3, 3,2)) 
  }
  model_list <- gamma_results
  #model_list <- lapply(gamma_results, function(x) summary(x)$estimates[,1]) 
  #est_gamma_aftermi <- colMeans(do.call(rbind, est_list))
  pooled_est <- Pooling(model = model_list, imputed_data = data_imp)
  return(pooled_est)
  #return(model_list)
  }, error = function(e) {
    message("Error occurred: ", e$message, "on attempt ", attempts)
    if (attempts == max_attempts) {
      message("Max attempts reached, abandoning this dataset")
      return(NULL)
    }
  })
  }
}
#MI_function_test <- MI_FamEvent_K_H0T(data = data)

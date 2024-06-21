####################################################################################
####################################################################################
####################################################################################
################# Multiple Imputations WITHOUT kinship using log(T) ################
####################################################################################
####################################################################################
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
      imp_model <- lm(newx ~ log(log(ageonset)) + mgene + poly(log(H0), 5) + gender + ageonset + naff + majorgene + generation + relation , data = data) 
      X <- model.matrix( ~ log(log(ageonset)) + mgene + poly(log(H0), 5) + gender + ageonset + naff + majorgene + generation + relation, data = data) 
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
        gamma_results[[i]] <- penmodel(survival::Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
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
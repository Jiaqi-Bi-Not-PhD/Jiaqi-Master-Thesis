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
MI_FamEvent_noK_logT <- function(data, M = 5) {
  ## Step 1 - Empirical estimates
  imp_model <- lm(newx ~ ageonset + log(time) + status , data = data) 
  X <- model.matrix( ~ ageonset + log(time) + status, data = data)
  summary(imp_model)
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
    #newx_Imp <- sapply(newx_Imp, find_closest, Y_obs)  # Comment out for regular MI
    
    ## Step 7
    data_imp[[i]] <- data |>
      mutate(newx_Imp_temp = newx_Imp) |>
      mutate(newx_I = ifelse(is.na(newx), newx_Imp_temp, newx))
    
    #data_imp[[i]] <- data |>
    #  mutate(newx_Imp = newx_Imp) |>
    #  mutate(newx_I = ifelse(is.na(newx), newx_Imp, newx))
  }
}

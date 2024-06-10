####################################################################################
############# Multiple Imputations WITHOUT kinship using log(T) - draft ############
####################################################################################
################################# Imputation Step ##################################
####################################################################################

## Step 1 - Empirical estimates
imp_model <- lm(newx ~ gender + currentage + log(time) + status + mgene , data = miss50_famx) 
X <- model.matrix( ~ gender + currentage + log(time) + status + mgene , data = miss50_famx)
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
  newx_Imp <- sapply(newx_Imp, find_closest, Y_obs)
  
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

num_cores <- detectCores() - 2 # 6 cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

gamma_results <- foreach(i = 1:5, .packages = c("survival")) %dopar% {
  penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
           gvar = "mgene", design = "pop", base.dist = "Weibull", 
           frailty.dist = "gamma", agemin = 20, 
           data = miss50_famx_imp[[i]], parms = c(0.016,3,0, 0, 0,2)) 
}
MI_noK_logT <- gamma_results 
est_list <- lapply(MI_noK_logT, function(x) summary(x)$estimates[,1]) 
est_gamma_aftermi <- colMeans(do.call(rbind, est_list))
stopCluster(cl)


#################### Check Zone ####################
check_model <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", 
         gvar = "mgene", design = "pop", base.dist = "Weibull", 
         frailty.dist = "gamma", agemin = 20, 
         data = miss50_famx_imp[[3]], parms = c(0.45,1,1,1,1, 1)) 
summary(check_model)

miss50_famx_imp[[5]] |> ggplot(aes(x = newx_I)) + geom_density()



######################################################
################## Simulation 1 ######################
######################################################
####### Missingness generated from mice package ######
######################################################

## Generate kinship data - Number of family = 10
install.packages("FamEvent_3.1.tar.gz", type = "source")
library(FamEvent)
library(mice)

set.seed(123)
famx <- simfam(N.fam = 20, design = "pop", variation = "frailty", 
               base.dist = "Weibull", frailty.dist = "gamma", interaction = FALSE,
               add.x = TRUE, x.dist = "normal", x.parms = c(0, 1),  depend = 2, 
               base.parms = c(0.016,3), vbeta = c(1, 3, 3)) # Number of family = 10
IBDmatrix <- diag(1, dim(famx)[1])
kinship_mat <- with(famx, kinship(id = indID, dadid = fatherID, momid = motherID,
                                         sex = gender))
kinship_mat <- 2*kinship_mat
inputdata <- famx[, !names(famx) %in%c("time", "status", "ageonset", "fsize", "naff")]

fam2 <- simfam2(design = "pop", variation = c("kinship", "IBD"),
                depend = c(2, 2), base.dist="Weibull", base.parms =c(0.016, 3),
                var_names = c("gender", "mgene", "newx"), vbeta = c(1, 3, 3),
                agemin=20, inputdata=inputdata, IBD=kinship_mat) # covariates are gender, mgene, newx

fam2 |>
  ggplot(aes(x = newx)) + 
  geom_density() # Check the normality

fam2 <- fam2 |>
  mutate(I_Tp_j.ap_j = ifelse(proband == 1 & time < currentage, 1, 0))

## Gamma frailty optimization test (No missing data)
initial_params <- c(1/41.41327,1,0,0,1)
X <- as.matrix(fam2[,c("mgene", "newx")], ncol = 2)
Y <- as.matrix(fam2[,c("time", "status")], ncol = 2)
nomiss_gamma <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                        data = fam2, X = X, Y = Y, nbase = 2,
                        design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                        agemin = 20, 
                        control = list(maxit = 2000), hessian = TRUE) # True parameters = -4.48759935 -0.41795578  1.46231306  0.66128631 -0.05276222

nomiss_gamma_hessian <- nomiss_gamma$hessian
nomiss_gamma_hessian <- as.matrix(nearPD(nomiss_gamma_hessian)$mat)
vcov_mat <- solve(nomiss_gamma_hessian)
param_est <- nomiss_gamma$par
standard_errors <- sqrt(diag(vcov_mat))
z_scores <- param_est/standard_errors
p_values <- 2*(pnorm(abs(z_scores), lower.tail = FALSE))
final_results <- list( Estimates = param_est, 
                       Std.Error = standard_errors, 
                       Z = z_scores, 
                       p_val = p_values) # p-value

################# Suppose 20% missing - Missing at Random ####################
################# Suppose a reasonable weights ###############################
miss_pattern <- matrix(c(rep(1, 11), 0, rep(1, 8)), nrow = 1, byrow = TRUE)
ampute_test <- ampute(data = fam2, prop = 0.2, patterns = miss_pattern)
md.pattern(ampute_test$amp, rotate.names = TRUE) # Missing pattern
ampute_test$mech
ampute_test$weights
reasonable_weights <- ampute_test$weights
reasonable_weights[1,] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0)
ampute_test <- ampute(data = fam2, prop = 0.4, patterns = miss_pattern, 
                      mech = "MAR", weights = reasonable_weights) # Missing proportions
miss20_fam2 <- ampute_test$amp
## CCA 
initial_params <- c(1/41.41327,1,0, 0, 1)
X <- as.matrix(miss20_fam2[,c("mgene", "newx")], ncol = 2)
Y <- as.matrix(miss20_fam2[,c("time", "status")], ncol = 2)
miss20_gamma_CCA <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                        data = miss20_fam2, X = X, Y = Y, nbase = 2,
                        design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                        agemin = 20, 
                        control = list(maxit = 2000), hessian = TRUE) # CCA parameters = -3.3481794 -0.4696807  0.6456010  0.1215593  2.2131349

miss20_gamma_hessian <- miss20_gamma_CCA$hessian
miss20_gamma_hessian <- as.matrix(nearPD(miss20_gamma_hessian)$mat)
vcov_mat <- solve(miss20_gamma_hessian)
param_est <- miss20_gamma_CCA$par
standard_errors <- sqrt(diag(vcov_mat))
z_scores <- param_est/standard_errors
p_values <- 2*(pnorm(abs(z_scores), lower.tail = FALSE))
final_results <- list( Estimates = param_est, 
                       Std.Error = standard_errors, 
                       Z = z_scores, 
                       p_val = p_values) # p-value

## MI
## Step 1 - Kinship matrix
kinship_mat <- with(miss20_fam2, kinship(id = indID, dadid = fatherID, momid = motherID,
                                       sex = gender))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(miss20_fam2)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step 2 - empirical estimates
model_test <- relmatLmer(newx ~ gender + proband + mgene + log(time) + status + (1|indID), data = miss20_fam2, relmat = list(indID = kinship_mat))
summary(model_test)
X <- model.matrix( ~ gender + proband + mgene + log(time) + status, data = miss20_fam2) 

#betas <- coef(model_test)
betas <- c(1.96808, -0.04012, 0.05765, -0.25335, -0.61984, 1.22873)
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
conditional_variances <- parSapply(cl, 1:nrow(miss20_fam2), cond_var)
conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))

stopCluster(cl)

## Imputation steps
miss20_fam2 <- miss20_fam2 |>
  mutate(cond_var = conditional_variances_temp)

newx <- miss20_fam2$newx


fam2_newx_imp_fam <- list() 
for (m in 1:20) {
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
  conditional_expectations <- parSapply(cl, 1:nrow(miss20_fam2), E_cond)
  conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
  
  conditional_vars <- miss20_fam2$cond_var
  
  ## Step 8
  w2i <- rnorm(n = nrow(miss20_fam2), mean = 0, sd = 1)
  
  ## Step 9
  newx_star <- conditional_expectations + w2i * conditional_vars
  fam2_newx_imp_fam[[m]] <- miss20_fam2 |>
    mutate(newx_I = newx_star) |>
    mutate(newx_I = ifelse(!is.na(newx), newx, newx_I))
  
  stopCluster(cl)
}



## Analysis Gamma
gamma_results <- list()
for (i in 1:20) {
  initial_params <- c(1/41.41327,1, 0, 0, 1)
  X <- as.matrix(fam2_newx_imp_fam[[i]][,c("mgene", "newx_I")], ncol = 2)
  Y <- as.matrix(fam2_newx_imp_fam[[i]][,c("time", "status")], ncol = 2)
  gamma_forgraph <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                          data = fam2_newx_imp_fam[[i]], X = X, Y = Y, nbase = 2,
                          design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                          agemin = 20, 
                          control = list(maxit = 2000))
  gamma_results[[i]] <- gamma_forgraph$par
}

results1 <- colMeans(do.call(rbind, gamma_results))
results_coxme <- coxme(Surv(time, status) ~ mgene + newx , varlist = )


###############################################################################################
########### Multiple Imputation without Considering Family Structure ##########################
###############################################################################################

## Step 1 - Empirical estimates
imp_model <- lm(newx ~ gender + proband + mgene + log(time) + status, data = miss20_fam2) 
summary(imp_model)
betas <- coef(imp_model)
sigmahat <- sigma(imp_model)
V <- vcov(imp_model)
SSE <- sum((imp_model$residuals)^2)

fam2_imp <- list()
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
  fam2_imp[[i]] <- miss20_fam2 |>
    mutate(newx_I = ifelse(is.na(newx), betastar[,1] + betastar[,2]*gender + betastar[,3]*proband + betastar[,4]*mgene + betastar[,5]*log(time) + betastar[,6]*status + rnorm(n = 1, mean = 0, sd = 1)*sigmastar, newx))
}

gamma_results_nofamstruc <- list()
for (i in 1:20) {
  initial_params <- c(1/41.41327,1, 0, 0, 1)
  X <- as.matrix(fam2_imp[[i]][,c("mgene", "newx_I")], ncol = 2)
  Y <- as.matrix(fam2_imp[[i]][,c("time", "status")], ncol = 2)
  gamma_forgraph <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                          data = fam2_imp[[i]], X = X, Y = Y, nbase = 2,
                          design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                          agemin = 20, 
                          control = list(maxit = 2000))
  gamma_results_nofamstruc[[i]] <- gamma_forgraph$par
}

results2 <- colMeans(do.call(rbind, gamma_results_nofamstruc))
nomiss_gamma$par
miss20_gamma_CCA$par
results1

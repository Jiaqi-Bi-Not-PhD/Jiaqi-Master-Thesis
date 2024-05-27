#### Multiple Imputation using log(T) ####
## MI
baseline_gammafr <- as.vector(summary(miss50_gamma)$estimates[1:2,1])
famx_newx_imp_fam <- list() 

baseline_gammafr <- as.vector(summary(miss50_gamma)$estimates[1:2,1])
## Step 1 - Kinship matrix
kinship_mat <- with(miss50_famx, kinship(id = indID, dadid = fatherID, momid = motherID,
                                         sex = gender))
kinship_mat_sparse <- Matrix(kinship_mat, sparse = TRUE)
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

Iden_mat <- diag(nrow = nrow(miss50_famx)) # Identity matrix n*n
Iden_mat_sparse <- Matrix(Iden_mat, sparse = TRUE)

## Step 2 - empirical estimates
model_test <- relmatLmer(newx ~ log(time) * status  + (1|indID), data = miss50_famx, relmat = list(indID = kinship_mat))
#summary(model_test)
X <- model.matrix(model_test) # imputation model design matrix

#betas <- coef(model_test)
betas <- as.vector(summary(model_test)$coefficients[,1]) # beta coefficients
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
conditional_variances <- parSapply(cl, 1:nrow(miss50_famx), cond_var)
conditional_variances_temp <- as.vector(do.call(rbind, conditional_variances))

stopCluster(cl)

## Imputation steps
miss50_famx <- miss50_famx |>
  mutate(cond_var = conditional_variances_temp)

newx <- miss50_famx$newx

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
  conditional_expectations <- parSapply(cl, 1:nrow(miss50_famx), E_cond)
  conditional_expectations <- as.vector(do.call(rbind, conditional_expectations))
  
  conditional_vars <- miss50_famx$cond_var
  
  ## Step 8
  w2i <- rnorm(n = nrow(miss50_famx), mean = 0, sd = 1)
  
  ## Step 9
  newx_star <- conditional_expectations + w2i * sqrt(conditional_vars)
  #newx_star <- conditional_expectations + conditional_vars
  famx_newx_imp_fam[[m]] <- miss50_famx |>
    mutate(newx_I = newx_star) |>
    mutate(newx_I = ifelse(!is.na(newx), newx, newx_I))
  
  stopCluster(cl)
}
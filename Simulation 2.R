######################################################
################## Simulation 2 ######################
######################################################
####### Missingness generated from logistic reg ######
######################################################

## Generate kinship data - Number of family = 80
install.packages("FamEvent_3.1.tar.gz", type = "source")
library(FamEvent)
library(mice)

set.seed(123123)
famx <- simfam(N.fam = 200, design = "pop", variation = "frailty", 
               base.dist = "Weibull", frailty.dist = "gamma", interaction = FALSE,
               add.x = TRUE, x.dist = "normal", x.parms = c(0, 1),  depend = 2, 
               base.parms = c(0.016,3), vbeta = c(1.5, 1, 1.5)) # Number of family = 100
IBDmatrix <- diag(1, dim(famx)[1])
kinship_mat <- with(famx, kinship(id = indID, dadid = fatherID, momid = motherID,
                                  sex = gender))
kinship_mat <- 2*kinship_mat
inputdata <- famx[, !names(famx) %in%c("time", "status", "ageonset", "fsize", "naff")]

fam2 <- simfam2(design = "pop", variation = "kinship", 
                depend = 2, base.dist="Weibull", base.parms =c(0.016, 3),
                var_names = c("gender", "mgene", "newx"), vbeta = c(1.5, 1, 1.5),
                agemin=20, inputdata=inputdata, IBD=kinship_mat) # covariates are gender, mgene, newx

fam2 |>
  ggplot(aes(x = newx)) + 
  geom_density() # Check the normality

## Gamma frailty optimization test (No missing data)
initial_params <- c(1/41.41327,1,0,0, 0, 1)
X <- as.matrix(fam2[,c("gender", "mgene", "newx")], ncol = 3)
Y <- as.matrix(fam2[,c("time", "status")], ncol = 2)
nomiss_gamma <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                      data = fam2, X = X, Y = Y, nbase = 3,
                      design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                      agemin = 20, 
                      control = list(maxit = 2000)) # True parameters = -6.712394e-01 -1.498446e+00  2.155013e+00 -1.027516e+00  2.522497e+00  3.845956e-18

################# Suppose 20% missing - Missing at Random ####################
################# Suppose a reasonable weights ###############################
set.seed(123)
simulate_mar <- function(data, target_variable, missing_rate = 0.2) {
  data_mar <- data
  n <- nrow(data_mar)
  formula <- as.formula("as.numeric(runif(n) < 0.5) ~ mgene + gender + proband + log(time)*status")
  logit_model <- glm(formula, data = data_mar, family = binomial)
  
  probabilities <- predict(logit_model, type = "response")
  num_missing <- round(missing_rate * n)
  missing_indices <- sample(1:n, num_missing, prob = probabilities, replace = FALSE)
  
  data_mar[missing_indices, target_variable] <- NA
  return(data_mar)
}
miss20_fam2 <- simulate_mar(data = fam2, "newx")
sum(is.na(miss20_fam2$newx))/nrow(miss20_fam2)

## CCA 
initial_params <- c(1/41.41327,1,0,0, 0, 1)
X <- as.matrix(miss20_fam2[,c("gender", "mgene", "newx")], ncol = 3)
Y <- as.matrix(miss20_fam2[,c("time", "status")], ncol = 2)
miss20_gamma_CCA <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                          data = miss20_fam2, X = X, Y = Y, nbase = 3,
                          design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                          agemin = 20, 
                          control = list(maxit = 2000)) # CCA parameters = -6.017403e-01 -1.638124e+00  2.634569e+00 -1.029436e+00  2.432744e+00  3.857695e-18
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
model_test <- relmatLmer(newx ~ gender + proband + mgene + log(time) * status + (1|indID), data = miss20_fam2, relmat = list(indID = kinship_mat))
summary(model_test)
X <- model.matrix( ~ gender + proband + mgene + log(time) * status, data = miss20_fam2) 

#betas <- coef(model_test)
betas <- c(0.90222, -0.10804, 0.23882, -0.02422, -0.30160, 4.64159, -1.08143)
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

num_cores <- detectCores() - 2 # 6 cores
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("Sigma", "newx", "mu_star")) 
clusterEvalQ(cl, library(Matrix))


fam2_newx_imp_fam <- list() 
for (m in 1:20) {
  ## Step 4
  w_1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  
  ## Step 5
  betastar <- betas + w_1 %*% chol(V)
  
  ## Step 6
  mu_star <- X %*% as.vector(betastar)
  
  ## Step 7
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
}

stopCluster(cl)

## Analysis Gamma
gamma_results <- list()
for (i in 1:20) {
  initial_params <- c(1/41.41327,1,0,0, 0, 1)
  X <- as.matrix(fam2_newx_imp_fam[[i]][,c("gender", "mgene", "newx_I")], ncol = 2)
  Y <- as.matrix(fam2_newx_imp_fam[[i]][,c("time", "status")], ncol = 2)
  gamma_forgraph <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                          data = fam2_newx_imp_fam[[i]], X = X, Y = Y, nbase = 2,
                          design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                          agemin = 20, 
                          control = list(maxit = 2000))
  gamma_results[[i]] <- gamma_forgraph$par
}

colMeans(do.call(rbind, gamma_results))

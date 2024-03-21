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

## Multiple imputation 
for (m in 1:20) {
  brca1_prs_MI_completed <- brca1_prs_MI
  for (i in missing_indices) {
    indID <- brca1_prs_MI_completed$indID[brca1_prs_MI_completed$indID == i]
    mu_PRS <- brca1_prs_MI_completed$mu_PRS[brca1_prs_MI_completed$indID == i]
    Sigma_index <- match(indID, kinship_ids)
    Sigma_i <- Sigma[Sigma_index, Sigma_index]
  
    imputed_value <- mvrnorm(n = 1, mu = mu_PRS, Sigma = matrix(Sigma_i, nrow = 1)) # Make only 1 draw?
    brca1_prs_MI_completed$PRS_I[brca1_prs_MI_completed$indID == i] <- imputed_value
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

check <- imputed_data_list[[2]]

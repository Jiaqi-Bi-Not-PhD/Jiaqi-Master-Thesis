##########################################################
########### Multiple Imputation ##########################
##########################################################
## Kinship matrix
library(kinship2)
library(Matrix)
library(lme4qtl)
library(lmerTest)
library(parallel)

kinship_mat <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID,
                                       sex = rep(2,nrow(brca1_prs))))
kinship_mat_sparse <- as(kinship_mat, "sparseMatrix")
kinship_mat_sparse <- 2 * kinship_mat_sparse
kinship_mat <- as.matrix(kinship_mat_sparse)

## Linear mixed effect model using kinship matrix
model_test <- relmatLmer(PRS ~ proband + mgeneI + currentage + timeBC * BC + (1|indID), 
                         data = brca1_prs, relmat = list(indID = kinship_mat))
summary(model_test)


## Constructing kinship and covariance
sigma_g <- attr(VarCorr(model_test)$indID, "stddev")
sigma_e <- attr(VarCorr(model_test), "sc")
K <- sigma_g^2 * as(kinship_mat, "sparseMatrix") + sigma_e^2
K <- as.matrix(K)

## Imputation step
numCores <- detectCores() - 2

impute_PRS_dataset <- function(dataset, mean_PRS, Sigma, seed) {
  set.seed(123)  # Ensure reproducibility
  
  ## PRS has NAs
  ## Logic

  return(imputed_dataset)
}
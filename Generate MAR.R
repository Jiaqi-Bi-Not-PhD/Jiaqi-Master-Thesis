library(parallel)
library(dplyr)
library(survival)
library(parallelly)
library(FamEvent)
library(mice)

source("cumhaz.R")
source("dlaplace.R")
source("find_closest.R")
source("gh_REVISED; May 27.R")
source("hazards.R")
source("hermite.R")
source("laplace.R")
source("loglik_frailty_REVISED; May 27.R")
source("penmodel_REVISED; May 27.R")
source("Bias, MSE, Coverage Functions.R")
source("Rubins Rule.R")


## True thetas
true_beta1 <- 1
true_beta2 <- 3
true_beta3 <- 3
true_kappa <- log(2)
true_alpha <- log(0.035)
true_lambda <- log(2.3)

true_value2 <- c(0.035, 2.3, 1, 3, 3, 0.6931472)
true_value <- c(-3.3524072, 0.8329091, 1.0000000, 3.0000000, 3.0000000, 0.6931472)

n_simulations <- 1000
n_cores <- parallelly::availableCores()

#######################################################################################
################################ Generate MAR #########################################
#######################################################################################
nomiss_datasets_test <- readRDS("1000_NoMissing_Gamma_LogNormal_50and500.RData")
names(nomiss_datasets_test) <- c("simulated_dataset_1000_Complete_lognormal_50_fam", "simulated_dataset_1000_Complete_lognormal_500_fam", "simulated_dataset_1000_Complete_gamma_50_fam", "simulated_dataset_1000_Complete_gamma_500_fam")
simulated_dataset_1000_Complete_gamma_50_fam <- nomiss_datasets_test$simulated_dataset_1000_Complete_gamma_50_fam
simulated_dataset_1000_Complete_gamma_500_fam <- nomiss_datasets_test$simulated_dataset_1000_Complete_gamma_500_fam
simulated_dataset_1000_Complete_lognormal_50_fam <- nomiss_datasets_test$simulated_dataset_1000_Complete_lognormal_50_fam
simulated_dataset_1000_Complete_lognormal_500_fam <- nomiss_datasets_test$simulated_dataset_1000_Complete_lognormal_500_fam

proportions <- c(0.20, 0.40, 0.60, 0.80)

miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)

#parallel::clusterEvalQ(cl, library(mice))
## Missingness & proportion function
introduce_missingness <- function(index, dataset_list, prop, miss_pattern) {
  famx <- dataset_list[[index]]
  
  ampute_test <- mice::ampute(data = famx, prop = prop, patterns = miss_pattern)
  reasonable_weights <- ampute_test$weights
  reasonable_weights[1,] <- c(0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 0)
  ampute_test <- mice::ampute(data = famx, prop = prop, patterns = miss_pattern, 
                              mech = "MAR", weights = reasonable_weights)
  miss_famx <- ampute_test$amp
  
  return(miss_famx)
}

## Function to run the missingness introduction for all proportions
run_for_all_proportions <- function(i, dataset_list, proportions, miss_pattern) {
  results <- lapply(proportions, function(prop) {
    introduce_missingness(i, dataset_list, prop, miss_pattern)
  })
  names(results) <- as.character(proportions)
  return(results)
}


#n_cores <- parallelly::availableCores()
#cl <- parallel::makeCluster(n_cores)


## List of complete datasets
complete_datasets <- list(
  gamma_50_fam = simulated_dataset_1000_Complete_gamma_50_fam,
  gamma_500_fam = simulated_dataset_1000_Complete_gamma_500_fam,
  lognormal_50_fam = simulated_dataset_1000_Complete_lognormal_50_fam,
  lognormal_500_fam = simulated_dataset_1000_Complete_lognormal_500_fam
)

#parallel::clusterExport(cl, c("introduce_missingness", "run_for_all_proportions", "miss_pattern", "proportions",
#                    "complete_datasets"))

mar_datasets <- list()

## Iterate over each complete dataset and apply the missing data mechanism
for (dataset_name in names(complete_datasets)) {
  dataset <- complete_datasets[[dataset_name]]
  # parallel::clusterExport(cl, "dataset")
  ## Apply the function in parallel
  results_parallel <- parallel::mclapply(1:n_simulations, function(i) {
    run_for_all_proportions(i, dataset, proportions, miss_pattern)
  }, mc.cores=n_cores)
  
  ## Store the results in the mar_datasets list
  for (prop in proportions) {
    mar_datasets[[paste0("simulated_dataset_1000_MAR", prop * 100, "_", dataset_name)]] <- parallel::mclapply(1:n_simulations, function(i) results_parallel[[i]][[as.character(prop)]], mc.cores = n_cores)
  }
}

#stopCluster(cl)

### mar_datasets contains 24 lists of 1000 datasets -> 2400 datasets
### nomiss_datasets contains 6 lists of 1000 datasets -> 600 datasets
### mar_datasets
#names(mar_datasets)
#saveRDS(complete_datasets, file = "100 complete datasets lists.RData")
saveRDS(mar_datasets, file = "1000_mar_lists.RData")
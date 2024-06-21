library(parallelly)
library(kinship2)
library(survival)
library(FamEvent)
library(Matrix)
library(lme4)
library(parallel)
library(mice)
library(dplyr)
library(MASS)
library(lme4qtl)
#MI_FamEvent_K_H0T
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
source("MI function; K log(H0(T)).R") # MI_FamEvent_K_H0T
source("MI function; No K log(H0(T)).R") # MI_FamEvent_noK_H0T
source("MI function; No K Log(T).R") # MI_FamEvent_noK_logT
source("MI function; with K log(T).R") # MI_FamEvent_K_logT
source("Rubins Rule.R")

simulated_dataset_1000_MAR20_gamma_50_fam <- readRDS("simulated_dataset_1000_MAR20_gamma_50_fam.RData")
true_value <- c(-3.3524072, 0.8329091, 1.0000000, 3.0000000, 3.0000000, 0.6931472)
n_cores <- parallelly::availableCores()

## Kinship + H0
process_datasets_K_H0T <- function(datasets, true_value) {
  results <- parallel::mclapply(datasets, function(data) {
    imputed_results <- MI_FamEvent_K_H0T(data, option = "General", M=5, frailty.dist = "gamma") # For MI no K log T Option General, adjust accordingly for other options
    metrics <- calculate_metrics(imputed_results, true_value)
    return(metrics)
  }, mc.cores = n_cores) # Specify cores
  results <- dplyr::bind_rows(results)
  return(results)
}

sim_results <- process_datasets_K_H0T(datasets = simulated_dataset_1000_MAR20_gamma_50_fam, 
                                      true_value = true_value)
sim_results_toprint <- sim_results |> dplyr::group_by(Parameters) |> dplyr::summarise(Bias = mean(Bias, na.rm = TRUE),
                                                                        RMSE = sqrt(mean(MSE, na.rm = TRUE)),
                                                                        MSE = mean(MSE, na.rm = TRUE),
                                                                        Coverage = mean(Coverage, na.rm = TRUE))
scenario <- "simulated_dataset_1000_MAR20_gamma_50_fam"
final_results <- list(scenario, sim_results_toprint)
print(final_results)

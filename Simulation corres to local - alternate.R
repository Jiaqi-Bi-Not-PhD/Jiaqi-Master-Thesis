###########################################################
############## Simulation parallel on local 2 #############
###########################################################
process_datasets <- function(datasets, true_value) {
  results <- mclapply(datasets, function(data) {
    imputed_results <- MI_FamEvent_K_H0T(data, option = "General") # For MI no K log T Option General, adjust accordingly for other options
    metrics <- calculate_metrics(imputed_results, true_value)
    return(metrics)
  }, mc.cores = 6) # Specify cores
  results <- bind_rows(results)
  return(results)
}
start_time <- Sys.time()
scenario_names <- names(simulated_dataset_100_MAR20_Gamma_50)
sim_results <- process_datasets(datasets = simulated_dataset_100_MAR20_Gamma_50$simulated_dataset_100_MAR20_Gamma_50, 
                                true_value = true_value)
sim_results |> group_by(Parameters) |> summarise(Bias = mean(Bias, na.rm = TRUE),
                                                 RMSE = sqrt(mean(MSE, na.rm = TRUE)),
                                                 MSE = mean(MSE, na.rm = TRUE),
                                                 Coverage = mean(Coverage, na.rm = TRUE))
end_time <- Sys.time()
running_time <- end_time - start_time
#names(mar_datasets_gamma)
#mar_datasets_gamma_50_500 <- mar_datasets_gamma[-(5:8)]
#mar_datasets_gamma_test <- mar_datasets_gamma[c(1,12)]
#mar_datasets_gamma_test$simulated_dataset_100_MAR20_Gamma_50 <- mar_datasets_gamma_test$simulated_dataset_100_MAR20_Gamma_50[c(1:2)]
#mar_datasets_gamma_test$simulated_dataset_100_MAR80_Gamma_500 <- mar_datasets_gamma_test$simulated_dataset_100_MAR80_Gamma_500[c(1:2)]
#mar_datasets_gamma_50_500
#subset_test <- simulated_dataset_100_MAR20_Gamma_50$simulated_dataset_100_MAR20_Gamma_50[c(1:2)]
#names(subset_test) <- names(simulated_dataset_100_MAR20_Gamma_50$simulated_dataset_100_MAR20_Gamma_50)[c(1:2)]
start_time <- Sys.time()
scenario_names <- names(simulated_dataset_100_MAR20_Gamma_50) # Gamma frailty, adjust accordingly for lognormal frailty

all_scenarios_results <- mclapply(scenario_names, function(scenario) {
  scenario_data <- simulated_dataset_100_MAR20_Gamma_50[[scenario]] # Gamma frailty, adjust accordingly for lognormal frailty
  scenario_results <- process_datasets(scenario_data, true_value)
  scenario_results <- scenario_results |> mutate(Scenario = as.character(scenario))
  return(scenario_results)
}, mc.cores = 6)

all_scenarios_results <- bind_rows(all_scenarios_results)

summary_results <- all_scenarios_results |>
  group_by(Parameters, Scenario) |>
  summarize(
    Bias = mean(Bias, na.rm = TRUE),
    RMSE = sqrt(mean(MSE, na.rm = TRUE)),
    MSE = mean(MSE, na.rm = TRUE),
    Coverage = mean(Coverage, na.rm = TRUE)
  )
saveRDS(summary_results, file = "test_results.RData")
test_results <- readRDS("test_results.RData")
end_time <- Sys.time()
running_time <- end_time - start_time

###########################################################
############## Simulation parallel on local 2 #############
###########################################################
process_datasets <- function(datasets, true_value) {
  results <- lapply(datasets, function(data) {
    imputed_results <- MI_FamEvent_noK_logT(data) # For MI no K log T Option General, adjust accordingly for other options
    metrics <- calculate_metrics(imputed_results, true_value)
    return(metrics)
  })
  results <- bind_rows(results)
  return(results)
}

#names(mar_datasets_gamma)
#mar_datasets_gamma_50_500 <- mar_datasets_gamma[-(5:8)]
#mar_datasets_gamma_test <- mar_datasets_gamma[c(1,12)]
#mar_datasets_gamma_test$simulated_dataset_100_MAR20_Gamma_50 <- mar_datasets_gamma_test$simulated_dataset_100_MAR20_Gamma_50[c(1:2)]
#mar_datasets_gamma_test$simulated_dataset_100_MAR80_Gamma_500 <- mar_datasets_gamma_test$simulated_dataset_100_MAR80_Gamma_500[c(1:2)]

scenario_names <- names(mar_datasets_gamma_test) # Gamma frailty, adjust accordingly for lognormal frailty

all_scenarios_results <- mclapply(scenario_names, function(scenario) {
  scenario_data <- mar_datasets_gamma_test[[scenario]] # Gamma frailty, adjust accordingly for lognormal frailty
  scenario_results <- process_datasets(scenario_data, true_value)
  scenario_results <- scenario_results |> mutate(Scenario = as.character(scenario))
  return(scenario_results)
}, mc.cores = 6)

all_scenarios_results <- bind_rows(all_scenarios_results)

summary_results <- all_scenarios_results |>
  group_by(Scenario, Parameter) |>
  summarize(
    Bias = mean(Bias),
    MSE = mean(MSE),
    Coverage = mean(Coverage)
  )

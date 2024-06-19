###########################################################
################### Simulation Functions ##################
###########################################################
true_value <- c(true_alpha, true_lambda, true_beta1, true_beta2, true_beta3, true_kappa)

calculate_metrics <- function(imp_results, true_value = true_value) {
  Parameters <- imp_results$Parameters
  parameter_estimates <- imp_results$Estimates
  total_variance <- imp_results$`Total Var.`
  lower_ci <- imp_results$`Lower.CI`
  upper_ci <- imp_results$`Upper.CI`
  imp_results$true_value <- true_value
  True_Value <- imp_results$true_value
  Estimates <- imp_results$Estimates
  
  bias <- parameter_estimates - true_value
  mse <- (parameter_estimates - true_value)^2
  coverage <- (lower_ci <= true_value & upper_ci >= true_value)
  
  return(tibble(
    Parameters = Parameters,
    Estimates = Estimates,
    #True_Value = True_Value,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    Bias = bias,
    MSE = mse,
    Coverage = coverage
  ))
}
calculate_metrics(MI_function_test, true_value)

###### Using simfam ######
library(FamEvent)
surv_data <- simfam(N.fam = 110, base.dist = "Weibull", frailty.dist = "lognormal",
                    base.parms = c(0.5,0.5), vbeta = c(2, 2), variation = "frailty", 
                    depend = 3)
family_events <- surv_data |>
  dplyr::group_by(famID) |>
  dplyr::summarise(df1 = sum(status))
surv_data <- merge(surv_data, family_events, by = "famID")

##### Simulate the log-normal frailty log-likelihood #####
X <- as.matrix(data.frame(surv_data$gender, surv_data$majorgene), 
               nrow=nrow(surv_data), 
               ncol = 2)
Y <- as.matrix(data.frame(surv_data$time, surv_data$status), 
               nrow = nrow(surv_data), 
               ncol = 2)

initial_params <- c(0.2,1, 0, 0, 1)
optim(par = initial_params, fn = lognormal_single,
      data = surv_data, X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
      agemin = 18, control = list(maxit = 2000)) # Log-Normal

optim(par = initial_params, fn = loglik_frailty_single_gamma,
      data = surv_data, X = X, Y = Y, nbase = 2,
      design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
      agemin = 18, 
      control = list(maxit = 2000)) # Gamma

######


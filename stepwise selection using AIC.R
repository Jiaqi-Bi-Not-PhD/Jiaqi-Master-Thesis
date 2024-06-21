miss50_gamma_cca <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                             design = "pop", base.dist = "Weibull", frailty.dist = frailty.dist, 
                             agemin = min(data$currentage[data$status == 1]), data = famx,
                             parms = c(0.035,2.3,1, 3, 3,2))
baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
famx <- famx |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0
full_model <- lm(newx ~ .  + poly(log(H0), 5) + log(log(ageonset)), data = famx)
null_model <- lm(newx ~ 1, data = famx)
library(MASS)
stepwise_model <- stepAIC(null_model, 
                          scope = list(lower = null_model, upper = full_model),
                          direction = "both", 
                          trace = FALSE)

summary(stepwise_model)

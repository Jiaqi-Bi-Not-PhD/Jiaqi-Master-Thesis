######################################################
################## Simulation 1 ######################
######################################################
####### Missingness generated from mice package ######
######################################################
brca1_prs <- brca1_prs |>
  mutate(weight = 1)
brca1_prs_cca <- brca1_prs |>
  filter(!is.na(PRS))
## Generate kinship data - Number of family = 10
install.packages("FamEvent_3.1.tar.gz", type = "source")
library(FamEvent)
library(mice)
library(latex2exp)

set.seed(123123)
famx <- simfam(N.fam = 200, design = "pop", variation = "frailty", 
               base.dist = "Weibull", frailty.dist = "gamma", interaction = FALSE,
               add.x = TRUE, x.dist = "normal", x.parms = c(0, 1),  depend = 2, 
               base.parms = c(0.016,3), vbeta = c(1, 3, 3)) # Number of family = 10

## Test if Gamma frailty model gives the correct estimates ##
test_model <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(test_model)
## Test if Gamma frailty model gives the correct estimates ##

## Estimates when there are no missing data
nomiss_gamma <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(nomiss_gamma)

###################################### Missing Data Plot ########################################
## Initial baseline cumulative hazard
miss50_gamma_cca <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                             design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx,
                             parms = c(1/41.41327,1,0,0,0, 1))
summary(miss50_gamma_cca)
baseline_gammafr <- as.vector(summary(miss50_gamma_cca)$estimates[1:2,1])
logalpha <- baseline_gammafr[1]
loglambda <- baseline_gammafr[2]
miss50_famx <- miss50_famx |>
  mutate(H0 = (exp(logalpha)^exp(loglambda)) * (exp(loglambda)^2) * (time^exp(loglambda))) # Generate H0

miss50_famx |>
  ggplot(aes(x = factor(status), y = newx)) +
  geom_boxplot() 

miss50_famx |>
  ggplot(aes(x = log(time), y = newx)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = TeX('$log(T)$ vs. newx'), 
       x = TeX('$log(T)$'), 
       y = "newx") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

miss50_famx |>
  ggplot(aes(x = H0, y = newx)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = TeX("$H_{0}(T)$ vs. newx"), 
       x = TeX("$H_0(T)$"), 
       y = "newx") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## CCA Gamma
miss50_famx_CCA <- miss50_famx |>
  filter(!is.na(newx)) 
miss50_gamma <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx_CCA,
                         parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(miss50_gamma)
## CCA Log-Normal
miss50_lognorm <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = miss50_famx_CCA,
                         parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(miss50_lognorm)





###############################################################################################
########### Multiple Imputation without Considering Family Structure ##########################
###############################################################################################

## Step 1 - Empirical estimates
imp_model <- lm(newx ~ gender + proband + mgene + log(time) + status, data = miss50_famx) 
summary(imp_model)
betas <- coef(imp_model)
sigmahat <- sigma(imp_model)
V <- vcov(imp_model)
SSE <- sum((imp_model$residuals)^2)

fam2_imp <- list()
for (i in 1:20) {
  ## Step 2 - g ~ chi^2 nobs-p
  g <- rchisq(n = 1, df = 1651)
  
  ## Step 3 - sigma star
  sigmastar <- sigmahat/(SSE/g)
  
  ## Step 4 - u1
  u1 <- c()
  u1 <- rnorm(n = length(betas), mean = 0, sd = 1)
  
  ## Step 5
  betastar <- betas + (sigmastar/sigmahat) * u1 %*% chol(V)
  
  ## Step 6
  #num_miss <- sum(is.na(brca1_prs$PRS))
  #u2i <- rnorm(n = 1, mean = 0, sd = 1)
  
  ## Step 7
  fam2_imp[[i]] <- miss20_fam2 |>
    mutate(newx_I = ifelse(is.na(newx), betastar[,1] + betastar[,2]*gender + betastar[,3]*proband + betastar[,4]*mgene + betastar[,5]*log(time) + betastar[,6]*status + rnorm(n = 1, mean = 0, sd = 1)*sigmastar, newx))
}

gamma_results_nofamstruc <- list()
for (i in 1:20) {
  gamma_results_nofamstruc[[i]] <- penmodel(Surv(time, status) ~ gender + mgene + newx_I, cluster = "famID", gvar = "mgene", 
                                 design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = fam2_imp[[i]],
                                 parms = c(1/41.41327,1,0,0,0, 1)) 
}

results2 <- colMeans(do.call(rbind, gamma_results_nofamstruc))
nomiss_gamma$par
miss20_gamma_CCA$par
results1

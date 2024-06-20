######################################################
################## Simulation 1 ######################
######################################################
####### Missingness generated from mice package ######
######################################################

## Generate kinship data - Number of family = 200
install.packages("FamEvent_3.1.tar.gz", type = "source")
library(FamEvent)
library(mice)

set.seed(123)
famx <- simfam(N.fam = 200, design = "pop", variation = "frailty", 
               base.dist = "Weibull", frailty.dist = "gamma", interaction = FALSE,
               add.x = TRUE, x.dist = "normal", x.parms = c(0, 1),  depend = 2, 
               base.parms = c(0.016,3), vbeta = c(1, 3, 3)) 

## Test if Gamma frailty model gives the correct estimates ##
test_model <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(test_model)
## Test if LogNormal frailty gives the correct estimates ##
test_model <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                       design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = famx,
                       parms = c(1/41.41327,1,0,0,0, 1)) # True params
summary(test_model)

################# Suppose 50% missing - Missing at Random ####################
################# Suppose a reasonable weights ###############################
miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)
ampute_test <- ampute(data = famx, prop = 0.5, patterns = miss_pattern)
md.pattern(ampute_test$amp, rotate.names = TRUE) # Missing pattern
ampute_test$mech
ampute_test$weights
reasonable_weights <- ampute_test$weights
reasonable_weights[1,] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0)
ampute_test <- ampute(data = famx, prop = 0.5, patterns = miss_pattern, 
                      mech = "MAR", weights = reasonable_weights) # Missing proportions
miss50_famx <- ampute_test$amp
## CCA Gamma
miss50_famx_CCA <- miss50_famx |>
  filter(!is.na(newx)) |>
  group_by(famID) |>
  filter(any(proband == 1)) |>
  ungroup() # This deletes families that contain the proband being missing
miss50_gamma <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx_CCA,
                         parms = c(1/41.41327,1,0,0,0, 1)) 
summary(miss50_gamma)
## CCA Log-Normal
miss50_lognorm <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                           design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = miss50_famx_CCA,
                           parms = c(1/41.41327,1,0,0,0, 1)) 
summary(miss50_lognorm)

################# If you want to reproduce the error sent via email ###############
miss_pattern <- matrix(c(rep(1, 14), 0, rep(1, 4)), nrow = 1, byrow = TRUE)
ampute_test <- ampute(data = famx, prop = 0.5, patterns = miss_pattern)
md.pattern(ampute_test$amp, rotate.names = TRUE) # Missing pattern
ampute_test$mech
ampute_test$weights
reasonable_weights <- ampute_test$weights
reasonable_weights[1,] <- c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0)
ampute_test <- ampute(data = famx, prop = 0.5, patterns = miss_pattern, 
                      mech = "MAR", weights = reasonable_weights) # Missing proportions
miss50_famx <- ampute_test$amp

## CCA
miss50_famx_CCA <- miss50_famx |>
  filter(!is.na(newx)) 

## CCA Gamma
miss50_gamma <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                         design = "pop", base.dist = "Weibull", frailty.dist = "gamma", agemin = 20, data = miss50_famx_CCA,
                         parms = c(1/41.41327,1,0,0,0, 1)) 
summary(miss50_gamma)
## CCA Log-Normal
miss50_lognorm <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                           design = "pop", base.dist = "Weibull", frailty.dist = "lognormal", agemin = 20, data = miss50_famx_CCA,
                           parms = c(1/41.41327,1,0,0,0, 1)) 
summary(miss50_lognorm)


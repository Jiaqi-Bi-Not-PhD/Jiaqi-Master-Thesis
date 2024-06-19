############################################################
############# CCA Functions used for simulation ############
############################################################
CCA_test <- penmodel(Surv(time, status) ~ gender + mgene + newx, cluster = "famID", gvar = "mgene", 
                                   design = "pop", base.dist = "Weibull", frailty.dist = "gamma", 
                                   agemin = min(famx$currentage[famx$status == 1]), data = miss50_famx,
                                   parms = c(1/41.41327,1,0,0,0, 1))

summary(CCA_test)

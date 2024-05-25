## Terminal window: cd Dropbox/FamEvent_update_project/March\ 2024/
#R CMD build FamEvent
#R CMD check --as-cran FamEvent_3.1.tar.gz

# In R 
  
setwd("~/Dropbox/FamEvent_update_project/March 2024/")
#remove.packages("FamEvent")
install.packages("~/Dropbox/FamEvent_update_project/March 2024/FamEvent_3.1.tar.gz", repos = NULL, type="source")
library(FamEvent)
library(MASS)
library(truncnorm)
library(kinship2)
library(matrixcalc)


# kinship
fam1 <- simfam(N.fam = 10, design = "noasc", variation = "kinship",
               base.dist = "Weibull", frailty.dist = "lognormal", depend = 0.5,
               base.parms = c(0.0165,3), vbeta = c(1.5, 0.5))

summary(fam1)
fam1 <- simfam(N.fam = 10, design = "noasc", variation = "kinship",
               add.x = FALSE, x.dist = NULL, x.parms = NULL, depend = 0.5,
               base.dist = "Weibull", frailty.dist = "lognormal", 
               base.parms = c(0.0165,3), vbeta = c(1.5, 0.5))
summary(fam1)
## no ascertainment correction
fam0 <- simfam(N.fam = 50, design = "noasc", variation = "frailty",
                       base.dist = "Weibull", frailty.dist = "gamma", depend = 0.5, 
                       base.parms = c(0.0165,3), vbeta = c(1.5, 0.5))

### Additional covariate
## continuous covariate generated from normal distribution 
famx <- simfam(N.fam = 10, design = "noasc", variation = "none",
               base.dist = "Weibull", frailty.dist = NULL, interaction = FALSE,
               add.x = TRUE, x.dist = "normal", x.parms = c(0, 1),  depend = NULL, 
               base.parms = c(0.016,3), vbeta = c(1.5, 1, 1.5))

summary(famx)
## binary covariate generated from binomial distribution 
famx <- simfam(N.fam = 10, design = "pop", variation = "kinship",
               base.dist = "Weibull", frailty.dist = "lognormal", interaction = FALSE,
               add.x = TRUE, x.dist = "binomial", x.parms = c(1, 0.5),  depend = 0.5,
               base.parms = c(0.016,3), vbeta = c(1.5, 1, 1.5))


### Frailties from two sources of variation (kinship  + IBD)

IBDmatrix <- diag(1, dim(famx)[1])
inputdata <- famx[, !names(famx) %in%c("time", "status", "ageonset", "fsize", "naff")]

#inputdata should contain variables:
#  famID, indID, gender, motherID, fatherID, proband, currentage, proband and other variables being used in the model
fam2 <- simfam2(design = "pop", variation = c("kinship","IBD"), 
                depend = c(1, 1), base.dist="Weibull", base.parms =c(0.016, 3),
                var_names = c("gender", "mgene", "newx"), vbeta = c(1,1,1),
                agemin=20, inputdata=inputdata, IBD=IBDmatrix) 

summary(fam2)
famibd <- simfam2(design = "pop", variation = c("IBD"), 
                  depend = c(1), base.dist="Weibull", base.parms =c(0.016,3),
                  var_names = c("gender", "mgene"), vbeta = c(1, 1),
                  agemin=20, inputdata=inputdata, IBD=IBDmatrix) 

famkin <- simfam2(design = "pop", variation = "kinship", 
                  depend = 1, base.dist="Weibull", base.parms =c(0.016,3),
                  var_names = c("gender", "mgene", "newx"), vbeta = c(1, 1, 1),
                  agemin=20, inputdata=inputdata, IBD=IBDmatrix) 

summary(famkin)
source("FamEvent/R/simfam.R")
source("FamEvent/R/familyStructure.R")
source("FamEvent/R/familyDesign.R")
# kinship and IBD
source("FamEvent/R/simfam2.R")
source("FamEvent/R/familyStructure2.R")
source("FamEvent/R/familyDesign2.R")

source("FamEvent/R/surv.dist.R")
source("FamEvent/R/survp.dist.R")
source("FamEvent/R/inv.surv.R")
source("FamEvent/R/inv.survp.R")
source("FamEvent/R/inv2.surv.R")
source("FamEvent/R/fgene.R")
source("FamEvent/R/fgeneZ.R")
source("FamEvent/R/fgeneZX.R")
source("FamEvent/R/Pgene.R")
source("FamEvent/R/kids.g.R")
source("FamEvent/R/parents.g.R")
source("FamEvent/R/laplace.R")

# competing risk models
source("FamEvent/R/simfam_c.R")
source("FamEvent/R/familyStructure_c.R")
source("FamEvent/R/familyDesign_c.R")
source("FamEvent/R/fgene_c.R")
source("FamEvent/R/surv_dist_c.R")
source("FamEvent/R/inv_surv_c.R")
source("FamEvent/R/cumhaz.R")
source("FamEvent/R/hazards.R")


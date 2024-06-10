library(tidyverse)
library(ggplot2)
library(dplyr)

## Read dataset
prs_data <- read_csv("PRS_March2021.csv") # PRS data
brca1 <- read_csv("brca1.csv")
brca2 <- read_csv("brca2.csv")
brca_combined <- rbind(brca1, brca2)

## Check if there are observations in each cluster appear in two datasets
brca1_test <- brca1 %>%
  mutate(dataset = 1)
brca2_test <- brca2 %>%
  mutate(dataset = 2)
brca_combined_test <- rbind(brca1_test, brca2_test)

brca_combined_test <- brca_combined_test %>%
  group_by(famID) %>%
  summarize(repeated = any(length(unique(dataset)) != 1)) 
# Check if any obs from each family appeared in both brca1 & brca2
sum(brca_combined_test$repeated) # No repeat appeared
################################################################

## Plot missing pattern 
missings <- brca_combined %>%
  gather(key = "key", value = "val") %>%
  mutate(is.missing = is.na(val)) %>%
  group_by(key, is.missing) %>%
  summarise(num.missing = n()) %>%
  mutate(prop.missing = num.missing/nrow(brca_combined)) %>%
  filter(is.missing==T) %>%
  arrange(desc(num.missing)) 

missings %>%
  ggplot(aes(x = key, y = num.missing)) +
  geom_bar(stat = "identity")

missings %>%
  ggplot(aes(x = key, y = prop.missing)) + 
  geom_bar(stat = "identity")

missing.values <- brca_combined %>%
  gather(key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, isna) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)

levels <- (missing.values  %>% filter(isna == T) %>% arrange(desc(pct)))$key

missing.values %>%
  ggplot(aes(x = reorder(key, desc(pct)), y = pct, 
             fill=isna)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(limits = levels) +
  scale_fill_manual(name = "", 
                    values = c('blue', 'red'), 
                    labels = c("Present", "Missing"))  

brca_combined %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
  geom_raster(alpha=0.8) +
  scale_fill_manual(name = "",
                    values = c('blue', 'red'),
                    labels = c("Present", "Missing")) +
  scale_x_discrete(limits = levels) +
  labs(x = "Variable",
       y = "Row Number", title = "Missing values in rows") +
  coord_flip()

############### Data Merging ########################

### May 4 Meeting: 
### Likelihood of model in the thesis (without competing risk) (BC only)
### with family correlation 
### Write down the likelihood with only 2 covariate of mutation status (binary) and
### Consider PRS as another covariate (cts)
### Outcome Time-to-BC
### Write down the model and the likelihood


### Summarize the family member by different relationship - reference for
### those who have PRS non-missing

#total_data %>%
#  group_by(famID) %>%
#  any(!is.na(total_data$PRS))

#### Data merge using R script provided in OneDrive
#### Merged to BRCA data, some obs not in BRCA from PRS are omitted

brca.PRS<-merge(brca_combined, prs_data[,c("Person_ID","PRS")], by.x="indID", by.y="Person_ID", all.x=T)
sum(!is.na(brca.PRS$PRS)) #808 non-missing individuals
summary(brca.PRS$PRS)

sum(duplicated(brca.PRS$famID[!is.na(brca.PRS$PRS)])) #253 pairs non-missing

## Merge data by individual id
prs_data2 <- prs_data %>%
  mutate(indID = Person_ID) %>%
  dplyr::select(-"Person_ID")

brca1_prs <- left_join(brca1, prs_data2, by = "indID")
total_data <- left_join(brca_combined, prs_data2, by = "indID")
brca2_prs <- left_join(brca2, prs_data2, by = "indID")
nrow(brca2)-sum(is.na(brca2_prs$PRS))

length(unique(total_data$indID))==nrow(total_data) # Merged correctly
length(unique(brca1$indID)) == nrow(brca1) # Merged correctly

### Pedigree tree plot: (Need to fix the error)
### mgene imputation 0 1 *NA*
### PRS missing imputation / make inference without imputation
### Pedigree tree - who is present who is missing (PRS)
brca1_pedtree_test <- brca1_prs[1:6,] %>%
  mutate(gender = "female", 
         PRS_exist = ifelse(is.na(PRS), 0, 1),
         indID = as.character(indID), 
         famID = as.character(famID),
         fatherID = as.character(fatherID),
         motherID = as.character(motherID))

unique_fatherID <- unique(na.omit(brca1_pedtree_test$fatherID))
rows_to_add <- data.frame(indID = unique_fatherID, gender = "male",
                          famID = brca1_pedtree_test$famID, proband = 0,
                          affect = 0, PRS_exist = 0)
first_family <- bind_rows(brca1_pedtree_test, rows_to_add)[-1]
first_family <- first_family %>%
  distinct(indID, .keep_all = TRUE) %>%
  dplyr::select(c(indID, gender, famID, proband, affect, PRS_exist, fatherID, motherID)) 

  
# first_family <- data.frame(first_family)
library(kinship2)
ped_fam1 <- pedigree(id = first_family$indID,
                     dadid = first_family$fatherID,
                     momid = first_family$motherID,
                     sex = first_family$gender,
                     famid = first_family$famID, 
                     affected = cbind(first_family$affect, first_family$PRS_exist, first_family$proband))

# ped_toplot <- ped_fam1['1']

ped_toplot <- ped_fam1['10000001']
plot.pedigree(ped_toplot)

# Use pdf to save the pedigree plot
# pdf("pedigree.pdf", width = 11, height = 8)
#   plot.pedigree.shrink(ped_fam1)
#   dev.off()

# debug(plot.pedigree(ped_fam1))

### ** Examples of pedigree tree plot
library(kinship2)
data(sample.ped)
head(sample.ped)
head(first_family)

sample.ped1 <- sample.ped %>%
  filter(ped == 1)

pedAll <- pedigree(sample.ped1$id, sample.ped1$father, sample.ped1$mother,
                   sample.ped1$sex,  #affected=sample.ped$affected,
                   affected=cbind(sample.ped1$affected, sample.ped1$avail),
                   famid=sample.ped1$ped)

ped2 <- pedAll['2']

plot.pedigree(ped2)




length(unique(brca1_prs$famID))

brca1_test2 <- brca1_prs %>%
  mutate(sex = "F", 
         fatherID = ifelse(is.na(fatherID), 0, fatherID),
         motherID = ifelse(motherID %in% indID, motherID, 0)) 

unique_fathers <- unique(family_data$fatherID[family_data$fatherID != 0])
unique_mothers <- unique(family_data$motherID[family_data$motherID != 0])


plot_family_pedigree <- function(family_id, data) {
  # Filter the data for the specific family
  family_data <- data[data$famID == family_id, ]
  
  # Create a list of unique fathers and mothers
  unique_fathers <- unique(family_data$fatherID[family_data$fatherID != 0])
  unique_mothers <- unique(family_data$motherID[family_data$motherID != 0])
  
  # Create a data frame for dummy male individuals
  dummy_male_individuals <- data.frame(
    indID = unique_fathers,
    fatherID = 0,
    motherID = 0,
    sex = 1,
    proband = 0,
    PRS = NA,
    famID = family_id
  )
  
  # Combine the family data and dummy male individuals
  family_data_ext <- rbind(family_data, dummy_male_individuals)
  
  # Create a pedigree object using the kinship2 package
  pedigree <- with(family_data_ext, kinship2::pedigree(
    id = indID,
    dadid = fatherID,
    momid = motherID,
    sex = sex,
    affected = as.logical(proband),
    famid = famID
  ))
  
  # Define a function to add PRS labels
  add_prs_labels <- function(x, y, id, ...) {
    text(x, y + 0.5, labels = paste("PRS:", family_data_ext$PRS[match(id, family_data_ext$indID)]), ...)
  }
  
  # Plot the pedigree
  plot(pedigree, symbolsize = 1.5, label = pedigree$id,
       labelfont = 2, labelsize = 1.2, cexid = 0.8,
       missing = which(is.na(family_data_ext$PRS)),
       postplot = add_prs_labels)
}

plot_family_pedigree(10000001, brca1_test2)

ped1 <- with(brca1_test2, pedigree(id = indID, dadid = fatherID, momid = motherID, 
                                   affected = proband, famid = famID, sex = sex))



### For May 18 meeting:
###* Fit the model for the complete data - brca1 and PRS
### Likelihood function (FamEvent) - integrate the Z_j, treating Z_j ~ N
###* Write the likelihood with the integration, b/c PRS missing, frailty term needs to be integrated
### Extending Gamma loglik function from FamEvent
### Look at the code of penmodelem.R, learn the structure of FamEvent, check some examples
### See how the above function implements EM
### Model -> Likelihood -> How to maximize (EM)
library(tidyverse)
library(survival)
library(FamEvent)
library(EnvStats)
brca1_prs_fitting <- brca1_prs %>%
  filter(!is.na(mgene) & !is.na(PRS))
surv_object <- Surv(brca1_prs_fitting$timeBC, brca1_prs_fitting$BC)
model1 <- coxph(surv_object ~ mgene + PRS + frailty(famID, distribution = "gaussian"), 
                data = brca1_prs_fitting)

model2 <- coxph(surv_object ~ mgene + PRS + frailty(famID, distribution = "gamma"), 
                data = brca1_prs_fitting)
summary(model1)
summary(model2)

sum(is.na(brca1_prs$PRS))
sum(is.na(brca1_prs$mgene))
sum(is.na(brca1_prs$PRS) & is.na(brca1_prs$mgene))

# penmodelEM example without frailty
set.seed(4321)
fam <- simfam(N.fam = 100, design = "pop+", base.dist = "Weibull", base.parms = c(0.01, 3),
              vbeta = c(1, 2), agemin = 20, allelefreq = 0.02, mrate = 0.2)

fit <- penmodelEM(Surv(time, status) ~ gender + mgene, cluster = "famID", gvar = "mgene", 
                  parms = c(0.01, 3, 1, 2), data = fam, design="pop+", robust = TRUE, 
                  base.dist = "Weibull", method = "mendelian", mode = "dominant", q = 0.02)

summary(fit)
plot(fit)

# penmodelEM example with frailty
set.seed(4321)
fam <- simfam(N.fam = 100, design = "pop+", base.dist = "Weibull", base.parms = c(0.01, 3),
              vbeta = c(1, 2), agemin = 20, allelefreq = 0.02, mrate = 0.2, frailty.dist = "gamma", 
              variation = "frailty", depend = 2)

fam2 <- fam %>%
  filter(!is.na(mgene))
reg_model <- coxph(Surv(time, status) ~ gender + mgene + frailty(famID, distribution = "gamma"), 
                   data = fam2)
summary(reg_model)

# Estimate the baseline hazard weibull params
wei_model <- survreg(Surv(time, status) ~ 1, data = fam, dist = "weibull")
summary(wei_model)
lambda <- exp(-wei_model$coefficients)
rho <- 1/wei_model$scale

fit <- penmodelEM(Surv(time, status) ~ gender + mgene + frailty(famID, distribution = "gamma"), cluster = "famID", gvar = "mgene", 
                  parms = c(0.01, 3, 1, 2, 0), data = fam, design="pop+", robust = TRUE, 
                  base.dist = "Weibull", method = "mendelian", mode = "dominant", q = 0.02)

summary(fit)
plot(fit)
# penmodelEM example ends

# For brca1
brca1_prs_tofit <- brca1_prs %>%
  filter(!is.na(PRS) & !is.na(mgene) & !is.na(BC) & !is.na(timeBC))
library(parfm)
wei_model2 <- parfm(Surv(timeBC, BC) ~ mgene, cluster = "famID",
                    data = brca1_prs_tofit, dist = "weibull", frailty = "gamma")
lambda <- exp(-wei_model2$coefficients)
rho <- 1/wei_model2$scale

surv_object2 <- Surv(brca1_prs$timeBC, brca1_prs$BC)
fit_brca1 <- penmodelEM(surv_object2 ~ mgene + PRS ,
                 data = brca1_prs, cluster = "famID", gvar = "mgene", 
                 agemin = 18, parms = c(lambda, rho, 0.81, 0.13))

summary(fit_brca1)
plot(fit_brca1)

### Write down the likelihood in R
### Try to maximize the likelihood
### Read the PRS paper, to check whether brca1 & brca2 are included for prs computation
### Derive the likelihood from the FamEvent paper by relating to the topic
### (one more condition on the PRS - need to integrate)

brca1_prs_fitting %>%
  ggplot(aes(x = PRS)) +
  geom_density() +
  theme_minimal() +
  ylab("Density") +
  ggtitle("PRS Distribution in the Incomplete Data") +
  theme(plot.title = element_text(hjust = 0.5))
shapiro.test(brca1_prs_fitting$PRS)

### June 15 - Check if BRCA1 and BRCA2 include SNPs of the PRS paper
prs_vars <- colnames(prs_data)
prs_vars[-c(1:4)]
tables7 <- readxl::read_excel("prs paper Mavaddat et al table s7.xlsx")
tables8 <- readxl::read_excel("prs paper Mavaddat et al table s8.xlsx")
tables7_vars <- tables7$SNPa 
tables7_vars <- paste0("_", tables7_vars)
tables8_vars <- tables8$SNPa
tables8_vars <- paste0("_", tables8_vars)
common_SNPs_s7 <- intersect(tables7_vars, prs_vars) # All Table S7 elements are in PRS calculation
common_SNPs_s8 <- intersect(tables8_vars, prs_vars) # 289/3820 SNPs of S8 in PRS calculation
common_SNPs_s7s8 <- intersect(common_SNPs_s7, common_SNPs_s8)

### June 29
library(naniar)
brca1_prs <- brca1_prs %>%
  mutate(Missing_PRS = ifelse(is.na(PRS), 1, 0),
         Missing_mgene = ifelse(is.na(mgene), 1, 0))
brca1_prs %>%
  ggplot(aes(x = mgene, y = PRS)) +
  geom_miss_point()
sum(brca1_prs$Missing_PRS==1 & brca1_prs$Missing_mgene==1)
sum(brca1_prs$Missing_mgene == 0 & brca1_prs$mgene == 1 & brca1_prs$Missing_PRS == 0)

library(norm)
brca1_prs_mgene <- brca1_prs %>%
  select(famID, indID, PRS, mgene) %>%
  mutate(bothmissing = ifelse(is.na(PRS) & is.na(mgene), 1, 0),
         PRSmissing = ifelse(is.na(PRS) & !is.na(mgene), 1, 0),
         mgenemissing = ifelse(is.na(mgene) & !is.na(PRS), 1, 0))
pre <- prelim.norm(as.matrix(brca1_prs_mgene))
thetahat <- em.norm(pre)
getparam.norm(pre, thetahat)
imp <- imp.norm(pre, thetahat, brca1_prs_mgene)

## EM algorithm with other packages
data(mdata)
s <- prelim.norm(mdata)
thetahat <- em.norm(s)
rngseed(1234567)
ximp <- imp.norm(s, thetahat, mdata)

################# Try to write the em ####################
em_miss <- function(df, max_iter = 100, tol = 1e-3) {
  df_impute <- df
  families <- unique(df$famID)
  
  for (family in families) {
    
    df_family <- df[df$famID == family,]
  
  
    ## Initial params
    mu <- colMeans(df_family[, c("PRS", "mgene")], na.rm = TRUE)
    Sigma <- var(df_family[, c("PRS", "mgene")], na.rm = TRUE)
  
    for (i in 1:max_iter) {
      mu_old <- mu
      Sigma_old <- Sigma
    
     ## E-step
      df_family_impute <- df_family
      for (j in 1:nrow(df_family)) {
        ## Scenario 1: PRS missing, mgene observed
        if (is.na(df_family$PRS[j]) && !is.na(df_family$mgene[j])) {
          Sigma_22 <- Sigma[2, 2]
          Sigma_12 <- Sigma[1, 2]
          Sigma_21 <- Sigma[2, 1]
          Sigma_11 <- Sigma[1, 1]
          df_family_impute$PRS[j] <- mu[1] + Sigma_12 * (df_family$mgene[j] - mu[2])/Sigma_22
      } 
        ## Scenario 2: mgene missing, PRS observed
        else if (!is.na(df_family$PRS[j]) && is.na(df_family$mgene[j])) {
          Sigma_22 <- Sigma[2, 2]
          Sigma_12 <- Sigma[1, 2]
          Sigma_21 <- Sigma[2, 1]
          Sigma_11 <- Sigma[1, 1]
          df_family_impute$mgene[j] <- mu[2] + Sigma_21 * (df_family$PRS[j] - mu[1])/Sigma_11
      }
        ## Scenario 3: Both missing, initialize the mean of both covariates
        else if (is.na(df_family$PRS[j]) && is.na(df_family$mgene[j])) {
          df_family_impute$PRS[j] <- mu[1]
          df_family_impute$mgene[j] <- mu[2]
      }
    }
    ## M-step, update the params
    mu <- colMeans(df_family_impute[, c("PRS", "mgene")])
    Sigma <- var(df_family_impute[, c("PRS", "mgene")])
    
    ## mu or Sigma missing
    if(any(is.na(mu), is.na(Sigma))) {
      break
    }
    
    ## Convergence
    if (sum(abs(mu-mu_old), na.rm = TRUE) < tol && sum(abs(Sigma - Sigma_old), na.rm = TRUE) < tol) {
      break
    }
    }
    df_impute[df_impute$famID == family, c("PRS", "mgene")] <- df_family_impute[, c("PRS", "mgene")]
  }
  return(df_impute)
}
df <- brca1_prs_mgene
df_imputed <- em_miss(df)

########### Debug Zone #############
colMeans(brca1_prs_mgene, na.rm = TRUE)


### July. 23 Trying to find f(m|x) and f(x|m) ###
library(dplyr)
mu <- brca1_prs %>%
  mutate(mgene = factor(mgene)) %>%
  group_by(mgene) %>%
  filter(!is.na(PRS) & !is.na(mgene)) %>%
  dplyr::summarize(mean = mean(PRS))
brca1_prs %>%
  filter(!is.na(PRS) & !is.na(mgene)) %>%
  ggplot(aes(x = PRS, color = factor(mgene))) +
  geom_density() +
  geom_vline(data = mu, aes(xintercept = mean, color = mgene))

### MICE for coxph with frailty
library(mice)
brca1_prs_mice$na <- nelsonaalen(data = brca1_prs_mice, timevar = "timeBC",
                                 statusvar = "BC")
micesurv0 <- mice(brca1_prs_mice, maxit = 0, method = "rf")
micesurvmethod <- micesurv0$method
micesurvpred <- micesurv0$predictorMatrix
micesurvmethod[c("PRS")] <- "rf"
micesurvmethod[c("mgene")] <- "rf"
micesurvpred[, "indID"] <- 0

micesurv <- mice(brca1_prs_mice, method = micesurvmethod, predictorMatrix = micesurvpred,
                 m = 5, seed = 123)
results_mice <- with(micesurv, coxph(Surv(timeBC, BC) ~ mgene + PRS + frailty(famID, distribution = "gamma")))
summary(pool(results_mice), exponentiate = TRUE) # Here the MI does not incorporate random effect
# 

### parfm
library(parfm)

parfm_model <- parfm(Surv(timeBC, BC) ~ mgeneI, inip = c(2.5, 0.01, 2), iniFpar = 1, cluster = "famID", data = brca1, dist = "weibull", frailty = "gamma")
sum(is.na(brca1$BC))
################## Not FamEvent Log-Normal ##############
initial_params <- c(1, 1, 0, 0)
lower_bounds <- c(0.00001, 0.00001, -Inf, -Inf)
upper_bounts <- c(Inf, Inf, Inf, Inf)
result <- optim(par = initial_params, fn = log_likelihood, data = brca1_prs_cca, 
                event = "BC", 
                time = "timeBC", 
                cluster = "famID", 
                covariates = c("PRS", "mgeneI"), 
                N_p = 15, sigma = 1)

########## Gamma Frailty Llhd test ##########
brca1_prs_cca <- brca1_prs %>%
  filter(!is.na(mgene) & !is.na(PRS))
brca1_prs_cca <- brca1_prs %>%
  filter(!is.na(mgene) & !is.na(PRS)) %>%
  group_by(famID) %>%
  filter(any(proband == 1)) %>%
  ungroup()
source("~/Desktop/MSc Biostatistics Western/Thesis related articles/PRS data/GammaFrailty_llhd_Jiaqi.R")
source("~/Desktop/MSc Biostatistics Western/Thesis related articles/FamEvent/R/cumhaz.R")
source("~/Desktop/MSc Biostatistics Western/Thesis related articles/FamEvent/R/hazards.R")
initial_params <- c(log(2), log(2), 0, 0)
X <- as.matrix(data.frame(brca1_prs_cca$PRS, brca1_prs_cca$mgeneI), 
            nrow=nrow(brca1_prs_cca), 
               ncol = 2)
Y <- as.matrix(data.frame(brca1_prs_cca$timeBC, brca1_prs_cca$BC), 
            nrow = nrow(brca1_prs_cca), 
               ncol = 2)
lower_bounds <- c(-Inf, -Inf, -Inf, -Inf, 0)
upper_bounds <- c(Inf, Inf, Inf, Inf, Inf)
result <- optim(par = initial_params, fn = loglik_frailty_single_gamma,
                data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
                design = "pop", frailty.dist = "gamma", base.dist = "Weibull",
                agemin = 18, k = 1, 
                control = list(maxit = 2000))

coxph(Surv(timeBC, BC) ~ PRS + mgeneI + frailty(famID, distribution = "gamma"), 
      data = brca1_prs_cca)


library(frailtyHL)
library(survival)
check <- frailtyHL(Surv(timeBC, BC) ~ PRS + mgeneI + (1|famID), data = brca1_prs_cca, varfixed = TRUE, RandDist = "Gamma")
summary(check)

# Debug
problematic_params <- list()
result2 <- optim(par = initial_params, fn = gamma_llhd_frailty,
                data = brca1_prs_cca, event = "BC", time = "timeBC", 
                cluster = "famID", covariates = c("mgeneI", "PRS"), v = 100, method = "L-BFGS-B", 
                lower = lower_bounds, upper = upper_bounds)

global_results <- c()
global_results <- c(global_results, results2)

print(problematic_params)

fit2 <- parfm(formula = Surv(timeBC, BC) ~ mgeneI + PRS, data = brca1_prs_cca,
              cluster = "famID", dist = "weibull", frailty = "none")


fit2 <- coxph(Surv(timeBC, BC) ~ mgeneI + PRS + frailty(famID, distribution = "gamma"), data = brca1_prs_cca)

######### Log-Normal Frailty ############
source("~/Desktop/MSc Biostatistics Western/Thesis related articles/PRS data/LogNormalFrailty_llhd_Jiaqi.R")
initial_params <- c(log(2), log(2), 0,0)
X <- as.matrix(data.frame(brca1_prs_cca$PRS, brca1_prs_cca$mgeneI), 
               nrow=nrow(brca1_prs_cca), 
               ncol = 2)
Y <- as.matrix(data.frame(brca1_prs_cca$timeBC, brca1_prs_cca$BC), 
               nrow = nrow(brca1_prs_cca), 
               ncol = 2)
lower_bounds <- c(0.0001, 0.0001, -Inf, -Inf)
upper_bounds <- c(Inf, Inf, Inf, Inf)
result <- optim(par = initial_params, fn = lognormal_single,
                data = brca1_prs_cca, X = X, Y = Y, nbase = 2,
                design = "pop", frailty.dist = "lognormal", base.dist = "Weibull",
                agemin = 18, control = list(maxit = 2000))

coxph(Surv(timeBC, BC) ~ PRS + mgeneI + frailty(famID, distribution = "gaussian"), 
      data = brca1_prs_cca)

loglik_frailty_single_lognormal(theta = initial_params,
                                data = brca1_prs_cca, X = X, Y = Y, 
                                nbase = 2,
                                design = "pop", 
                                frailty.dist = "lognormal", base.dist = "Weibull",
                                agemin = min(brca1_prs_cca$timeBC[brca1_prs_cca$BC == 1]))

check_model <- coxph(Surv(timeBC, BC) ~ mgeneI + PRS + frailty(famID, distribution = "gamma"), data = brca1_prs_cca)
library(survival)
# Use log-transformation to get rid of baseline parameter lower bounds

############## MCEM ################

#############################################
#------Asthma dataset------
data(asthma)
head(asthma)
# type 'help(asthma)' for a description of the data set

asthma$time <- asthma$End - asthma$Begin
parfm(Surv(time, Status) ~ Drug, cluster = "Patid", data = asthma,
      dist = "weibull", frailty = "gamma")
### CCA
brca1_prs_cca <- brca1_prs %>%
  filter(!is.na(mgene) & !is.na(PRS))
results_cca <- coxph(Surv(timeBC, BC) ~ mgene + PRS + frailty(famID, distribution = "gamma"), 
                     data = brca1_prs_cca)
summary(results_cca)

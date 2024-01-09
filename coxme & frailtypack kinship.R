###### coxme & kinship matrix ######
library(coxme)
library(kinship2)
## Split by families, and kinship matrix for each family, return a list
families <- split(brca1_prs, brca1_prs$famID)
kinship_mat <- lapply(families, function(data) {
  with(data, kinship(id = indID, dadid = fatherID, momid = motherID, sex = rep(2, nrow(data))))
})
## a large kinship matrix for all
kinship_mat <- with(brca1_prs, kinship(id = indID, dadid = fatherID, momid = motherID,
                                       sex = rep(2,nrow(brca1_prs))))

## coxme
model_coxme <- coxme(Surv(timeBC, BC) ~ PRS + mgeneI + (1|indID),
                     varlist = list(kinship_mat),
                     data = brca1_prs)
summary(model_coxme)
###### frailtypack package ######
library(frailtypack)

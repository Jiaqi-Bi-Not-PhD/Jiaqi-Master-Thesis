## Load packages
library(tidyverse)
library(ggplot2)
library(dplyr)

## Read dataset
prs_data <- read_csv("PRS_March2021.csv") # PRS data
brca1 <- read_csv("brca1.csv")
brca2 <- read_csv("brca2.csv")

## Merge data by individual id
prs_data2 <- prs_data %>%
  mutate(indID = Person_ID) %>%
  dplyr::select(-"Person_ID")
brca1_prs <- left_join(brca1, prs_data2, by = "indID")

## For MCEM, if one family contains 0 observed PRS, delete this family
brca1_prs_MCEM <- brca1_prs |>
  mutate(PRS_obs = ifelse(!is.na(PRS), 1, 0)) |>
  group_by(famID) |>
  summarise(n_obs = sum(PRS_obs)) 

brca1_prs_MCEM <- merge(brca1_prs_MCEM, brca1_prs, by = "famID")
brca1_prs_MCEM <- brca1_prs_MCEM |>
  filter(n_obs != 0)

## If CCA is preferred (471 obs)
brca1_prs_cca <- brca1_prs %>%
  filter(!is.na(PRS)) %>%
  group_by(famID) %>%
  filter(any(proband == 1)) %>%
  ungroup()
## or (503 obs)
brca1_prs_cca <- brca1_prs %>%
  filter(!is.na(PRS)) 
##################################################################
####### Pedigree tree - who is present who is missing (PRS) ######
##################################################################
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

library(kinship2)
ped_fam1 <- pedigree(id = first_family$indID,
                     dadid = first_family$fatherID,
                     momid = first_family$motherID,
                     sex = first_family$gender,
                     famid = first_family$famID, 
                     affected = cbind(first_family$affect, first_family$PRS_exist, first_family$proband))

## First family pedigree tree to visualize the missing PRS information
ped_toplot <- ped_fam1['10000001']
plot.pedigree(ped_toplot)
##################################################################
####### ---------------End of Pedigree tree--------------- #######
##################################################################

## Missing mechanism specification
brca1_prs <- brca1_prs |>
  mutate(miss_index = ifelse(is.na(PRS) == TRUE, 1, 0))
model_miss <- glm(miss_index ~ mgeneI + log(timeBC) + currentage + proband, data = brca1_prs, family = binomial)
summary(model_miss)


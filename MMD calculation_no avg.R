## R-4.1.1 (RStudio-2021.09.0-351)
library(plyr)                          # version 1.8.6 
library(dplyr)                         # version 1.0.7

rm(list = ls())

source("Functions for indices calculation.R")

#### Load data ####
clean.HD <- read.csv("final working HD set.csv")
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", "sodium", "rbc")

# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) clean.HD$pltlt <- sqrt(clean.HD$pltlt)
clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(clean.HD[ , which(colnames(clean.HD) %in% c("wbc", "rdw", "Glucose"))])

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)

# Scale variables
dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, scale)

# Caclculate covariance matrix based on all individuals
covars = cov(dat[ , selected.biomarkersHD])
dat.mmd <- ddply(dat, .(id_no), mmd_calc, covars, selected.biomarkersHD, ma = T)

# Apply log transformation to MMD
dat.mmd$MMD <- log(dat.mmd$MMD)

write.csv(dat.mmd[ complete.cases(dat.mmd$MMD), 
                   c("id_no", "date_visit", "age_visit", "sex", "diabetes", "date_death", "date_last", "fu_length", "MMD")], 
          "MMD no avg.csv", 
          row.names = F)


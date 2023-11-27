## R-4.1.1 (RStudio-2021.09.0-351)
library(plyr)                          # version 1.8.6
library(dplyr)                         # version 1.0.7

rm(list = ls())

#######################################################################################################################
#########         With variables transformed as needed and by 6-month time windows, for the main text         #########
#######################################################################################################################
source("Functions for indices calculation.R")

#### Load data ####
clean.HD <- read.csv("final working HD set.csv")
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", "sodium", "rbc")

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)

months <- 6

##############################################################################################################
###  For CVPC1, MMD, and NMV, variables not normally distributed are first log or square root-transformed  ###
##############################################################################################################
dat.trsf <- dat
# Transform variables not normally distributed
if ("pltlt" %in% selected.biomarkersHD) dat.trsf$pltlt <- sqrt(dat.trsf$pltlt)
dat.trsf[ , which(colnames(dat.trsf) %in% c("wbc", "rdw", "Glucose"))] <- 
  log(dat.trsf[ , which(colnames(dat.trsf) %in% c("wbc", "rdw", "Glucose"))])


#### 1. CVPC1
dat.cvs <- ddply(dat.trsf, .(id_no), cv.4mth, selected.biomarkersHD = selected.biomarkersHD, 
                 nb.month = months, min.vis.nb = 4, SDs = apply(dat[ , selected.biomarkersHD] , 2, sd), 
                 covars = c("id_no", "sex", "date_death", "date_birth", "date_last", 
                            "diabetes", "fu_length"))

# Apply log transformation on CVs
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, log)

# Correct CVs for the number of visits included
dat.cvs[ , selected.biomarkersHD] <- apply(dat.cvs[ , selected.biomarkersHD], 2, function(x) {
  residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b, data = dat.cvs, start = list(a = 1, b = 1)))
})

# Add "status" variable (for cox regressions)
dat.cox <- ddply(dat.cvs, .(id_no), status_fct)

# Calculate PCA on CVs
pca_cv <- prcomp(dat.cox[ , selected.biomarkersHD], center = T, scale. = T)
dat1 <- cbind(dat.cox, pca_cv$x)
dat1 <- pc_sign(dat1, "PC1", pca_cv)
colnames(dat1)[which(colnames(dat1)=="PC1")] = "CVPC1"

##############################################################################################################
###          For all indices except CVPC1, variables are z-score transformed before calculation            ###
##############################################################################################################
# Scale variables
dat.trsf[ , selected.biomarkersHD] <- apply(dat.trsf[ , selected.biomarkersHD], 2, scale)
dat[ , selected.biomarkersHD] <- apply(dat[ , selected.biomarkersHD], 2, scale)

### 2. MMD
# Caclculate covariance matrix based on all individuals
covars = cov(dat.trsf[ , selected.biomarkersHD])
dat.mmd <- ddply(dat.trsf, .(id_no), mmd_calc, covars, selected.biomarkersHD, ma = T)

# Calculate average by 6 months
dat.mmd2 <- dat.mmd[complete.cases(dat.mmd$MMD), ]
dat.mmd3 <- ddply(dat.mmd2, .(id_no), index.tw, vars = selected.biomarkersHD, 
                  time.window = months, min.vis.nb = 4, index.nm = "MMD", mean = T)

dat2 <- merge(dat1, dat.mmd3, by = c("id_no", "years"), all = T)

### 3. Other indices from Weinans
# Degenerative fingerprinting
pca <- prcomp(dat[ , selected.biomarkersHD])
dat$PC1 <- pca$x[ , 1]
dat.df <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                time.window = months, min.vis.nb = 4, index.nm = "Df", index.fct = df.fct)
dat3 <- merge(dat2, dat.df, by = c("id_no", "years"), all = T)

# MAF autocorrelation
dat.maf.ac <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                    time.window = months, min.vis.nb = 13, index.nm = "MAF_ac", index.fct = maf.autocorr)
dat4 <- merge(dat3, dat.maf.ac, by = c("id_no", "years"), all = T)

# MAF eigenvalue
dat.maf.ev <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                    time.window = months, min.vis.nb = 4, index.nm = "MAF_ev", index.fct = maf.eigenval)
dat5 <- merge(dat4, dat.maf.ev, by = c("id_no", "years"), all = T)

# Mutual information
dat.mi <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                time.window = months, min.vis.nb = 4, index.nm = "MI", index.fct = mut.inf)
dat6 <- merge(dat5, dat.mi, by = c("id_no", "years"), all = T)

# Average autocorrelation
dat.av.ac <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                   time.window = months, min.vis.nb = 4, index.nm = "Av_Ac", index.fct = avg.autocorr)
dat7 <- merge(dat6, dat.av.ac, by = c("id_no", "years"), all = T)

# Node maximum autocorrelation
dat.nma <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                 time.window = months, min.vis.nb = 4, index.nm = "NMA", index.fct = node.max.autocorr)
dat8 <- merge(dat7, dat.nma, by = c("id_no", "years"), all = T)

# MAF variance
dat.maf.var <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                     time.window = months, min.vis.nb = 13, index.nm = "MAF_var", index.fct = maf.var)
dat9 <- merge(dat8, dat.maf.var, by = c("id_no", "years"), all = T)

# Node maximum variance
dat.nmv <- ddply(dat.trsf, .(id_no), index.tw, vars = selected.biomarkersHD,
                 time.window = months, min.vis.nb = 4, index.nm = "NMV", index.fct = node.max.var)
dat10 <- merge(dat9, dat.nmv, by = c("id_no", "years"), all = T)

# Average variance
dat.av.var <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                    time.window = months, min.vis.nb = 4, index.nm = "Av_Var", index.fct = avg.var)
dat11 <- merge(dat10, dat.av.var, by = c("id_no", "years"), all = T)

# PCA variance
dat.pc.var <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                    time.window = months, min.vis.nb = 4, index.nm = "PC_var", index.fct = pca.var)
dat12 <- merge(dat11, dat.pc.var, by = c("id_no", "years"), all = T)

# Maximum value of covariance matrix
dat.mvcm <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                  time.window = months, min.vis.nb = 4, index.nm = "Max_cov", index.fct = max.cov)
dat13 <- merge(dat12, dat.mvcm, by = c("id_no", "years"), all = T)

# Explained variance
dat.ev <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                time.window = months, min.vis.nb = 4, index.nm = "Ex_var", index.fct = expl.var)
dat14 <- merge(dat13, dat.ev, by = c("id_no", "years"), all = T)

# Average absolute cross-correlation
dat.aacc <- ddply(dat, .(id_no), index.tw, vars = selected.biomarkersHD,
                  time.window = months, min.vis.nb = 4, index.nm = "Av_ab_cc", index.fct = avg.abs.cc)
dat15 <- merge(dat14, dat.aacc, by = c("id_no", "years"), all = T)


indices <- c("CVPC1", "MMD", "Df", "MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA", "MAF_var", "NMV", "Av_Var",
             "PC_var", "Max_cov", "Ex_var", "Av_ab_cc")


####### Apply transformations to indices #######
dat.all <- dat15
sqrt.ind <- c("MMD", "Max_cov", "NMV", "PC_var")
new.trf.ind <- c("MAF_ac", "MAF_var")

new.trf <- function(x) {
  max.val <- max(x, na.rm = T)
  log(max.val + 0.01 - x) * -1
}

dat.all$Av_Var <- log(dat.all$Av_Var)
dat.all[ , sqrt.ind] <- apply(dat.all[ , sqrt.ind], 2, sqrt)
dat.all[ , new.trf.ind] <- apply(dat.all[ , new.trf.ind], 2, new.trf)

## Correct indices for the number of observations included in their calculation
blue <- c("Av_ab_cc")
red <- c("MAF_ev", "MAF_var")
green <- c("MMD", "Av_Ac", "Df", "MI", "Max_cov", "MAF_ac")
lin <- c("Av_Var", "Ex_var", "NMA", "NMV", "PC_var")

dat.all[!is.na(dat.all$Av_ab_cc), "Av_ab_cc"] <-
  residuals(nls(Av_ab_cc ~ (a/nb_vis) + b, data = dat.all, start = list(a = 1, b = 100)))

red.fct <- function(x) {
  residuals(nls(x ~ I(a / sqrt(nb_vis)) + b, data = dat.all, start = list(a = 1, b = 100)))
}
for (i in 1:length(red)) {
  dat.all[!is.na(dat.all[ , red[i]]), red[i]] <- red.fct(dat.all[ , red[i]])
}

green.fct <- function(x) {
  residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b * sqrt(nb_vis), data = dat.all, start = list(a = 1, b = 1)))
}
for (i in 1:length(green)) {
  dat.all[!is.na(dat.all[ , green[i]]), green[i]] <- green.fct(dat.all[ , green[i]])
}

dat.all[ , lin] <- apply(dat.all[ , lin], 2, function(x) {
  x[!is.na(x)] <- residuals(lm(x ~ nb_vis, data = dat.all))
  return(x)
})


write.csv(dat.all, "Indices.csv", row.names = F)

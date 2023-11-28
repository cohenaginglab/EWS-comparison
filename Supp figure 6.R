## R-4.3.1 (RStudio-2023.06.0)
library(plyr)                          # version 1.8.8 
library(dplyr)                         # version 1.1.2
library(survival)                      # version 3.5-5
library(splines)                       # version 4.3.1
library(viridis)                       # version 0.6.3


rm(list = ls())

### Load indices 
# Calculated by 6-month time windows
dat_6mth <- read.csv("Indices.csv")
# Calculated using various time windows
dats <- readRDS("HD indices_diff time windows.rds")

# Include the 6-month version with the others, between the 4-month and the 1-year versions
dats2 <- dats[1:3]
dats2[[4]] <- dat_6mth
dats2[[5]] <- dats[[4]]
for (j in 1:length(dats2)) dats2[[j]]$tw <- c("2 months", "3 months", "4 months", "6 months", "12 months")[j]

indices <- c("CVPC1", "MMD", "Df", "MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA", "MAF_var", "NMV", "Av_Var", 
             "PC_var", "Max_cov", "Ex_var", "Av_ab_cc")
# Order by categories of indices (variance-based, etc.)
ind.cat <- list("Variance" = c("CVPC1", "NMV", "Av_Var", "PC_var", "Max_cov", "Ex_var"),
                "Var + AC" = c("MMD", "Df", "MAF_var"),
                "Auto-correlation" = c("MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA"),
                "Cross-correlation" =  "Av_ab_cc")

## Run Cox models for mortality prediction
cnms <- c("Index", "model", "tw", "n", "event", "HR", "LCI", "UCI", "p", "HR95", "HR95_LCI", "HR95_UCI")
# Exclude MAF_var and MAF_ac 
indices0 <- setdiff(unlist(ind.cat), c("MAF_var", "MAF_ac"))

res.tw <- list()
for (k in 1:length(dats2)) {
  
  dat <- dats2[[k]]
  
  res.all <- list()
  for (i in 1:length(indices0)) {
    res <- as.data.frame(matrix(NA, 2, length(cnms)))
    colnames(res) <- cnms
    res$Index <- indices0[i]
    res$model <- 1:2
    res$tw <- c("2 months", "3 months", "4 months", "6 months", "12 months")[k]
    
    ## model 1: the index only
    dat_var <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                                "nb_vis", indices0[i])])
    colnames(dat_var)[ncol(dat_var)] <- "predictor"
    res$n[1] <- nrow(dat_var)
    res$event[1] <- sum(dat_var$status==1)
    mod <- coxph(Surv(years, status) ~ predictor + bs(age_visit, df = 5) + sex +  
                   diabetes + fu_length, data = dat_var, weights = sqrt(nb_vis),
                 cluster = id_no)
    res[1, c("HR", "LCI", "UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)]
    res[1, "p"] <- summary(mod)$coefficients[1, 6]
    res[1, c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)] ^ 
      (quantile(na.omit(dat_var[ , "predictor"]), 0.975) - quantile(na.omit(dat_var[ , "predictor"]), 0.025))
    
    ## model 2: all 13 indices
    dat_var <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                                "nb_vis", indices0)])
    res$n[2] <- nrow(dat_var)
    res$event[2] <- sum(dat_var$status==1)
    fmla <- formula(paste("Surv(years, status) ~ ", 
                          paste(indices0, collapse = " + "),
                          "+ bs(age_visit, df = 5) + sex + diabetes + fu_length", sep = " "))
    mod <- coxph(fmla, data = dat_var, weights = sqrt(nb_vis), cluster = id_no)
    mod.tab <- summary(mod)$conf.int
    res[2, c("HR", "LCI", "UCI")] <- mod.tab[which(rownames(mod.tab)==indices0[i]), c(1, 3, 4)]
    res[2, "p"] <- summary(mod)$coefficients[which(rownames(summary(mod)$coefficients)==indices0[i]), 6]
    res[2, c("HR95", "HR95_LCI", "HR95_UCI")] <- mod.tab[which(rownames(mod.tab)==indices0[i]), c(1, 3, 4)] ^ 
      (quantile(na.omit(dat_var[ , indices0[i]]), 0.975) - quantile(na.omit(dat_var[ , indices0[i]]), 0.025))
    res.all[[i]] <- res
  }
  res.tw[[k]] <- do.call(rbind, res.all)
}
res <- do.call(rbind, res.tw)

# Add unaveraged MMD
dat_MMD <- read.csv("MMD no avg.csv")
dat_MMD[ , c("date_death", "date_last", "date_visit")] <- 
  lapply(dat_MMD[ , c("date_death", "date_last", "date_visit")], as.Date)
## Add "status" variable
dat1 <- ddply(dat_MMD, .(id_no), function(dat.id) {
  dat.id$status <- ifelse(!is.na(dat.id$date_death) & dat.id$date_visit==max(dat.id$date_visit), 1, 0)
  dat.id$time <- as.numeric(unique(dat.id$date_last) - dat.id$date_visit) / 365.25
  return(dat.id)
})
res0 <- as.data.frame(matrix(NA, 1, length(cnms)))
colnames(res0) <- cnms
res0$Index <- "MMD_all"
res0$model <- 1
dat_var <- na.omit(dat1[ , c("id_no", "time", "status", "age_visit", "sex", "diabetes", "fu_length","MMD")])
res0$n <- nrow(dat_var)
res0$event <- sum(dat_var$status==1)
mod <- coxph(Surv(time, status) ~ MMD + bs(age_visit, df = 5) + sex + diabetes + fu_length, data = dat_var, cluster = id_no)
res0[ , c("HR", "LCI", "UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)]
res0[ , "p"] <- summary(mod)$coefficients[1, 6]
res0[ , c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)] ^ 
  (quantile(na.omit(dat_var[ , "MMD"]), 0.975) - quantile(na.omit(dat_var[ , "MMD"]), 0.025))
res <- rbind(res, res0)

indices2 <- setdiff(c(unlist(ind.cat)[1:7], "MMD_all", unlist(ind.cat)[8:15]), c("MAF_var", "MAF_ac"))

# Manual jitter
ypos <- c(-0.2, -0.1, 0, 0.1, 0.2)

# Set color palette
myPal <- rev(viridis(6)[-6])
pchs <- c(15, 17:20)

# Plot results
tiff("Fig_S6.tiff", width = 200, height = 110, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(1, 2), mar = c(3, 4.5, 1.5, 0.2), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 0.8)

## Forest plots (one per model)
for (j in 1:2) {
  res.mod <- res[res$model==j, ]
  
  # Set x-axis range 
  if (j==1) x_range <- c(0, 30) else x_range <- c(0, 90)
  
  nn <- length(levels(factor(res.mod$Index)))
  plot(1:10, 1:10, type = "n", xlim = x_range, ylim = c(0, nn - 0.5), 
       ylab = "", xlab = "HR95", yaxt = "n", frame.plot = F, font.lab = 2, cex.axis = 0.8)
  
  abline(v = 1, lty = 3)
  mtext(c("A", "B")[j], side = 1, line = -23, cex = 1.1, font = 2, at = ifelse(j==1, -8, -23))
  
  # Add variable names
  if (j==1) txt.x <- -7 else txt.x <- -20
  if (j==1) txt_labs <- indices2 else txt_labs <- setdiff(indices2, "MMD_all")
  text(x = rep(txt.x, nn), y = (nn - 1):0, labels = txt_labs, xpd = T, adj = 0)
  
  for (i in 1:length(txt_labs)) {
    res.ind <- res.mod[res.mod$Index==txt_labs[i], ]
    
    for (k in 1:nrow(res.ind)) {
      points(x = res.ind[k, "HR95"], y = nn - i + ypos[k], col = ifelse(j==1&i==8, "black", myPal[k]), 
             pch = pchs[k], cex = 1.2)
      if (i==3 & j==2 & k > 4) {
        arrows(x0 = res.ind[k, "HR95_LCI"], x1 = 85, y0 = nn - i + ypos[k], y1 = nn - i + ypos[k], 
               length = 0.1, col = myPal[k], lwd = ifelse(res.ind$p[k] < 0.05, 2, 1))
      } else {
        segments(x0 = res.ind[k, "HR95_LCI"], x1 = res.ind[k, "HR95_UCI"], y0 = nn - i + ypos[k], y1 = nn - i + ypos[k], 
                 col = ifelse(i==8, "black", myPal[k]), lwd = ifelse(res.ind$p[k] < 0.05, 2, 1))
      }
    }
  }
  
  # Legend
  if (j==2) legend(x = 50, y = 6, legend = c("2 months", "3 months", "4 months", "6 months", "12 months"),
                   col = myPal, pch = pchs, lwd = 2, bty = "n")
}

dev.off()
par(opar)


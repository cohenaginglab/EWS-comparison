## R-4.3.1 (RStudio-2023.06.0)
library(plyr)                          # version 1.8.8 
library(dplyr)                         # version 1.1.2
library(pROC)                          # version 1.18.2
library(survival)                      # version 3.5-5
library(splines)                       # version 4.3.1
library(pals)                          # version 1.7


rm(list = ls())

### Load indices 
dat <- read.csv("Indices.csv")

indices <- c("CVPC1", "MMD", "Df", "MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA", "MAF_var", "NMV", "Av_Var", 
             "PC_var", "Max_cov", "Ex_var", "Av_ab_cc")

# Order by categories of indices (variance-based, etc.)
ind.cat <- list("Variance" = c("CVPC1", "NMV", "Av_Var", "PC_var", "Max_cov", "Ex_var"),
                "Var + AC" = c("MMD", "Df", "MAF_var"),
                "Autocorrelation" = c("MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA"),
                "Cross-correlation" =  "Av_ab_cc")

##### 4. Mortality prediction
cnms <- c("Index", "model", "n", "event", "HR", "LCI", "UCI", "p", "HR95", "HR95_LCI", "HR95_UCI")

res.all <- list()
for (i in 1:length(indices)) {
  if (unlist(ind.cat)[i] %in% c("MAF_ac", "MAF_var")) nn <- 2 else nn <- 3
  res <- as.data.frame(matrix(NA, nn, length(cnms)))
  colnames(res) <- cnms
  res$Index <- unlist(ind.cat)[i]
  res$model <- 1:nn
  
  ## model 1: the index only
  dat_var <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                              "nb_vis", unlist(ind.cat)[i])])
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
  
  ## model 2: all 15 indices
  dat_var <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                              "nb_vis", indices)])
  res$n[2] <- nrow(dat_var)
  res$event[2] <- sum(dat_var$status==1)
  fmla <- formula(paste("Surv(years, status) ~ ", 
                        paste(indices, collapse = " + "),
                        "+ bs(age_visit, df = 5) + sex + diabetes + fu_length", sep = " "))
  mod <- coxph(fmla, data = dat_var, weights = sqrt(nb_vis), cluster = id_no)
  mod.tab <- summary(mod)$conf.int
  res[2, c("HR", "LCI", "UCI")] <- mod.tab[which(rownames(mod.tab)==unlist(ind.cat)[i]), c(1, 3, 4)]
  res[2, "p"] <- summary(mod)$coefficients[which(rownames(summary(mod)$coefficients)==unlist(ind.cat)[i]), 6]
  res[2, c("HR95", "HR95_LCI", "HR95_UCI")] <- mod.tab[which(rownames(mod.tab)==unlist(ind.cat)[i]), c(1, 3, 4)] ^ 
    (quantile(na.omit(dat_var[ , unlist(ind.cat)[i]]), 0.975) - quantile(na.omit(dat_var[ , unlist(ind.cat)[i]]), 0.025))
  
  ## model 3: all indices except MAF_ac and MAF_var (smaller sample size)
  if (nn > 2) {
    indices2 <- setdiff(indices, c("MAF_ac", "MAF_var"))
    dat_var <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                                "nb_vis", indices2)])
    res$n[3] <- nrow(dat_var)
    res$event[3] <- sum(dat_var$status==1)
    fmla <- formula(paste("Surv(years, status) ~ ", 
                          paste(indices2, collapse = " + "),
                          "+ bs(age_visit, df = 5) + sex + diabetes + fu_length", sep = " "))
    mod <- coxph(fmla, data = dat_var, weights = sqrt(nb_vis), cluster = id_no)
    mod.tab <- summary(mod)$conf.int
    res[3, c("HR", "LCI", "UCI")] <- mod.tab[which(rownames(mod.tab)==unlist(ind.cat)[i]), c(1, 3, 4)]
    res[3, "p"] <- summary(mod)$coefficients[which(rownames(summary(mod)$coefficients)==unlist(ind.cat)[i]), 6]
    res[3, c("HR95", "HR95_LCI", "HR95_UCI")] <- mod.tab[which(rownames(mod.tab)==unlist(ind.cat)[i]), c(1, 3, 4)] ^ 
      (quantile(na.omit(dat_var[ , unlist(ind.cat)[i]]), 0.975) - quantile(na.omit(dat_var[ , unlist(ind.cat)[i]]), 0.025))
  }
  
  res.all[[i]] <- res
}
res <- do.call(rbind, res.all)

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
dat_var <- na.omit(dat1[ , c("id_no", "time", "status", "age_visit", "sex", "diabetes", "fu_length", "MMD")])
res0$n <- nrow(dat_var)
res0$event <- sum(dat_var$status==1)
mod <- coxph(Surv(time, status) ~ MMD + bs(age_visit, df = 5) + sex + diabetes + fu_length, data = dat_var, cluster = id_no)
res0[ , c("HR", "LCI", "UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)]
res0[ , "p"] <- summary(mod)$coefficients[1, 6]
res0[ , c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod)$conf.int[1, c(1, 3, 4)] ^ 
  (quantile(na.omit(dat_var[ , "MMD"]), 0.975) - quantile(na.omit(dat_var[ , "MMD"]), 0.025))

res <- rbind(res, res0)

# Set color palette
pal1 <- brewer.blues(7)[-1]
pal2 <- brewer.purples(4)[-1]
pal3 <- brewer.reds(6)[-1]
pal4 <- "#00955e"

l.pos <- c(-7.5, -7.5, -12.7)

# Plot results
tiff("Fig_4.tiff", width = 210, height = 180, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(2, 2), mar = c(3, 4.5, 1.5, 0.2), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 0.8)

## Forest plots (one per model)
for (j in 1:3) {
  res.mod <- res[res$model==j, ]
  if (j==1)  res.mod <- res.mod[match(c(ind.cat$Variance, "MMD", "MMD_all", ind.cat$`Var + AC`[-1],
                                        unlist(ind.cat[3:4])), res.mod$Index), ]
  if (j==2) res.mod <- res.mod[match(unlist(ind.cat), res.mod$Index), ]
  if (j==3) res.mod <- res.mod[match(setdiff(unlist(ind.cat), c("MAF_ac", "MAF_var")), res.mod$Index), ]

  # Set x-axis range 
  if (j < 3) x_range <- c(0, 30) else x_range <- c(0, 50)
  
  plot(1:10, 1:10, type = "n", xlim = x_range, ylim = c(0, nrow(res.mod) - 0.5), 
       ylab = "", xlab = "HR95", yaxt = "n", frame.plot = F, font.lab = 2, cex.axis = 0.8)
  
  abline(v = 1, lty = 3)
  mtext(c("A", "B", "C")[j], side = 1, line = -18.5, cex = 1.1, font = 2, at = l.pos[j])
  
  # Set color palette
  if (j==1) myPal <- c(pal1, c("#CBC9E2", "#CBC9E2", "#9E9AC8", "#6A51A3"), pal3, pal4) 
  if (j==2) myPal <- c(pal1, pal2, pal3, pal4) 
  if (j==3) myPal <- c(pal1, pal2, pal3, pal4)[-c(9:10)]
  
  # Add variable names
  if (j < 3) txt.x <- -7 else txt.x <- -12
  text(x = rep(txt.x, nrow(res.mod)), y = (nrow(res.mod) - 1):0, labels = res.mod$Index, xpd = T, adj = 0, col = myPal)
  
  for (i in 1:nrow(res.mod)) {
    if (i==3 & j > 1) {
      arrows(x0 = res.mod[i, "HR95_LCI"], x1 = ifelse(j==3, 49, 29), 
             y0 = nrow(res.mod) - i, y1 = nrow(res.mod) - i, col = myPal[i], length = 0.08, lwd = 2)
    } else {
      segments(x0 = res.mod[i, "HR95_LCI"], x1 = res.mod[i, "HR95_UCI"], 
               y0 = nrow(res.mod) - i, y1 = nrow(res.mod) - i, col = myPal[i],
               lwd = ifelse(res.mod$p[i] < 0.05, 2, 1))
      segments(x0 = res.mod[i, "HR95_UCI"], x1 = res.mod[i, "HR95_UCI"], 
               y0 = nrow(res.mod) - i + 0.15, y1 = nrow(res.mod) - i - 0.15, 
               col = myPal[i], lwd = ifelse(res.mod$p[i] < 0.05, 2, 1))
    }
    points(x = res.mod[i, "HR95"], y = nrow(res.mod) - i, col = myPal[i], pch = 16, cex = 1.2)
    segments(x0 = res.mod[i, "HR95_LCI"], x1 = res.mod[i, "HR95_LCI"], 
             y0 = nrow(res.mod) - i + 0.15, y1 = nrow(res.mod) - i - 0.15, 
             col = myPal[i], lwd = ifelse(res.mod$p[i] < 0.05, 2, 1))
  }
}

## ROC curves with models including indices the most predictive indices
cbn <- list("CVPC1", "Av_Var")
cbn[[3]] <- c("CVPC1", "Av_Var")
cbn[[4]] <- c("CVPC1", "Av_Var", "PC_var")
cbn[[5]] <- "MMD_all"

rocs <- list()
for (j in 0:length(cbn)) {
  # First model: only demographic variables
  if (j==0) {
    mod <- coxph(Surv(years, status) ~ bs(age_visit, df = 5) + sex + diabetes + fu_length, data = dat, 
                 weights = sqrt(nb_vis), cluster = id_no)
    rocs[[j + 1]] <- roc(dat$status ~ predict(mod))
  } else {
    if (j==5) {
      mod <- coxph(Surv(time, status) ~ MMD + bs(age_visit, df = 5) + sex + diabetes + fu_length, data = dat1, 
                   cluster = id_no)
      rocs[[j + 1]] <- roc(dat1$status ~ predict(mod))
    } else {
      fmla <- formula(paste("Surv(years, status) ~ ", paste(cbn[[j]], collapse = " + "), 
                            "+ bs(age_visit, df = 5) + sex + diabetes + fu_length", sep = ""))
      mod <- coxph(fmla, data = dat, weights = sqrt(nb_vis), cluster = id_no)
      rocs[[j + 1]] <- roc(dat$status ~ predict(mod))
    }
 }
}

## Add model with all indices (13)
indices2 <- setdiff(indices, c("MAF_ac", "MAF_var"))
dat_var <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                            "nb_vis", indices2)])
fmla <- formula(paste("Surv(years, status) ~ ", paste(indices2, collapse = " + "), 
                      "+ bs(age_visit, df = 5) + sex + diabetes + fu_length", sep = ""))
mod <- coxph(fmla, data = dat_var, weights = sqrt(nb_vis), cluster = id_no)
rocs[[length(rocs) + 1]] <- roc(dat_var$status ~ predict(mod))

par(mar = c(2.5, 2.5, 0.5, 0.5), cex = 0.9)
cols <- c("black", "#df2b5a", "#193da5", "#dfb02b", "#008856", "#b02bdf", "#6dcce9")
plot(1:10, 1:10, type = "n", xlim = c(1, 0), ylim = c(0, 1), xlab = "Specificity",
     ylab = "Sensitivity", frame.plot = F, cex.axis = 0.8)
mtext("D", side = 1, line = -17, cex = 1.1, font = 2, at = 1.15)
segments(x0 = 1, y0 = 0, x1 = 0, y1 = 1, col = "grey")
for (i in 1:length(rocs)) plot(rocs[[i]], add = T, col = cols[i])
legend(x = 0.65, y = 0.37, 
       legend = paste(c("Control variables", unlist(lapply(cbn, function(x) paste(x, collapse = " + "))),
                        "13 indices"),
                      paste("(", unlist(lapply(rocs, function(x) round(x$auc, 3))), ")", sep = "")),
       cex = 0.8, col = cols, lwd = 2, bty = "n", xpd = T)

dev.off()
par(opar)


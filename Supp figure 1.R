## R-4.1.3 (RStudio-2023.06.0)
library(plyr)                          # version 1.8.7
library(dplyr)                         # version 1.0.9

rm(list = ls())

source("Functions for indices calculation.R")

#### Load data ####
clean.HD <- read.csv("final working HD set.csv")
clean.HD$date_visit <- as.Date(clean.HD$date_visit)
clean.HD$date_death <- as.Date(clean.HD$date_death)

## variable list
selected.biomarkersHD <- c("wbc", "hb", "hct", "mch", "mchc", "mcv", "pltlt", "rdw", "k", "sodium", "rbc")

# Add date of last contact
dat <- ddply(clean.HD, .(id_no), lastvis)

indices <- c("CVPC1", "MMD", "Df", "MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA", "MAF_var", "NMV", "Av_Var", 
             "PC_var", "Max_cov", "Ex_var", "Av_ab_cc")

# Load indices calculated for the main text to get the "nb_vis" variable
dat0 <- read.csv("Indices.csv")

dat.i <- list()
for (i in 1:2) {
  
  data <- dat
  
  if (i==2) {
    data$pltlt <- sqrt(data$pltlt)
    data[ , which(colnames(data) %in% c("wbc", "rdw", "Glucose"))] <-
      log(data[ , which(colnames(data) %in% c("wbc", "rdw", "Glucose"))])
  }
  
  data[ , selected.biomarkersHD] <- apply(data[ , selected.biomarkersHD], 2, scale)
  
  # Calculate Av_Var
  dat.av.var <- ddply(data, .(id_no), index.tw, vars = selected.biomarkersHD,
                      time.window = 6, min.vis.nb = 4, index.nm = "Av_Var", index.fct = avg.var)
  dat1 <- merge(dat0, dat.av.var, by = c("id_no", "years"), all = T)
  
  # Calculate Av_Ac
  dat.av.ac <- ddply(data, .(id_no), index.tw, vars = selected.biomarkersHD,
                     time.window = 6, min.vis.nb = 4, index.nm = "Av_Ac", index.fct = avg.autocorr)
  dat2 <- merge(dat1, dat.av.ac, by = c("id_no", "years"), all = T)
  
  # Calculate Df
  pca <- prcomp(data[ , selected.biomarkersHD])
  data$PC1 <- pca$x[ , 1]
  dat.df <- ddply(data, .(id_no), index.tw, vars = selected.biomarkersHD,
                  time.window = 6, min.vis.nb = 4, index.nm = "Df", index.fct = df.fct)
  dat3 <- merge(dat2, dat.df, by = c("id_no", "years"), all = T)
  
  # Calculate Av_ab_cc
  dat.aacc <- ddply(data, .(id_no), index.tw, vars = selected.biomarkersHD,
                    time.window = 6, min.vis.nb = 4, index.nm = "Av_ab_cc", index.fct = avg.abs.cc)
  dat4 <- merge(dat3, dat.aacc, by = c("id_no", "years"), all = T)
  
  dat4$version <- ifelse(i==1, "raw", "trsf")
  
  ## Apply transformations on indices
  dat4$Av_Var <- log(dat4$Av_Var)
  dat4$Av_ab_cc <- residuals(nls(Av_ab_cc ~ (a/nb_vis) + b, data = dat4, start = list(a = 1, b = 100)))
  dat4[ , c("Av_Ac", "Df")] <- apply(dat4[ , c("Av_Ac", "Df")], 2, function(x) {
    residuals(nls(x ~ I(1 / sqrt(nb_vis) * a) + b * sqrt(nb_vis), data = dat4, start = list(a = 1, b = 1)))
  }) 
  dat4$Av_Var <- residuals(lm(Av_Var ~ nb_vis, data = dat4))
  
  dat.i[[i]] <- dat4
}
datt <- do.call(rbind, dat.i)

### Plot
myPal <- c("#d81b60", "#1e88e5")

# Manual jitter
xpos1 <- c(-0.03, 0.03)

txts <- c(0.75, 0.67)

tiff("Fig_S1.tiff", width = 198, height = 165, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(2, 2), mar = c(2.5, 2.5, 0.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), mgp = c(1.5, 0.4, 0), tcl = -0.3)

indices <- c("Av_Var", "Av_Ac", "Df", "Av_ab_cc")

for (k in 1:length(indices)) {
  
  plot(1:10, 1:10, type = "n", xlim = c(5, 0), ylim = c(-0.4, 0.9), frame.plot = F,
       xlab = "Years before death", ylab = indices[k], font.lab = 2, cex.lab = 1, cex.axis = 0.8)
  mtext(c("A", "B", "C", "D")[k], side = 1, line = -16.25, at = 5.9, cex = 1.2, font = 2, xpd = T, adj = 0)
  
  datt0 <- datt[ , c("id_no", "years", "version", "date_death", "status", "age_visit", "sex", "diabetes", "fu_length", 
                     "nb_vis", indices[k])]
  dat.s <- split(datt0, datt0$version)
  
  for (i in 1:length(dat.s)) {
    
    datt1 <- dat.s[[i]]
    # Scale indices
    datt1[ , indices[k]] <- scale(datt1[ , indices[k]])
    
    dat.5yr <- datt1[datt1$years <= (5 * (12/6)) & !is.na(datt1$date_death), ]
    dat.5yr$time <- dat.5yr$years / (12/6)
    time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
    
    # Calculate mean and 95%CI per time point for each variable
    means_raw <- tapply(dat.5yr[ , indices[k]], dat.5yr$time, function(x) mean(x, na.rm = T))
    cis <- tapply(dat.5yr[ , indices[k]], dat.5yr$time, 
                  function(x) qnorm(0.975) * sd(x, na.rm = T) / sqrt(length(x)))
    # Center means so that values at 5 years start at 0
    means <- means_raw - means_raw[which(names(means_raw)==5)]
    
    points(x = time + xpos1[i], y = means, pch = c(16:17)[i], col = myPal[i], cex = 1.2)
    segments(x0 = time + xpos1[i], y0 = means + cis, x1 = time + xpos1[i], y1 = means - cis, col = myPal[i])
    segments(x0 = time + xpos1[i] - 0.02, y0 = means + cis, 
             x1 = time + xpos1[i] + 0.02, y1 = means + cis, col = myPal[i])
    segments(x0 = time + xpos1[i] - 0.02, y0 = means - cis, 
             x1 = time + xpos1[i] + 0.02, y1 = means - cis, col = myPal[i])
    for (j in 1:(length(time) - 1)) {
      segments(x0 = time[j] + xpos1[i], y0 = means[j], x1 = time[j + 1] + xpos1[i], 
               y1 = means[j + 1], lty = 3, col = myPal[i])
    }
    
    # Add HR95
    dat_var <- na.omit(datt1[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                                  "nb_vis", indices[k])])
    colnames(dat_var)[ncol(dat_var)] <- "index"
    mod <- coxph(Surv(years, status) ~ index + bs(age_visit, df = 5) + sex +  
                   diabetes + fu_length, data = dat_var, weights = sqrt(nb_vis),
                 cluster = id_no)
    hr95 <- summary(mod)$conf.int[1, c(1, 3, 4)] ^ 
      (quantile(na.omit(dat_var$index), 0.975) - quantile(na.omit(dat_var$index), 0.025))
    text(x = 3, y = txts[i],  col = myPal[i], font = 2, adj = 0, cex = 0.9, 
         labels = paste("HR95 = ", round(hr95[1], 2), " [", round(hr95[2], 2), ", ", 
                        round(hr95[3], 2), "]", sep = ""))
    
  }
  
  # Add legend
  legend(x = 4.9, y = 0.75, legend = c("Raw", "Transformed"), col = myPal, pch = 16:17, bty = "n", cex = 1, pt.cex = 1.2)
  text(x = 4.8, y = 0.77, labels = "Variables", adj = 0, font = 2, cex = 1)
  
}

dev.off()
par(opar)

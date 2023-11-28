## R-4.3.1 (RStudio-2023.06.0)
library(plyr)                          # version 1.8.8 
library(dplyr)                         # version 1.1.2
library(ggplot2)                       # version 3.4.2
library(RColorBrewer)                  # version 1.1-3
library(cowplot)                       # version 1.1.1

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
  
  # Calculate Node maximum variance
  dat.nmv <- ddply(data, .(id_no), index.tw, vars = selected.biomarkersHD,
                   time.window = 6, min.vis.nb = 4, index.nm = "NMV", index.fct = node.max.var)
  dat1 <- merge(dat0, dat.nmv, by = c("id_no", "years"), all = T)
  
  dat1$version <- ifelse(i==1, "Raw", "Transformed")
  
  ## Apply transformation and correct for number of observations
  dat1$NMV <- sqrt(dat1$NMV)
  dat1$NMV <- residuals(lm(NMV ~ nb_vis, data = dat1))
  
  dat.i[[i]] <- dat1
}
dat <- do.call(rbind, dat.i)

panelA <- ggplot(dat, aes(x=version, y=NMV, fill=version)) +
  geom_violin(width=2.1, size=0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position="none", axis.text.y=element_text(size=12, face="bold"),
        axis.text.x=element_text(size=9), axis.title=element_text(size=12)) +
  coord_flip() + 
  xlab("") +
  ylab("NMV (z-scores)")


##### 2. Trend before death and mortality prediction
cnms <- c("HR", "LCI", "UCI", "p", "HR95", "HR95_LCI", "HR95_UCI")
res <- as.data.frame(matrix(NA, 2, length(cnms)))
colnames(res) <- cnms
dat_raw <- na.omit(dat.i[[1]][ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                                   "nb_vis", "NMV")])
dat_trsf <- na.omit(dat.i[[2]][ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                                    "nb_vis", "NMV")])
fmla <- formula("Surv(years, status) ~ NMV + bs(age_visit, df = 5) + sex + diabetes + fu_length")
mod_raw <- coxph(fmla, data = dat_raw, weights = sqrt(nb_vis), cluster = id_no)
mod_trsf <- coxph(fmla, data = dat_trsf, weights = sqrt(nb_vis), cluster = id_no)
res[2, c("HR", "LCI", "UCI")] <- summary(mod_raw)$conf.int[1, c(1, 3, 4)]
res[2, "p"] <- summary(mod_raw)$coefficients[1, 6]
res[2, c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod_raw)$conf.int[1, c(1, 3, 4)] ^ 
  (quantile(na.omit(dat_raw[ , "NMV"]), 0.975) - quantile(na.omit(dat_raw[ , "NMV"]), 0.025))
res[1, c("HR", "LCI", "UCI")] <- summary(mod_trsf)$conf.int[1, c(1, 3, 4)]
res[1, "p"] <- summary(mod_trsf)$coefficients[1, 6]
res[1, c("HR95", "HR95_LCI", "HR95_UCI")] <- summary(mod_trsf)$conf.int[1, c(1, 3, 4)] ^ 
  (quantile(na.omit(dat_trsf[ , "NMV"]), 0.975) - quantile(na.omit(dat_trsf[ , "NMV"]), 0.025))

### Plot trend before death
# Manual jitter
xpos1 <- c(-0.03, 0.03)

myPal <- rev(brewer.pal(3, "Set1")[-3])

opar <- par(mar = c(2.5, 2.5, 0.5, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), mgp = c(1.5, 0.4, 0), tcl = -0.3)
plot(1:10, 1:10, type = "n", xlim = c(5, 0), ylim = c(-0.4, 0.9), frame.plot = F,
     xlab = "Years before death", ylab = "NMV (z-scores)", font.lab = 2, cex.lab = 1, cex.axis = 0.8)
segments(x0 = 5, y0 = -0.45, x1 = 0, y1 = -0.45)

for (i in 1:2) {
  if (i==1) dat <- dat.i[[2]] else dat <- dat.i[[1]]
  
  # Scale indices
  dat$NMV <- scale(dat$NMV)
  
  dat.5yr <- dat[dat$years <= (5 * (12/6)) & !is.na(dat$date_death), ]
  dat.5yr$time <- dat.5yr$years / (12/6)
  time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]
  
  # Calculate mean and 95%CI per time point for each variable
  means_raw <- tapply(dat.5yr$NMV, dat.5yr$time, function(x) mean(x, na.rm = T))
  cis <- tapply(dat.5yr$NMV, dat.5yr$time, 
                function(x) qnorm(0.975) * sd(x, na.rm = T) / sqrt(length(x)))
  # Center means so that values at 5 years start at 0
  means <- means_raw - means_raw[which(names(means_raw)==5)]
  
  points(x = time + xpos1[i], y = means, pch = c(16, 17)[i], col = myPal[i], cex = 1.1)
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
  text(x = 3.5, y = c(0.37, 0.28)[i], col = myPal[i], font = 2, adj = 0, cex = 0.9, 
       labels = paste("HR95 = ", round(res$HR95[i], 2), " [", round(res$HR95_LCI[i], 2), ", ", 
                      round(res$HR95_UCI[i], 2), "]", sep = ""))
}

# Add legend
legend(x = 4.9, y = 0.75, legend = c("Transformed", "Raw"), col = myPal, cex = 0.9, pch = 16:17, bty = "n", pt.cex = 1.1,
       y.intersp = 0.7)
text(x = 4.8, y = 0.76, labels = "Variables", cex = 1.1, adj = 0, font = 2)
panelB <- recordPlot()

tiff("Fig_S2.tiff", width = 300, height = 125, "mm", res = 600, compression = "lzw")
plot_grid(plotlist = list(panelA, panelB), ncol = 2, labels = c("A", "B"), label_size = 16, scale = c(0.95, 0.95))
dev.off()
par(opar)








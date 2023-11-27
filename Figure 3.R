## R-4.3.1 (RStudio-2023.06.0)
library(plyr)                          # version 1.8.8 
library(dplyr)                         # version 1.1.2
library(ggplot2)                       # version 3.4.2
#library(pals)                          # version 1.7
library(survival)                      # version 3.5-5
library(splines)                       # version 4.3.1
library(RColorBrewer)                  # version 1.1-3

rm(list = ls())

### Load indices 
dat <- read.csv("Indices.csv")

indices <- c("CVPC1", "MMD", "Df", "MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA", "MAF_var", "NMV", "Av_Var", 
             "PC_var", "Max_cov", "Ex_var", "Av_ab_cc")

### Change point analyses
# Use 3-month versions (as in the iScience article)
dats <- readRDS("HD indices_diff time windows.rds")
dat3 <- dats$mths3
# Cut time at 5 years before death
dat.5yr <- dat3[dat3$years <= (5 * (12/3)) & !is.na(dat3$date_death), ]
dat.5yr$time <- dat.5yr$years / (12/3)

cp.mods <- list()
for (var in 1:length(indices)) {
  print(paste("indice", var, sep = " = "))
  model = list(
    formula(paste(indices[var], "~ 1 + time", sep = " ")),
    1 + (1|id_no) ~ 0 + time
  )
  dat.ordered <- dat.5yr[order(dat.5yr$time, decreasing = T), ]
  cp.mods[[var]] <- mcp(model, data = dat.ordered)
  
}


tiff("Fig_3.tiff", width = 230, height = 90, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2, 8), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 0.7)

# Scale indices
dat[ , indices] <- apply(dat[ , indices], 2, scale)

# Set color palette
pal1 <- brewer.blues(7)[-1]
pal2 <- brewer.purples(4)[-1]
pal3 <- brewer.reds(6)[-1]
pal4 <- "#00955e"
myPal <- c(pal4, pal1, pal2, pal3)

dat.5yr <- dat[dat$years <= (5 * (12/6)) & !is.na(dat$date_death), ]
dat.5yr$time <- dat.5yr$years / (12/6)

time <- unique(dat.5yr$time)[order(unique(dat.5yr$time))]

ylims <- c(-0.4, 0.9)

# Manual jitter
xpos1 <- c(-0.09, -0.06, -0.03, 0, 0.03, 0.06, 0.09, 
           -0.105, -0.075, -0.045, -0.015, -0.015, 0.045, 0.075, 0.105)

pchs <- c(16, 15, 17, 18, 8, 11, 13, 16, 15, 17, 18, 8, 11, 13, 20)

# Change order by categories of indices
ind.cat2 <- list("Cross-correlation" =  "Av_ab_cc",
                 "Variance" = c("CVPC1", "NMV", "Av_Var", "PC_var", "Max_cov", "Ex_var"),
                 "Var + AC" = c("MMD", "Df", "MAF_var"),
                 "Auto-correlation" = c("MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA"))

plot(1:10, 1:10, type = "n", xlim = rev(range(time)), ylim = ylims, frame.plot = F,
     xlab = "Years before death", ylab = "Index (z-scores)", font.lab = 2)

for (i in 1:7) {
  # Calculate mean and 95%CI per time point for each variable
  means_raw <- tapply(dat.5yr[ , unlist(ind.cat2)[i]], dat.5yr$time, function(x) mean(x, na.rm = T))
  cis <- tapply(dat.5yr[ , unlist(ind.cat2)[i]], dat.5yr$time, 
                function(x) qnorm(0.975) * sd(x, na.rm = T) / sqrt(length(x)))
  # Center means so that values at 5 years start at 0
  means <- means_raw - means_raw[which(names(means_raw)==5)]
  
  points(x = time + xpos1[i], y = means, pch = pchs[i], col = myPal[i], cex = 1.2)
  segments(x0 = time + xpos1[i], y0 = means + cis, x1 = time + xpos1[i], y1 = means - cis, col = myPal[i])
  segments(x0 = time + xpos1[i] - 0.02, y0 = means + cis, 
           x1 = time + xpos1[i] + 0.02, y1 = means + cis, col = myPal[i])
  segments(x0 = time + xpos1[i] - 0.02, y0 = means - cis, 
           x1 = time + xpos1[i] + 0.02, y1 = means - cis, col = myPal[i])
  for (j in 1:(length(time) - 1)) {
    segments(x0 = time[j] + xpos1[i], y0 = means[j], x1 = time[j + 1] + xpos1[i], 
             y1 = means[j + 1], lty = 3, col = myPal[i])
  }
  
  # Add change point
  fit <- cp.mods[[which(indices==unlist(ind.cat2)[i])]]
  abline(v = summary(fit)$mean[1], col = myPal[i], lty = 2)
}

# Add legend for left and right panels
legend(x = -0.35, y = 0.9, legend = unlist(ind.cat2)[1:7], col = myPal[1:7], pch = pchs[1:7], lty = 3, xpd = T, 
       pt.cex = 1.3)
arrows(x0 = -0.8, x1 = -0.35, y0 = 0.94, y1 = 0.94, length = 0.08, lwd = 2, xpd = T)
legend(x = -0.54, y = 0.21, legend = unlist(ind.cat2)[8:15], col = myPal[8:15], pch = pchs[8:15], lty = 3, xpd = T, 
       pt.cex = 1.3)
arrows(x0 = -1.65, x1 = -2.1, y0 = 0.25, y1 = 0.25, length = 0.08, lwd = 2, xpd = T)

plot(1:10, 1:10, type = "n", xlim = rev(range(time)), ylim = ylims, frame.plot = F,
     xlab = "Years before death", ylab = "Index (z-scores)", font.lab = 2)

for (i in 8:15) {
  # Calculate mean and 95%CI per time point for each variable
  means_raw <- tapply(dat.5yr[ , unlist(ind.cat2)[i]], dat.5yr$time, function(x) mean(x, na.rm = T))
  cis <- tapply(dat.5yr[ , unlist(ind.cat2)[i]], dat.5yr$time, 
                function(x) qnorm(0.975) * sd(x, na.rm = T) / sqrt(length(x)))
  # Center means so that values at 5 years start at 0
  means <- means_raw - means_raw[which(names(means_raw)==5)]
  
  points(x = time + xpos1[i], y = means, pch = pchs[i], col = myPal[i], cex = 1.2)
  segments(x0 = time + xpos1[i], y0 = means + cis, 
           x1 = time + xpos1[i], y1 = means - cis, col = myPal[i])
  segments(x0 = time + xpos1[i] - 0.02, y0 = means + cis, 
           x1 = time + xpos1[i] + 0.02, y1 = means + cis, col = myPal[i])
  segments(x0 = time + xpos1[i] - 0.02, y0 = means - cis, 
           x1 = time + xpos1[i] + 0.02, y1 = means - cis, col = myPal[i])
  for (j in 1:(length(time) - 1)) {
    segments(x0 = time[j] + xpos1[i], y0 = means[j], x1 = time[j + 1] + xpos1[i], y1 = means[j + 1], lty = 3, 
             col = myPal[i])
  }
  
  # Add change point
  fit <- cp.mods[[which(indices==unlist(ind.cat2)[i])]]
  abline(v = summary(fit)$mean[1], col = myPal[i], lty = 2)
  
}

dev.off()
par(opar)


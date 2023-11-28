## R-4.3.1 (RStudio-2023.06.0)
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

indices <- c("CVPC1", "MMD", "Df", "MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA", "MAF_var", "NMV", "Av_Var", 
             "PC_var", "Max_cov", "Ex_var", "Av_ab_cc")
# Order by categories of indices (variance-based, etc.)
ind.cat <- list("Variance" = c("CVPC1", "NMV", "Av_Var", "PC_var", "Max_cov", "Ex_var"),
                "Var + AC" = c("MMD", "Df", "MAF_var"),
                "Auto-correlation" = c("MAF_ac", "MAF_ev", "MI", "Av_Ac", "NMA"),
                "Cross-correlation" =  "Av_ab_cc")


# Put all time windows together
for (j in 1:length(dats2)) dats2[[j]]$tw <- c("2 months", "3 months", "4 months", "6 months", "12 months")[j]
for (j in 1:length(dats2)) {
  mths <- c(2, 3, 4, 6, 12)[j]
  dat <- dats2[[j]]
  dat <- dat[dat$years <= (5 * (12/mths)) & !is.na(dat$date_death), ]   # Cut at 5 years before death
  dat$time <- dat$years / (12/mths)
  dats2[[j]] <- dat[ , c("id_no", "years", "time", "tw", "date_death", indices)]
} 
dat <- do.call(rbind, dats2)
dat$tw <- factor(dat$tw, levels = c("2 months", "3 months", "4 months", "6 months", "12 months"))

# Scale indices
dat[ , indices] <- apply(dat[ , indices], 2, scale)

# Set color palette
myPal <- rev(viridis(6)[-6])
pchs <- c(15, 17:20)

tiff("Fig_S5.tiff", width = 240, height = 300, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(5, 4), mar = c(2.3, 2.3, 0.2, 0.5), oma = c(0.2, 0.2, 0.2, 0.2), 
            mgp = c(1.5, 0.4, 0), tcl = -0.3, cex = 0.7)
for (k in 1:length(indices)) {
  
  if (k %in% 9:10) ylims <- c(-2.5, 2.5) else ylims <- c(-0.7, 1.2)
  
  plot(1:10, 1:10, type = "n", xlim = rev(range(dat$time)), ylim = ylims, frame.plot = F, 
       xlab = "Years before death", ylab = paste(unlist(ind.cat)[k], "(z-scores)"), font.lab = 2)
  
  for (i in 1:(length(dats2))) {
    dat.tw <- dat[dat$tw==levels(dat$tw)[i], ]
    time <- unique(dat.tw$time)[order(unique(dat.tw$time))]
    # Calculate mean and 95%CI per time point for each variable
    means_raw <- tapply(dat.tw[ , unlist(ind.cat)[k]], dat.tw$time, function(x) mean(x, na.rm = T))
    cis <- tapply(dat.tw[ , unlist(ind.cat)[k]], dat.tw$time, 
                  function(x) qnorm(0.975) * sd(x, na.rm = T) / sqrt(length(x)))
    # Center means so that values at 5 years start at 0
    means <- means_raw - means_raw[which(names(means_raw)==5)]
    
    points(x = time, y = means, pch = pchs[i], col = myPal[i], cex = 1.4)
    segments(x0 = time, y0 = means + cis, x1 = time, y1 = means - cis, col = myPal[i])
    segments(x0 = time - 0.02, y0 = means + cis, x1 = time + 0.02, y1 = means + cis, col = myPal[i])
    segments(x0 = time - 0.02, y0 = means - cis, x1 = time + 0.02, y1 = means - cis, col = myPal[i])
    for (j in 1:(length(time) - 1)) {
      segments(x0 = time[j], y0 = means[j], x1 = time[j + 1], y1 = means[j + 1], lty = 3, col = myPal[i])
    }
  }
  
  # Plot legend 
  if (k %in% c(3, 6, 9, 12)) {
    plot(1:10, 1:10, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
    if (k==3) {
      legend(x = -1, y = 9, legend = levels(dat$tw), col = myPal, pch = pchs, lty = 3, xpd = T, pt.cex = 1.4, bty = "n")
      text(x = -1, y = 9.3, labels = "Time windows", cex = 1.2, font = 2, adj = 0, xpd = T)
    } 
  }
}
dev.off()
par(opar)


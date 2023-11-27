## R-4.3.1 (RStudio-2023.06.0)
library(corrplot)                      # version 0.92
library(corrr)                         # version 0.4.4

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


# Define colors for categories of indices
cols <- c("#0097f2", "#804C9A", "#ff0041", "#00955e")
col_cat = c(rep(cols[1], length(ind.cat[[1]])), rep(cols[2], length(ind.cat[[2]])), 
            rep(cols[3], length(ind.cat[[3]])), rep(cols[4], length(ind.cat[[4]])))

tiff("Fig_2.tiff", width = 180, height = 140, "mm", res = 300, compression = "lzw")
corrplot.mixed(cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs"), lower = "number", upper = "ellipse",
               p.mat = cor.mtest(dat[ , unlist(ind.cat)])$p, tl.pos = "lt", tl.cex = 0.8, pch.cex = 0.9,
               number.cex = 0.7, cl.cex = 0.6, cl.ratio = 0.1, mar = c(0.5, 0.5, 1.5, 7.5),
               tl.col = col_cat)

segments(x0 = 0.5, x1 = 9.5, y0 = 15.5, y1 = 15.5, col = cols[1], lwd = 2)
segments(x0 = 0.5, x1 = 0.5, y0 = 15.5, y1 = 6.5, col = cols[1], lwd = 2)
segments(x0 = 9.5, x1 = 9.5, y0 = 15.5, y1 = 6.5, col = cols[1], lwd = 2)
segments(x0 = 0.5, x1 = 9.5, y0 = 6.5, y1 = 6.5, col = cols[1], lwd = 2)

segments(x0 = 6.5, x1 = 14.5, y0 = 9.5, y1 = 9.5, col = cols[3], lwd = 2)
segments(x0 = 14.5, x1 = 14.5, y0 = 1.5, y1 = 9.5, col = cols[3], lwd = 2)
segments(x0 = 6.5, x1 = 6.5, y0 = 1.5, y1 = 9.5, col = cols[3], lwd = 2)
segments(x0 = 6.5, x1 = 14.5, y0 = 1.5, y1 = 1.5, col = cols[3], lwd = 2)

segments(x0 = 15.5, x1 = 14.5, y0 = 1.5, y1 = 1.5, col = cols[4], lwd = 2)
segments(x0 = 14.5, x1 = 14.5, y0 = 1.5, y1 = 0.5, col = cols[4], lwd = 2)
segments(x0 = 15.5, x1 = 15.5, y0 = 1.5, y1 = 0.5, col = cols[4], lwd = 2)
segments(x0 = 15.5, x1 = 14.5, y0 = 0.5, y1 = 0.5, col = cols[4], lwd = 2)

legend(x = 17, y = 16, legend = names(ind.cat), col = cols, lwd = 2, xpd = T, cex = 0.8, bty = "n")

dev.off()

### Get correlations among categories
var.ind <- cor(dat[ , ind.cat$Variance], use = "pairwise.complete.obs")
ac.ind <- cor(dat[ , ind.cat$Autocorrelation], use = "pairwise.complete.obs")
all <- cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")

mean(var.ind[lower.tri(var.ind)])
sd(var.ind[lower.tri(var.ind)])
mean(ac.ind[lower.tri(ac.ind)])
sd(ac.ind[lower.tri(ac.ind)])
mean(all[lower.tri(all)])
sd(all[lower.tri(all)])
cor.mmd <- cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")[ , "MMD"]
mean(cor.mmd[which(names(cor.mmd) %in% ind.cat$Variance)])
cor.mtest(dat[ , unlist(ind.cat)])$p
cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")[ind.cat$Variance, c("Df", "MAF_var")]
cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")[ind.cat$Autocorrelation, c("Df", "MAF_var")]
cor.df <- cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")[ , "Df"]
mean(cor.df[-which(names(cor.df) %in% c("Df", "Av_Var"))])
cor.maf <- cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")[ , "MAF_var"]
mean(cor.maf[-which(names(cor.maf)=="MAF_var")])
cor.cc <- cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs")[ , "Av_ab_cc"]
mean(cor.cc[-which(names(cor.cc)=="Av_ab_cc")])


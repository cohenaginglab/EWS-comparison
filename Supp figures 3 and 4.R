## R-4.3.1 (RStudio-2023.06.0)
library(ggplot2)                       # version 3.4.2
library(viridis)                       # version 0.6.3
library(corrplot)                      # version 0.92

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

#### 1. Summary plot of mean correlations by category of indices
# Extract correlation coefficients
get.cc <- function(x, category) {
  Mcor <- cor(x[ , category], use = "pairwise.complete.obs")
  Mcor[lower.tri(Mcor, diag = F)]
  
}
get.cc2 <- function(x, cat1, cat2) {
  c(cor(x = x[ , cat1], y = x[ , cat2], use = "pairwise.complete.obs"))
}

cor.var <- lapply(dats2, get.cc, category = ind.cat$Variance)
cor.var_AC <- lapply(dats2, get.cc, category = ind.cat$`Var + AC`)
cor.AC <- lapply(dats2, get.cc, category = ind.cat$`Auto-correlation`)
cor.varVSAC <- lapply(dats2, get.cc2, cat1 = ind.cat$Variance, cat2 = ind.cat$`Auto-correlation`)
cor.varVSCC <- lapply(dats2, get.cc2, cat1 = ind.cat$Variance, cat2 = ind.cat$`Cross-correlation`)
cor.ACvsCC <- lapply(dats2, get.cc2, cat1 = ind.cat$`Auto-correlation`, cat2 = ind.cat$`Cross-correlation`)
res <- as.data.frame(matrix(NA, 
                            (length(cor.var[[1]]) + length(cor.var_AC[[1]]) + length(cor.AC[[1]]) + 
                               length(cor.varVSAC[[1]]) + length(cor.varVSCC[[1]]) + length(cor.ACvsCC[[1]])) * 5,
                            3))
colnames(res) <- c("compare", "tw", "value")
cor.list <- list(cor.var, cor.var_AC, cor.AC, cor.varVSAC, cor.varVSCC, cor.ACvsCC)
compare_cat <- c("Variance", "Var + AC", "Auto-correlation", "Var vs AC", "Var vs CC", "AC vs CC")
for (k in 1:length(cor.list)) {
  cor.res <- cor.list[[k]]
  nb <- length(cor.res[[1]])
  start_row <- if (k==1) 1 else stop_row + 1
  stop_row <- if (k==1) nb * 5 else stop_row + nb * 5
  res[start_row:stop_row, "value"] <- unlist(cor.res)
  res[start_row:stop_row, "compare"] <- rep(compare_cat[k], nb * 5)
  res[start_row:stop_row, "tw"] <- c(rep(2, nb), rep(3, nb), rep(4, nb), rep(6, nb), rep(12, nb)) 
}
res$compare <- factor(factor(res$compare, levels = compare_cat), 
                      labels = c("Var", "Var + AC", "AC", "Var vs AC", "Var vs CC", "AC vs CC"))
res$tw <- factor(res$tw, labels = paste(levels(factor(res$tw)), "months"))

myPal <- viridis(6)[-1]

p <- ggplot(res, aes(x = compare, y = value, fill = tw)) + 
  geom_boxplot() +
  scale_fill_manual(values = myPal, name = "Time window") +
  ylab("Correlation coefficient") +
  xlab("") +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.title = element_text(size = 12), axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 12))

tiff("Fig_S3.tiff", width = 190, height = 110, "mm", res = 600, compression = "lzw")
print(p)
dev.off()


#### 2. Check correlations among indices
# Define colors for categories of indices
cols <- c("#0097f2", "#804C9A", "#ff0041", "#00955e")
col_cat = c(rep(cols[1], length(ind.cat[[1]])), rep(cols[2], length(ind.cat[[2]])), 
            rep(cols[3], length(ind.cat[[3]])), rep(cols[4], length(ind.cat[[4]])))

tiff("Fig_S4.tiff", width = 210, height = 200, "mm", res = 600, compression = "lzw")
opar <- par(mfrow = c(2, 2))

for (tw in 1:length(dats)) {
  dat <- dats[[tw]]
  corrplot.mixed(cor(dat[ , unlist(ind.cat)], use = "pairwise.complete.obs"), lower = "number", upper = "ellipse",
                 p.mat = cor.mtest(dat[ , unlist(ind.cat)])$p, tl.pos = "lt", tl.cex = 0.8, pch.cex = 0.9,
                 number.cex = 0.7, cl.cex = 0.6, cl.ratio = 0.1, tl.col = col_cat)
  
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
  
  mtext(c("A", "B", "C", "D")[tw], side = 1, line = -17, at = -2, adj = 0, cex = 1.3, font = 2)
}

dev.off()
par(opar)


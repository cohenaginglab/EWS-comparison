## R-4.3.1 (RStudio-2023.06.0)
library(plyr)                          # version 1.8.8 
library(dplyr)                         # version 1.1.2
library(ggplot2)                       # version 3.4.2
library(pROC)                          # version 1.18.2
library(survival)                      # version 3.5-5
library(splines)                       # version 4.3.1
library(corrr)                         # version 0.4.4
library(cowplot)                       # version 1.1.1

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

# Set color palette
cols <- c("#0097f2", "#804C9A", "#ff0041", "#00955e")
pal1 <- rep(cols[1], length(ind.cat$Variance))
pal2 <- rep(cols[2], length(ind.cat$`Var + AC`))
pal3 <- rep(cols[3], length(ind.cat$Autocorrelation))
pal4 <- cols[4]
myPal <- c(pal1, pal2, pal3, pal4)
names(myPal) <- unlist(ind.cat)

hc.plots <- list()
auc.plots <- list()
for (k in 1:2) {
  if (k==2) indices2 <- setdiff(indices, c("MAF_ac", "MAF_var")) else indices2 <- indices
  if (k==2) myPal2 <- myPal[-which(names(myPal) %in% c("MAF_ac", "MAF_var"))] else myPal2 <- myPal
  
  datt <- na.omit(dat[ , c("id_no", "years", "status", "age_visit", "sex", "diabetes", "fu_length", 
                           "nb_vis", indices2)])
  
  # Re-order indices based on hierarchical clustering
  p <- correlate(datt[ , indices2]) %>% rearrange(absolute = FALSE) %>% rplot()
  hc.plots[[k]] <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  new.indices <- colnames(rearrange(correlate(datt[ , indices2])))[-1]
  
  res_j <- list()
  for (j in 1:length(indices2)) {
    res <- as.data.frame(matrix(0, 1, length(new.indices) + 2))
    colnames(res) <- c(new.indices, "n", "AUC")
    res[ , 1:j] <- 1
    
    res$n <- length(unique(datt$id_no))
    fmla <- formula(paste("Surv(years, status) ~ ", 
                          paste(new.indices[1:j], collapse = " + "), 
                          "+ bs(age_visit, df = 5) + sex + diabetes + fu_length",
                          sep = ""))
    
    # Run cox model
    cox.mod <- coxph(fmla, data = datt, weights = sqrt(nb_vis), cluster = id_no)
    roc.mod <- roc(datt$status ~ predict(cox.mod))
    res$AUC <- as.numeric(auc(roc.mod))
    res_j[[j]] <- res
  }
  res1 <- bind_rows(res_j)
  
  # Run model with only control variables
  res <- as.data.frame(matrix(0, 1, length(new.indices) + 2))
  colnames(res) <- c(new.indices, "n", "AUC")
  res$n <- length(unique(datt$id_no))
  cox.mod <- coxph(Surv(years, status) ~ bs(age_visit, df = 5) + sex + diabetes + fu_length, 
                   data = datt, weights = sqrt(nb_vis), cluster = id_no)
  roc.mod <- roc(datt$status ~ predict(cox.mod))
  res$AUC <- as.numeric(auc(roc.mod))
  res1 <- rbind(res, res1)
  
  par(mar = c(4.5, 3.5, 0, 0.2), mgp = c(2, 0.4, 0), tcl = -0.3, oma = c(0, 0, 0, 0))
  plot(1:nrow(res1), 1:nrow(res1), type = "n", xlab = "", ylim = c(0.71, 0.84), ylab = "AUC", frame.plot = F, 
       cex.axis = 0.75, xaxt = "n", font.lab = 2, cex.lab = 0.85)
  ## Add grid 
  abline(h = c(0.75, 0.8), col = "grey85", lwd = 2)
  abline(h = c(0.725, 0.775, 0.825), col = "grey85")
  
  segments(x0 = 1, y0 = 0.71, x1 = nrow(res1), y1 = 0.71, xpd = T)
  for (i in 1:nrow(res1)) segments(x0 = i, y0 = 0.71, x1 = i, y1 = 0.707, xpd = T)
  for (i in 1:nrow(res1)) segments(x0 = i, y0 = 0.711, x1 = i, y1 = 0.84, xpd = T, col = "grey85")
  text(x = 1:nrow(res1), y = 0.703, labels = c("Controls", new.indices), 
       col = c("black", myPal[match(new.indices, names(myPal))]),
       srt = 90, cex = 0.75, font = 2, adj = 1, xpd = T)
  
  lines(x = 1:nrow(res1), y = res1$AUC, lwd = 2) 
  # Add legend
  legend(x = 8, y = 0.78, legend = names(ind.cat), col = cols, lwd = 3, cex = 0.7,
         bg = "white", box.col = "white")
  auc.plots[[k]] <- recordPlot() 
}

# Plot
plots <- list(auc.plots[[1]], hc.plots[[1]], auc.plots[[2]], hc.plots[[2]])
tiff("Fig_S7.tiff", width = 210, height = 170, "mm", res = 300, compression = "lzw")
plot_grid(plotlist = plots, ncol = 2, labels = c("A", "B", "C", "D"), label_size = 13, scale = c(0.86, 0.95, 0.86, 0.95))
dev.off()



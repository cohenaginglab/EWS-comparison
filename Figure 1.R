## R-4.3.1 (RStudio-2023.06.0)
library(ggplot2)                       # version 3.4.2
library(tidyverse)                     # version 2.0.0
library(cowplot)                       # version 1.1.1
library(flextable)                     # version 0.9.2
library(grid)                          # version 4.3.1
library(gridGraphics)                  # version 0.5-1


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

##### 1. Plot distributions for each index
dat.ind <- list()
for (i in 1:length(indices)) {
  datt <- dat[ , c("id_no", "years", indices[i])]
  datt$index <- indices[i]
  colnames(datt)[3] <- "value"
  dat.ind[[i]] <- datt
}
dat.ind <- do.call(rbind, dat.ind)
dat.ind$index <- factor(dat.ind$index)

plot.var <- ggplot(dat.ind[dat.ind$index %in% ind.cat$Variance, ], aes(x = value)) + 
  geom_histogram() + facet_grid(. ~ index, scales = "free") +
  background_grid(major = 'y', minor = "none") +  
  panel_border() +
  theme(text = element_text(size = 8), strip.text.x = element_text(size = 11, face = "bold")) 

plot.ac <- ggplot(dat.ind[dat.ind$index %in% ind.cat$Autocorrelation, ], aes(x = value)) + 
  geom_histogram() + facet_grid(. ~ index, scales = "free") +
  background_grid(major = 'y', minor = "none") +  
  panel_border() +
  theme(text = element_text(size = 8), strip.text.x = element_text(size = 11, face = "bold")) 

plot.var_ac <- ggplot(dat.ind[dat.ind$index %in% ind.cat$`Var + AC`, ], aes(x = value)) + 
  geom_histogram() + facet_grid(. ~ index, scales = "free") +
  background_grid(major = 'y', minor = "none") +  
  panel_border() +
  theme(text = element_text(size = 8), strip.text.x = element_text(size = 11, face = "bold")) 

plot.cc <- ggplot(dat.ind[dat.ind$index %in% ind.cat$`Cross-correlation`, ], aes(x = value)) + 
  geom_histogram() + facet_grid(. ~ index, scales = "free") +
  background_grid(major = 'y', minor = "none") +  
  panel_border() +
  theme(text = element_text(size = 8), strip.text.x = element_text(size = 11, face = "bold"))

nms_ind <- tibble("Abbreviation" = unlist(ind.cat),
                  "Index" = c("First principal component of coefficients of variation",
                              "Node maximum variance",
                              "Average variance",
                              "Principal component analysis of variance",
                              "Maximum value of covariance matrix",
                              "Explained variance",
                              "Multivariate moving distance",
                              "Degenerate fingerprinting",
                              "MAF variance",
                              "Maximum autocorrelation factor (MAF) autocorrelation",
                              "MAF eigenvalue",
                              "Mutual information",
                              "Average autocorrelation",
                              "Node maximum autocorrelation", 
                              "Average absolute cross-correlation"))

ft_raster <- nms_ind[order(nms_ind$Abbreviation), ] %>% flextable::flextable() %>% 
  as_raster()
plot_tbl <- ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(ft_raster), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

p <- ggdraw() +
  theme(plot.background = element_rect(fill="white", color = NA), plot.margin = margin(0, 0, 0, 0)) +
  draw_plot(plot.var, 0, .64, 1, .32) +
  draw_plot(plot.ac, 0, .32, .7, .32) +
  draw_plot(plot.var_ac, 0, 0, .5, .32) +
  draw_plot(plot.cc, .5, 0, .2, .32) +
  draw_plot(plot_tbl, .66, 0.02, .35, .63) +
  draw_plot_label(c("Variance", "Autocorrelation", "Variance + Autocorrelation", "Cross-correlation"), 
                  x = c(0.01, 0.01, 0.01, .5), y = c(0.98, .66, .34, .34), size = 11, hjust = 0)

ggsave2("Fig_1.tiff", p, width = 180, height = 200, units = "mm", dpi = 600, compression = "lzw")



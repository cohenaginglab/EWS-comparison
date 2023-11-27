## R-4.1.1 (RStudio-2021.09.0-351)
require(plyr)                          # version 1.8.6 
require(dplyr)                         # version 1.0.7
require(maf)                           # version 0.0.0.9000
require(DescTools)                     # version 0.99.44
require(survival)                      # version 3.2-11
require(splines)                       # version 4.1.1


## Calculate date of last contact
lastvis <- function (dat.id) {
  dat.id$date_last <- rep(max(c(max(dat.id$date_visit), max(dat.id$date_death)), na.rm = T), nrow(dat.id))
  dat.id$fu_length <- as.numeric(dat.id$date_last[1] - min(dat.id$date_visit)) / 365.25
  return (dat.id)
}


## Add "status" variable
status_fct <- function(dat.id) {
  dat.id$status <- ifelse(!is.na(dat.id$date_death) & dat.id$years==min(dat.id$years), 1, 0)
  return(dat.id)
}


## Function to calculate index by selected time window
index.tw <- function(dat.id, vars = selected.biomarkersHD, min.vis.nb = 4, time.window = 6,
                     index.nm = "Df", index.fct = df.fct, mean = F) {
  
  dat.id$years <- floor(as.numeric(dat.id$date_last - dat.id$date_visit) / (365 / (12/time.window)))
  # remove years with insufficient visit number
  counts <- table(dat.id$years)
  dat.id2 <- dat.id[dat.id$years %in% names(counts[which(counts >= min.vis.nb)]), ]
  
  if (nrow(dat.id2) > 0) {
    # calculate index per time window
    ID_by_yr <- split(dat.id2, dat.id2$years)
    ID_cv <- as.data.frame(matrix(NA, length(ID_by_yr), 1))
    colnames(ID_cv) <- "years"
    ID_cv[ , "years"] <- as.numeric(names(ID_by_yr))
    
    if (mean==F) ID_cv[ , index.nm] <- unlist(lapply(ID_by_yr, index.fct)) else 
      ID_cv[ , index.nm] <- unlist(lapply(ID_by_yr, function(z) mean(z$MMD)))
    
    return(ID_cv)
  }
}


## Degenerative fingerprinting (on PC1)
df.fct <- function(x) {
  x <- x[order(x$date_visit), ]
  x_t0 <- x$PC1[1:(nrow(x) - 1)]
  x_t1 <- x$PC1[2:nrow(x)]
  return(cor(x_t0, x_t1))
}


## MAF autocorrelation
maf.autocorr <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(maf(as.matrix(x))$autocor[1])
}


## MAF eigenvalue
maf.eigenval <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(min(maf(as.matrix(x))$rotation[ , 1]))
}


## Mutual information
mut.inf <- function(x) {
  x <- x[order(x$date_visit), selected.biomarkersHD]
  x_t0 <- x[1:(nrow(x) - 1), ]
  x_t1 <- x[2:nrow(x), ]
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x_t0 <- apply(x_t0, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  x_t1 <- apply(x_t1, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(MutInf(as.matrix(x_t0), as.matrix(x_t1)))
}


## Average autocorrelation
avg.autocorr <- function(x) {
  x <- x[order(x$date_visit), selected.biomarkersHD]
  x_t0 <- x[1:(nrow(x) - 1), ]
  x_t1 <- x[2:nrow(x), ]
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x_t0 <- apply(x_t0, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  x_t1 <- apply(x_t1, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  corr <- c()
  for (i in 1:length(selected.biomarkersHD)) corr[i] <- cor(x_t0[ , i], x_t1[ , i])
  return(mean(corr))
}


## Node maximum autocorrelation
node.max.autocorr <- function(x) {
  x <- x[order(x$date_visit), selected.biomarkersHD]
  x_t0 <- x[1:(nrow(x) - 1), ]
  x_t1 <- x[2:nrow(x), ]
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x_t0 <- apply(x_t0, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  x_t1 <- apply(x_t1, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  corr <- c()
  for (i in 1:length(selected.biomarkersHD)) corr[i] <- cor(x_t0[ , i], x_t1[ , i])
  return(max(corr))
}


## MAF variance
maf.var <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(var(maf(as.matrix(x))$rotation[ , 1]))
}


## Node maximum variance
node.max.var <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(max(apply(x, 2, var)))
}


## Average variance 
avg.var <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(mean(apply(x, 2, var)))
}


## PCA variance
pca.var <- function(x) var(x$PC1)


## Maximum value of covariance matrix
max.cov <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  xx <- cov(x)
  return(max(abs(xx[lower.tri(xx)])))
}


## Explained variance
expl.var <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  return(max(prcomp(cov(x))$sdev^2) / sum(prcomp(cov(x))$sdev^2))
}


## Average absolute cross-correlation
avg.abs.cc <- function(x) {
  # If there's no variation in one biomarker across all measures (due to lab rounding), 
  # introduce random noise
  x <- x[ , selected.biomarkersHD]
  x <- apply(x, 2, function(z) if (length(unique(z))==1) z + rnorm(length(z), sd = 0.001) else z)
  corr <- cor(x)
  return(mean(abs(corr[upper.tri(corr)])))
}


## For CVPC1
cv.4mth <- function(dat.id, selected.biomarkersHD, nb.month = 6, min.vis.nb = 4, SDs = SDs,
                    covars = c("id_no", "sex", "date_death", "date_birth", "date_last", "diabetes", 
                               "fu_length")) {
  
  dat.id$years <- floor(as.numeric(dat.id$date_last - dat.id$date_visit) / (365 / (12/nb.month)))
  # remove years with insufficient visit number, otherwise CVPC1 is biased
  counts <- table(dat.id$years)
  dat.id2 <- dat.id[dat.id$years %in% names(counts[which(counts >= min.vis.nb)]), ]
  
  if (nrow(dat.id2) > 0) {
    # calculate CVs
    ID_by_yr <- split(dat.id2, dat.id2$years)
    ID_cv <- as.data.frame(matrix(NA, length(ID_by_yr), length(selected.biomarkersHD) + length(covars) + 3))
    colnames(ID_cv) <- c(covars, "age_visit", "years", "nb_vis", selected.biomarkersHD)
    ID_cv[ , covars] <- rdply(length(ID_by_yr), ID_by_yr[[1]][1, covars], .id = NULL)
    ID_cv[ , "age_visit"] <- unlist(lapply(ID_by_yr, function(x) min(x$age_visit)))
    ID_cv[ , "years"] <- as.numeric(names(ID_by_yr))
    ID_cv[ , "nb_vis"] <- unlist(lapply(ID_by_yr, nrow))
    
    # Calculate CVs 
    cvs <- lapply(ID_by_yr, function(x) { 
      cvs.yr <- c()
      for (j in 1:length(selected.biomarkersHD)) {
        z <- x[ , selected.biomarkersHD[j]]
        # If there's no variation in one biomarker across all measures (due to lab rounding), 
        # introduce random noise
        if (length(unique(z))==1) z <- z + rnorm(length(z), sd = SDs[j] / 1000) 
        cvs.yr[j] <- sd(z) / mean(z)
      }
      return(cvs.yr)
    })
    ID_cv[ , selected.biomarkersHD] <- t(bind_cols(cvs))
    return(as.data.frame(ID_cv))
  }
}


## Change sign of PC
pc_sign <- function(data, pc.name = "PC1", pca) {
  datt <- data
  colnames(datt)[which(colnames(datt)==pc.name)] <- "predictor"
  mod <- coxph(Surv(years, status) ~ predictor + bs(age_visit, df = 5) + sex +  
                 diabetes + fu_length, data = datt, weights = sqrt(nb_vis),
               cluster = id_no)
  hr <- as.numeric(summary(mod)$conf.int[1, 1])
  if (hr < 1) {
    PC <- -1 * datt$predictor
    datt.2 <- cbind(datt[ , -which(colnames(datt)=="predictor")], PC)
    colnames(datt.2)[which(colnames(datt.2)=="PC")] <- pc.name
    return(datt.2)
  } else return(data)
}


# Function to calculate MMD
mmd_calc = function(dat.id, covars, vars, time.max = NULL, ma = F, ma.time = 0.5) {
  if (nrow(dat.id) > 1) {
    
    dat.id = dat.id[order(dat.id$age_visit), ]
    
    if (!is.null(time.max)) {
      # Select visits that fall in the maximum time interval 
      dat.id$time.diff = c(0, unlist(Map(function(x, y) as.numeric(y - x), 
                                         dat.id$age_visit[-length(dat.id$age_visit)], dat.id$age_visit[-1]))) 
      gaps = which(dat.id$time.diff > time.max) 
      # Split data according to gaps in visits (as defined by time.max)
      if (length(gaps) > 0) {
        for (j in 1:(length(gaps) + 1)) {
          if (j==1) {
            dat.id[1:gaps[1], "period"] = 1
          } else {
            if (j==(length(gaps) + 1)) {
              dat.id[gaps[j - 1]:nrow(dat.id), "period"] = j
            } else {
              dat.id[gaps[j - 1]:gaps[j], "period"] = j
            }
          }
        }
      } else {
        dat.id$period = 1
      }
      # Calculate moving multivariate distance (MMD) by "period"
      dat.id$MMD = NA
      dat.id = ddply(dat.id, .(period), function(ID_by_period) {
        if (nrow(ID_by_period) > 1) {
          ref = as.numeric(ID_by_period[1, vars])
          for (i in 2:nrow(ID_by_period)) {
            ID_by_period$MMD[i] = sqrt(mahalanobis(as.matrix(ID_by_period[i, vars]), ref, covars))
            ref = as.numeric(ID_by_period[i, vars])
          }
          return(ID_by_period)
        }
      })
      
    } 
    if (ma==T) {
      # Calculate moving multivariate distance (MMD) using a moving window
      nr = which((dat.id$age_visit - ma.time) >= min(dat.id$age_visit))
      dat.id$MMD = NA
      for (k in nr) {
        ages = which(dat.id$age_visit >= (dat.id[k, "age_visit"] - ma.time) & dat.id$age_visit <= dat.id[k, "age_visit"])
        if (length(ages) > 1) {
          ref = dat.id[ages, vars]
          if (nrow(ref) > 1) means = colMeans(ref) else means = as.numeric(ref)
          dat.id$MMD[k] = sqrt(mahalanobis(x = as.matrix(dat.id[k, vars]), center = means, cov = covars))
        }
      }
    } else {
      # Calculate moving multivariate distance (MMD)
      dat.id$MMD = NA
      for (i in 2:nrow(dat.id)) {
        dat.id$MMD[i] = sqrt(mahalanobis(x = as.matrix(dat.id[i, vars]), 
                                         center = as.numeric(dat.id[i - 1, vars]), 
                                         cov = covars))
      }
    }
    return(dat.id)
  }
}


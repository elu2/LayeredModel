library(tidyverse)
library(glmnet)
library(zeallot)
library(parallel)

select = dplyr::select

source("LASSOFunctions.R")

m1d3 <- read_csv("dataOut/m1d3.csv")
m1d7 <- read_csv("dataOut/m1d7.csv")
m2d3 <- read_csv("dataOut/m2d3.csv")
m2d7 <- read_csv("dataOut/m2d7.csv")

fpath <- "dataOut/l2d7final_nzc.csv"
epath <- "dataOut/l2d7lmErrs.csv"

# Raw delta2 values
#dm1 <- (m1d7 - m1d3) %>% select(-Public, -day)
#dm2 <- (m2d7 - m2d3) %>% select(-Public, -day)

# Raw day7 values
#dm1 <- m1d7 %>% select(-Public, -day)
#dm2 <- m2d7 %>% select(-Public, -day)

# Log2 delta2 values
#dm1 <- (log2(m1d7) - log2(m1d3)) %>% select(-Public, -day)
#dm2 <- (log2(m2d7) - log2(m2d3)) %>% select(-Public, -day)

# Log2 day7 values
dm1 <- log2(m1d7) %>% select(-Public, -day)
dm2 <- log2(m2d7) %>% select(-Public, -day)

# Add y column
dm1 <- dm1 %>% mutate(y=0)
dm2 <- dm2 %>% mutate(y=1)


# Combine both sets
xy <- rbind(dm1, dm2)

x <- as.matrix(xy %>% select(-y, -`Acuity max`))
y <- as.matrix(xy %>% select(y))


# Run model for 1000 iterations
layers <- 3
layer <- seq(1, layers)
trials <- seq(1, 1000)


run_trial <- function(){
  c(final_coefs, oErrs, fnrs, fprs) %<-% layered_model(x, y, xy, alpha1=1, alphai=1, layers)
  
  err_df <- data.frame(oErrs, fnrs,fprs, trial, layer)
  df <- t(data.frame(final_coefs))
  
  write.table(err_df, epath, sep=",", col.names=!file.exists(epath), row.names=F, append=T)
  write.table(df, fpath, sep=",", col.names=!file.exists(fpath), row.names=F, append=T)
}

print("Running 1000 trials.")
print(paste("Outputting coefficients to", fpath))
print(paste("Outputting error rates to", epath))
mclapply(trials, run_trial, mc.cores = detectCores())
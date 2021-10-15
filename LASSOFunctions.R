library(glmnet)
library(tidyverse)

run <- function(x, y, alpha=1, foldid=NULL){
  fit <- glmnet(x, y, family="binomial", alpha=alpha)
  cv_fit <- cv.glmnet(as.matrix(x), as.matrix(y), family="binomial", type.measure="class", alpha=alpha, foldid=foldid, nfolds=10)
  return(list(fit, cv_fit))
}


# From an x and y, extract the cross-validated error
cv_error <- function(run_x, run_y){
  # Ridge-to-LASSO parameter vector
  alpha_values <- seq(0, 1, length.out=100)
  fold_ids <- sample(1:10, nrow(run_x), replace=T)
  # Vector to store errors
  min_error <- vector()
  
  i <- 1
  for (alpha in alpha_values){
    # Verbose runtime confirmation
    print(paste("Running alpha=", alpha))
    
    run_fits <- run(run_x, run_y, alpha=alpha, foldid=fold_ids)
    
    # Cross validated LASSO at specified alpha value
    cv.fit <- run_fits[[2]]
    
    # Index of a lambda value (lambda min in this case)
    index <- which(cv.fit$lambda == cv.fit$lambda.min)
    
    # ...to index the error value
    mse.min <- cv.fit$cvm[index]
    min_error[i] <- mse.min
    
    i <- i + 1
  }
  print(plot(alpha_values, min_error))
  return(min_error)
}


# neit - dataframe with neither bistable or rebi
# secondary - dataframe with either bistable or rebi
# n - number of total rows
# prop - proportion of neit to secondary data rows

enrich_combine <- function(neit, secondary, n_rows, prop){
  n1 <- ceiling(n_rows * prop)
  n2 <- ceiling(n_rows * (1 - prop))
  
  rNeitIndices <- sample(nrow(neit), n1)
  rBiSubsetIndices <- sample(nrow(secondary), n2)
  
  df_neit <- neit[rNeitIndices, ]
  df_bi <- secondary[rBiSubsetIndices, ]
  
  combined <- rbind(df_neit, df_bi)
  
  return(combined)
}


# Function that takes data and partitions out x, y, and the data unsampled from

train_set <- function(data, n) {
  last_col <- ncol(data)
  train_indices <- sample(nrow(data), n)
  
  tr_x <- as.matrix(data[train_indices, 1:(last_col-1)])
  tr_y <- as.matrix(rev(data)[train_indices, 1])
  
  remainder <- data[-train_indices,]
  
  return(list(tr_x, tr_y, remainder))
}


# fit: fit object
# data: data to sample (randomly) from for newx
# s: specific lambda value

prediction_fun <- function(fit, data, s){
  test_indices <- sample(nrow(data), 10000)
  last_col <- ncol(data)
  te_x <- as.matrix(data[test_indices, 1:(last_col)])
  return(predict(fit, newx=te_x, s=s, type="class"))
}


# data: dataframe of data to sample from

fitter <- function(data, n){
  train_data <- train_set(data, n)
  tr_x <- train_data[[1]]
  tr_y <- train_data[[2]]
  
  fit <- glmnet(tr_x, tr_y, family="binomial")
  cv_fit <- cv.glmnet(tr_x, tr_y, family="binomial", type.measure="class")
  
  return(list(fit, cv_fit))
}


# Create a confusion matrix from
# x, y: data and labels for running a cross validation. (At time of creation, these are ordered)
# zflm: the zoom-fit lambda min. Previous cross-validation run's lambda.min
# iterations: number of times to run cross-validation
# folds: i.e. the proportion of data to hold out.
confusion_f <- function(x, y, zflm, iterations, folds){
  confusion <- matrix(rep(0, 4), nrow=2, ncol=2)
  
  for (i in 1:iterations){
    cutoff <- floor(nrow(x)/folds)
    
    shuffle_i <- sample(nrow(x))
    while (mean(y[shuffle_i,][1:cutoff]) == 1){
      shuffle_i <- sample(nrow(x))
    }
    
    held_szx <- x[shuffle_i,][1:cutoff,]
    held_szy <- y[shuffle_i,][1:cutoff]
    
    szx <- x[shuffle_i,][cutoff:nrow(x),]
    szy <- y[shuffle_i,][cutoff:nrow(y)]
    
    fit <- glmnet(szx, szy, family="binomial", type.measure="class", lambda=zflm)
    confusion_matrix <- confusion.glmnet(fit, newx=held_szx, newy=held_szy, family="binomial")
    
    if (nrow(confusion_matrix) == 1){
      confusion[1, 1] <- confusion[1, 1] + confusion_matrix[1, 1]
      confusion[1, 2] <- confusion[1, 2] + confusion_matrix[1, 2]
    }
    if (nrow(confusion_matrix) == 2){
      confusion[1, 1] <- confusion[1, 1] + confusion_matrix[1, 1]
      confusion[1, 2] <- confusion[1, 2] + confusion_matrix[1, 2]
      confusion[2, 1] <- confusion[2, 1] + confusion_matrix[2, 1]
      confusion[2, 2] <- confusion[2, 2] + confusion_matrix[2, 2]
    }
  }
  
  cdf <- data.frame(confusion)
  rownames(cdf) <- c("Predict A1", "Predict A2")
  colnames(cdf) <- c("True A1", "True A2")
  fpr <- (cdf[2,1] / (cdf[1,1] + cdf[2, 1]))
  fnr <- (cdf[1,2] / (cdf[2,2] + cdf[1, 2]))

  return(c(fnr, fpr))
}


# Extracts the non-zero coefficients from a cross-validated model and a lambda value
extract_nzc <- function(fit, lm){
  coefs <- coef(fit, s=lm)
  nzc <- coefs@Dimnames[[1]][coefs@i + 1]
  nzc <- nzc[2:length(nzc)]
  return(nzc)
}


# layers: How many iterations to run the model for. Minimum 1.
# min_params: downstream layers will perform poorly with too few parameters, adjust as needed
layered_model <- function(x, y, xy, layers, alpha1=1, alphai=1, min_params=15){
  # Store false rates
  fnrs <- rep(0, layers)
  fprs <- rep(0, layers)
  oErrs <- rep(0, layers)
  
  # Apparently the optimal alpha is 0.57
  # Generate pseudo-stratified fold IDs
  foldids <- sample(1:5, nrow(xy), replace=T)
  while (abs(mean(foldids[1:19]) - 3) > 0.05){
    foldids <- sample(1:5, nrow(xy), replace=T)
  }
  
  # ----- LAYER 1 -----
  l1_fit <- cv.glmnet(x, y, nfolds=5, foldid=foldids, family="gaussian", alpha=alpha1)
  l1_lm <- l1_fit$lambda.min
  l1_nzc <- extract_nzc(l1_fit, l1_lm)
  
  # Do not use models that have fit too few parameters
  while (length(l1_nzc) < min_params){
    l1_fit <- cv.glmnet(x, y, nfolds=5, family="gaussian", alpha=alpha1)
    l1_lm <- l1_fit$lambda.min
    l1_nzc <- extract_nzc(l1_fit, l1_lm)
  }
  #print(paste("Layer 1 MSE: ", min(l1_fit$cvm)))
  # plot(l1_fit)
  # print(paste(length(l1_nzc), "parameters left after first layer."))
  
  frs <- confusion_f(x, y, l1_lm, 1000, 5)
  oErrs[[1]] <- min(l1_fit$cvm)
  fprs[[1]] <- frs[[1]]
  fnrs[[1]] <- frs[[2]]
  
  # ----- LAYERS 2+ ------
  prev_nzc <- l1_nzc
  li_lm <- l1_lm
  
  optimized <- F
  for (i in seq(2:(layers))){
    li_xy <- xy[append(prev_nzc, c("y"))]
    li_x <- as.matrix(li_xy %>% select(-y))
    li_y <- as.matrix(li_xy %>% select(y))
    
    li_fit <- cv.glmnet(li_x, li_y, foldid=foldids, nfolds=5, family="binomial", type.measure="class", alpha=alphai)
    
    # If lambda.min is no longer changing, we know it's optimized.
    if (li_fit$lambda.min == li_lm & optimized == F){
      #print(paste("Optimized at layer", i))
      optimized <- T
    }
    
    li_lm <- li_fit$lambda.min
    li_nzc <- extract_nzc(li_fit, li_lm)
    prev_nzc <- li_nzc
    
    #print(paste("Layer", i+1, "MCE:", min(li_fit$cvm)))
    # plot(li_fit)
    
    frs <- confusion_f(li_x, li_y, li_lm, 1000, 5)
    oErrs[[i+1]] <- min(li_fit$cvm)
    fprs[[i+1]] <- frs[[1]]
    fnrs[[i+1]] <- frs[[2]]
  }
  
  #print(fnrs)
  #print(fprs)
  return(list(prev_nzc, oErrs, fnrs, fprs))
}

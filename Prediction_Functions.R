VecsToKernel <- function(x.train, x.test, FUN = GaussianKernel, bw){
## This function takes as input a vector x 
## and a function that is a type kernel
## Inputs:
##   - x.train: a vector of size N
##   - x.test: a vector of size M
##   - FUN: a kernel function 
## Output: a matrix of size M -by- N with the values 
##         K(x[i], x[j]) in position (i,j)
   n.train <- length(x.train)
   n.test  <- length(x.test)
   kern.mat <- matrix(0, nrow = n.test, ncol = n.train)
   test.mat <- matrix(rep(x.test, n.train), nrow = n.test, ncol = n.train)
   train.mat <- matrix(rep(x.train, n.test), nrow = n.test, ncol = n.train, 
                       byrow = TRUE)
   kern.mat <- FUN(test.mat - train.mat, bw = bw)
   sum.kern.mat <- matrix(rep(rowSums(kern.mat), n.train), 
                          nrow = n.test, ncol = n.train, 
                          byrow = TRUE)
   return((kern.mat / sum.kern.mat))
}

MatricesToKernel <- function(X.train, X.test, FUN = GaussianKernel){
## This function takes  matrices of observations X.train and X.test
## and computes kernels by applying the function 
## VecsToKernel to each column of the matrices 
## Inputs:
##   - X.train: a matrix of dimension n.train -by- num.resp (K)
##        each column is the same type of predictors for all the responses
##   - FUN: a kernel function
## Output: the function returns a matrix of dimension 
##         num.obs - by- num_of_obs * num.resp (K)
##         where we have [kern.mat_1 |... |kern.mat_K]
  if(dim(X.train)[2] != dim(X.test)[2])
     stop("check the dimension of the matrices X.train and X.test");
  n.train <- dim(X.train)[1]
  n.test  <- dim(X.test)[1]
  n.resp  <- dim(X.train)[2]
 
  kernel.matrix <- matrix(0, nrow = n.test, ncol = (n.train * n.resp))
  beg.idx <- 1
  end.idx <- n.train
  for (resp.idx in seq_len(n.resp)){
     bth <- 0.6 * sd(X.train[,resp.idx]) * (n.train)^(-1/5)
     kernel.matrix[,beg.idx:end.idx] <- VecsToKernel(X.train[, resp.idx], X.test[, resp.idx], 
                                                     FUN, bth)
     beg.idx <- end.idx + 1
     end.idx <- end.idx + num.obs
  }
  return(kernel.matrix)
}

ArraysToKernel <- function(X.train, X.test, FUN = GaussianKernel){
# returns an array of kernel of dimension n.test -by- n.train * n.resp
# -by- n.pred
   dim.train <- dim(X.train)
   dim.test  <- dim(X.test)
   if ((dim.train[2] != dim.test[2]) || (dim.train[3] != dim.test[3]))
      stop("Check the dimensions of arrays X.train and X.test")
   n.train <- dim.train[1]
   n.test  <- dim.test[1]
   n.pred <- dim.train[2] 
   n.resp <- dim.train[3]
   kernel.array <- array(0, dim = c(n.test, n.train * n.resp, n.pred))
   for (pred.idx in seq_len(n.pred)){
      kernel.array[, , pred.idx] <- MatricesToKernel(as.matrix(X.train[,pred.idx,]), 
                                    as.matrix(X.test[,pred.idx,]), FUN)
   }
   return(kernel.array)
}

PredictAdditiveFunctions <- function(X.train, X.test, additive.models, pred.lam.idx,
                                     FUN = GaussianKernel){
## This function takes in the training data, the test data, 
## the fitted additive models for the multitask model
## Inputs:
##   - X.train: an array of dimension n.train -by- n.pred -by- n.resp 
##   - X.test: an array of dimension n.test -by- n.pred -by- n.resp 
##   - additive.models: output of the function MTSpAM_C with one lambda
##   - pred.lam.idx: lambda at which forecasting is done
## Output: returns an array of dimension n.test -by- n.pred -by- n.resp
##         with the fit values of the inputs in X.test 

  if ((dim(X.train)[2] != dim(X.test)[2]) || (dim(X.train)[3] != dim(X.test)[3]))
      stop("Check the dimensions of arrays X.train and X.test")  
  # Compute the smoothing kernel arrays for the test inputs
  # smoothers array[n.test, n.train* n.resp, n.pred]

  # add_fun array[n.train, n.pred, n.resp]
  add_fun  <- additive.models[[1]][[pred.lam.idx]] 
  n.pred <- dim(X.train)[2]
  n.test <- dim(X.test)[1]
  n.resp <- dim(X.train)[3]
  n.train <- dim(X.train)[1]
  add_fun_norm <- apply(add_fun, MARGIN = c(2, 3), FUN = Norm2)
  prediction.arr <- array(0, dim = c(n.test, n.pred, n.resp))
  for (pred.idx in seq_len(n.pred)){
    for (resp.idx in seq_len(n.resp)){
        if (add_fun_norm[pred.idx, resp.idx] != 0){
          covariate <- X.train[,pred.idx, resp.idx]
          add.fun <- add_fun[, pred.idx, resp.idx]
          dat.fra <- data.frame(cbind(covariate, add.fun))
          model   <- loess(add.fun ~ covariate, data = dat.fra, 
                          control=loess.control(surface="direct"))
          new.data <- data.frame(covariate = X.test[,pred.idx, resp.idx])
          prediction.arr[, pred.idx,resp.idx] <- predict(model, 
                                               newdata = new.data)
        }
    }
  }
  return(prediction.arr)
}


cross.validation <- function(cov.train, responses, FUN = GaussianKernel, max.it = 1000,
                             tol = 10^(-5), lam.scale = (1 - seq(0, 0.95, by = 0.025)),
                             k.fold = 10, max.lam = 4, max.lam.it = 50, max.lam.tol = 10^(-3), 
                             sm.type = sm.type, deg = deg){
## This function returns the cross validated Mean Squared Error 
## Inputs: 
##   - cov.train: array of dimension n.obs -by- n.pred -by- n.resp with the 
##                covariates of the multitask model
##   - responses: matrix of dimension n.obs -by- n.resp 
##   - max.it: maximum number of iterations for the model 
##   
   n.obs <- dim(cov.train)[1]
   n.pred <- dim(cov.train)[2]
   n.resp <- dim(cov.train)[3]
   test.seq <- matrix(sample(seq_len(n.obs)), nrow = k.fold)
   n.test <- dim(test.seq)[2]
   n.val <- n.obs - n.test
   cov.cv  <- array(0, dim = c(n.val, n.pred ,n.resp)) 
   resp.cv <- matrix(0, n.val, n.resp)
   cov.test  <- array(0, dim = c(n.test, n.pred ,n.resp)) 
   resp.test <- matrix(0, n.test, n.resp)
   # Compute the MaxLambda Value 
   if (sm.type == 1){
     smoothers <- ArrayToKernel(cov.train, FUN) 
   } else {
     smoothers <- ArrayToSmoothing(cov.train, deg = deg)
   }
   max.lam <- MTLambdaMaxR(responses, cov.train, smoothers, max.lam, 
                           max.it, tol, max.lam.it, max.lam.tol)
   lam.seq <- max.lam * lam.scale
   n.lam <- length(lam.seq)
   mse.lam <- rep(0, n.lam)
   for (k.idx in seq_len(k.fold)){
      cov.cv[,,] <- cov.train[-test.seq[k.idx,], , ]
      resp.cv[,] <- responses[-test.seq[k.idx,], ]
      additive.models <- MTSpAM_C(resp.mat = resp.cv, pred.arr = cov.cv,
                                  FUN = FUN, lambda.seq = lam.seq, max.it = max.it, tol = tol, 
                                   sm.type = sm.type, deg = deg)
      for (lam.idx in seq_len(n.lam)){
         cov.test[,,] <- cov.train[test.seq[k.idx, ], , ]
         resp.test[,] <- responses[test.seq[k.idx,], ]
         pred.fun.arr <- PredictAdditiveFunctions(cov.cv, cov.test, additive.models, lam.idx,
                                                  FUN)
         predicted.responses <- apply(pred.fun.arr, MARGIN = c(1, 3), FUN = sum)
         err.mat <- resp.test - predicted.responses
         mse <- mean(colMeans(err.mat^2))
         mse.lam[lam.idx] <- mse.lam[lam.idx] + mse
      }
   }
   mse.lam <- mse.lam / k.fold
   return(list(mse = mse.lam, lambda.seq = lam.seq))
}


cross.validation.smt <- function(cov.train, responses, FUN = GaussianKernel, max.it = 1000,
                                 tol = 10^(-5), alpha.seq = seq(0, 1, by = 0.1),
                                 lam.scale = (1 - seq(0, 0.95, by = 0.05)), lam.seq = NULL,
                                 k.fold = 10, max.lam = 4, max.lam.it = 50, max.lam.tol = 10^(-3), 
                                 sm.type = 1, deg = 1, bw = 0.8){
## This function returns the cross validated Mean Squared Error 
## Inputs: 
##   - cov.train: array of dimension n.obs -by- n.pred -by- n.resp with the 
##                covariates of the multitask model
##   - responses: matrix of dimension n.obs -by- n.resp 
##   - max.it: maximum number of iterations for the model 
##   
   n.obs <- dim(cov.train)[1]
   n.pred <- dim(cov.train)[2]
   n.resp <- dim(cov.train)[3]
   #test.seq <- matrix(sample(seq_len(n.obs)), nrow = k.fold)
   test.seq <- matrix(seq_len(n.obs), nrow = k.fold)
   n.test <- dim(test.seq)[2]
   n.val <- n.obs - n.test
   cov.cv  <- array(0, dim = c(n.val, n.pred ,n.resp)) 
   resp.cv <- matrix(0, n.val, n.resp)
   cov.test  <- array(0, dim = c(n.test, n.pred ,n.resp)) 
   resp.test <- matrix(0, n.test, n.resp)
   # Compute the MaxLambda Value 
   if (sm.type == 1){
     smoothers <- ArrayToKernel(cov.train, FUN, bw = bw) 
   } else {
     smoothers <- ArrayToSmoothing(cov.train, deg = deg, bw = bw)
   }
   if (is.null(lam.seq)){
     max.lam <- MTLambdaMaxR(responses, cov.train, smoothers, max.lam, 
                             max.it, tol, max.lam.it, max.lam.tol)
     lam.seq <- max.lam * lam.scale
   }
   n.alpha <- length(alpha.seq)
   n.lam <- length(lam.seq)
   mse.lam <- matrix(rep(0, n.lam), nrow = n.lam,
                     ncol = n.alpha)
   for (k.idx in seq_len(k.fold)){
      cov.cv[,,] <- cov.train[-test.seq[k.idx,], , ]
      resp.cv[,] <- responses[-test.seq[k.idx,], ]
      for (alpha.idx in seq_len(n.alpha)){
        alpha <- alpha.seq[alpha.idx]
        additive.models <- SMTSpAM_C(resp.mat = resp.cv, pred.arr = cov.cv, alpha = alpha,
                                     lambda.seq = lam.seq, FUN = FUN,
                                     max.it = max.it, tol = tol, sm.type = sm.type, deg = deg, bw = bw)
        for (lam.idx in seq_len(n.lam)){
           cov.test[,,] <- cov.train[test.seq[k.idx, ], , ]
           resp.test[,] <- responses[test.seq[k.idx,], ]
           pred.fun.arr <- PredictAdditiveFunctions(cov.cv, cov.test, additive.models, lam.idx,
                                                  FUN)
           predicted.responses <- apply(pred.fun.arr, MARGIN = c(1, 3), FUN = sum)
           err.mat <- resp.test - predicted.responses
           mse <- median(colMeans(err.mat^2))
           mse.lam[lam.idx, alpha.idx] <- mse.lam[lam.idx, alpha.idx] + mse
        }
      }
   }
   mse.lam <- mse.lam / k.fold
   return(list(mse = mse.lam, lambda.seq = lam.seq))
}



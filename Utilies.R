Norm2 <- function(x){
## This function is a vector of length num.obs
## returns the norm 2 of a vector 
  return(sqrt(sum(x^2)))
}

TriCubeWeights <- function(u, bw = 0.05){
# This function computes the tricube function 
# Output: return the tricube value
  EX <- ifelse(abs(u/bw) < 1, (1 - abs(u/bw)^3)^3, 0)
  return (EX)
 #return((1 - abs(u/bw)^2)^3 * as.numeric(abs(u/bw) < 1))
}

VecToSmoothing <- function(x, bw = 0.8, deg = 1){
  num.elts <- length(x) 
  weight.mat <- matrix(0, nrow = num.elts, ncol = num.elts)
  for (idx in seq_len(num.elts)){
       des.mat <- cbind(rep(1, num.elts), (x - x[idx]))
       h_i <- median(x - x[idx])
       tw <- TriCubeWeights(x - x[idx], bw = bw)
       S_1 <- sum(tw * (x - x[idx]))
       S_2 <- sum(tw * (x - x[idx])^2)
       Bi  <- tw *(S_2 - (x - x[idx])*S_1)
       if (sum(Bi) != 0) {
         weight.mat[idx,] <- Bi / sum(Bi) 
       } else {
         weight.mat[idx,] <- Bi
       }
  }
 
  return(weight.mat)
}

VecToSmoothing_ <- function(x, bw = 1, deg = 1){

  num.elts <- length(x) 
  weight.mat <- matrix(0, nrow = num.elts, ncol = num.elts)
  if (deg != 1){
    for (idx in seq_len(num.elts)){
       des.mat <- cbind(rep(1, num.elts), (x - x[idx]), (x - x[idx])^2)
       temp.W <- TriCubeWeights(x - x[idx], bw = bw) / sum(TriCubeWeights(x - x[idx], bw = bw))
       #W <- diag(TriCubeWeights(x - x[idx], bw = bw))
        W <- diag(temp.W)
       weight.mat[idx,]  <- c(1, 0, 0)%*% solve(t(des.mat) %*% W %*% des.mat)%*% 
                          t(des.mat) %*%  W  
    } 
  } else {
  for (idx in seq_len(num.elts)){
       des.mat <- cbind(rep(1, num.elts), (x - x[idx]))
       tw <- TriCubeWeights(x - x[idx], bw = bw)
       #if (sum(tw) == 0) {
       #    W <- diag(tw)
       #} else {
       #W <- diag(tw / sum(tw)) }
       W <- diag(tw)
       tp <- t(des.mat) %*% diag(tw) %*% des.mat
       tp.inv <- matrix(c(tp[2,2],-tp[2, 1], - tp[1,2], tp[1,1]), 2, 2) /
                 (tp[2,2] * tp[1, 1] - tp[2, 1] * tp[1,2])
      # weight.mat[idx,]  <- c(1, 0)%*% solve(t(des.mat) %*% W %*% des.mat)%*% 
      #                    t(des.mat) %*%  W  
       weight.mat[idx,]  <- c(1, 0)%*% tp.inv%*% 
                          t(des.mat) %*%  W  
       weight.mat[idx,] <- weight.mat[idx,] / sum(weight.mat[idx,])
    } 
  }
  
  return(weight.mat)
}

MatrixToSmoothing <- function(X, deg = 1, bw = NULL){
## This function takes a matrix of observations X
## and computes kernels by applying the function 
## VecToKernel to each column of the matrix
## Inputs:
##   - X: a matrix of dimension num.obs -by- num.resp (K)
##        each column is the same type of predictors for all the responses
##   - FUN: a kernel function
## Output: the function returns a matrix of dimension 
##         num.obs - by- num_of_obs * num.resp (K)
##         where we have [kern.mat_1 |... |kern.mat_K]
 
   num.obs  <- nrow(X)
   num.resp <- ncol(X)
   kernel.matrix <- matrix(0, nrow = num.obs, ncol = (num.obs * num.resp))
   beg.idx <- 1
   end.idx <- num.obs
   for (resp.idx in seq_len(num.resp)){
      if (is.null(bw)){
        bth <- 0.8 * sd(X[,resp.idx]) * (num.obs)^(-1/5) 
      } else {
        bth <- bw * sd(X[,resp.idx]) * (num.obs)^(-1/5) 
      }
      kernel.matrix[ ,beg.idx:end.idx] <- VecToSmoothing(X[,resp.idx], bth, deg)
      beg.idx <- end.idx + 1
      end.idx <- end.idx + num.obs  
   }
   return(kernel.matrix)
}

ArrayToSmoothing <- function(X, deg = 1, bw = NULL){
## This function takes a 3d-array of observations X
## and computes kernels by applying the function 
## MatrixToKernel to each matrix of the 3d-array X[,j,]
## Inputs:
##   - X: an array of dimension num.obs -by- num.pred (P)
##        -by- num_of_resp (K)
##        each matrix X[,,k] is the set of predictors for a 
##        regression problem
##   - FUN: a kernel function
## Output: the function returns a 3d array of dimension 
##         num.obs - by- num.obs * num.resp (K) -by- 
##         num.pred
##         where we have [kern.mat_1 |... |kern.mat_K]_j
##         with j is an index of the predictors

   dim.cov <- dim(X)
   num.obs <- dim.cov[1]
   num.pred <- dim.cov[2]
   num.resp <- dim.cov[3] 
   kernel.array <- array(0, dim = c(num.obs, 
                   num.obs * num.resp, num.pred))
   for (pred.idx in seq_len(num.pred)){
      kernel.array[ , ,pred.idx] <- MatrixToSmoothing(as.matrix(X[ ,pred.idx, ]), deg = deg, 
                                    bw = bw)
   } 
   return(kernel.array)
}

GaussianKernel <- function(u, bw = 0.05){
## This function computes the Gaussian Kernel
## Input:
##   - u: a number used to compute the Kernel function 
## Output: the value of the gaussian kernel
  gauss.kern <- (1/sqrt(2 *pi)) * exp(-0.5 *(u/bw)^2)
  return(gauss.kern)
}


VecToKernel <- function(x, FUN = GaussianKernel, bw){
## This function takes as input a vector x 
## and a function that is a type kernel
## Inputs:
##   - x: a vector of size N
##   - FUN: a kernel function 
## Output: a matrix of size N -by- N with the values 
##         K(x[i], x[j]) in position (i,j)
   num.elts <- length(x)
   kern.mat <- matrix(0, nrow = num.elts, ncol = num.elts)
   row.idx <- row(kern.mat)
   col.idx <- col(kern.mat)
   row.elts <- x[row.idx]
   col.elts <- x[col.idx] 
   diff.all.x <- FUN(row.elts - col.elts, bw = bw)
   kern.mat <- matrix(diff.all.x, num.elts, num.elts)
   sum.kern.mat <- matrix(rep(rowSums(kern.mat), num.elts), 
                          nrow = num.elts, ncol = num.elts, 
                          byrow = TRUE)
   return((kern.mat / sum.kern.mat))
}

MatrixToKernel <- function(X, FUN = GaussianKernel, bw = .8){
## This function takes a matrix of observations X
## and computes kernels by applying the function 
## VecToKernel to each column of the matrix
## Inputs:
##   - X: a matrix of dimension num.obs -by- num.resp (K)
##        each column is the same type of predictors for all the responses
##   - FUN: a kernel function
## Output: the function returns a matrix of dimension 
##         num.obs - by- num_of_obs * num.resp (K)
##         where we have [kern.mat_1 |... |kern.mat_K]
 
   num.obs  <- nrow(X)
   num.resp <- ncol(X)
   kernel.matrix <- matrix(0, nrow = num.obs, ncol = (num.obs * num.resp))
   beg.idx <- 1
   end.idx <- num.obs
   for (resp.idx in seq_len(num.resp)){
      bth <- bw * sd(X[,resp.idx]) * (num.obs)^(-1/5)
      kernel.matrix[ ,beg.idx:end.idx] <- VecToKernel(X[,resp.idx], FUN, bth)
      beg.idx <- end.idx + 1
      end.idx <- end.idx + num.obs  
   }
   return(kernel.matrix)
}

ArrayToKernel <- function(X, FUN = GaussianKernel, bw = .8){
## This function takes a 3d-array of observations X
## and computes kernels by applying the function 
## MatrixToKernel to each matrix of the 3d-array X[,j,]
## Inputs:
##   - X: an array of dimension num.obs -by- num.pred (P)
##        -by- num_of_resp (K)
##        each matrix X[,,k] is the set of predictors for a 
##        regression problem
##   - FUN: a kernel function
## Output: the function returns a 3d array of dimension 
##         num.obs - by- num.obs * num.resp (K) -by- 
##         num.pred
##         where we have [kern.mat_1 |... |kern.mat_K]_j
##         with j is an index of the predictors

   dim.cov <- dim(X)
   num.obs <- dim.cov[1]
   num.pred <- dim.cov[2]
   num.resp <- dim.cov[3] 
   kernel.array <- array(0, dim = c(num.obs, 
                   num.obs * num.resp, num.pred))
   for (pred.idx in seq_len(num.pred)){
      kernel.array[ , ,pred.idx] <- MatrixToKernel(as.matrix(X[ ,pred.idx, ]), FUN, 
                                    bw = bw)
   } 
   return(kernel.array)
}

UpdateMTAdditiveComponent <- function(kernel.mat, part.res.mat, lambda){
## This function takes as input a matrix with kernel entries
## of dimension num.obs -by- num.obs * num.pred and 
## part.res.mat a partial residual matrix of dimension 
## num.obs -by- num.resp(K) and returns a matrix where each 
## column is needed to compute the updated additive functions 
## (f_j^[1], ..., f_j^[K])
## Inputs:
##   - kernel.mat: a matrix of dimension num.obs -by- num.obs * num.resp
##                 all the entries are related to same type of predictors
##                 kernel.mat <- kernel.array[ , ,j]
##   - part.res.mat: a matrix of dimension num.obs -by- num.resp (K)
##                   that contains the partial residuals for all the responses
## Output: returns the matrix of updated additive functions
##         [f_j(1)|...| f_j(K)] of dimension num.obs -by- num.resp
##         also returns the l_2\l_1 functional norm
   
   num.obs  <- dim(kernel.mat)[1]
   num.resp <- dim(kernel.mat)[2] / (num.obs)
   add.comp.mat <- matrix(0, num.obs, num.resp)
   beg.idx <- 1
   end.idx <- num.obs
   for (resp.idx in seq_len(num.resp)){
      add.comp.mat[,resp.idx] <- kernel.mat[ ,beg.idx:end.idx] %*% 
                                 part.res.mat[ ,resp.idx]   
      beg.idx <- end.idx + 1
      end.idx <- end.idx + num.obs 
   }   
   Wj <- sqrt(sum((apply(add.comp.mat, MARGIN = 2, FUN = Norm2)^2) /
                    num.obs))
   add.comp.mat <- max(1 - (lambda * sqrt(num.resp))/Wj, 0) * 
                   add.comp.mat   
   add.mean <- rep(colMeans(add.comp.mat), nrow = num.obs, ncol = num.resp,
                   byrow = TRUE)
   add.comp.mat <- add.comp.mat - add.mean
   add.comp.norms <- apply(add.comp.mat, MARGIN = 2, FUN = Norm2)
   f.group.norm <- sqrt(sum(add.comp.norms^2))
   return(list(fj = add.comp.mat, f.group.norm = f.group.norm,
               f.norms = add.comp.norms))   
}

UpdateMTPartialResiduals <- function(part.res, add.fun.arr, pred.idx){
## This function is used to update the partial residual 
## Inputs: 
##   - part.res: matrix of partial residuals of dimension 
##   - add.fun.arr: 3D-array of dimension num.obs -by- num.pred -by- num.resp
##                  containing all the additive functions  
##   - pred.idx: index of the predictor of interest
## Output: returns an updated partial residual matrix of dimension 
##         num_of_obs -by- num.resp
##    
   add.mat <- add.fun.arr[ ,pred.idx, ]
   part.res <- part.res - add.mat
   return(part.res)
}

ComputePartialResiduals <- function(resp.mat, add.fun.arr, pred.idx){
## This function returns the partial residuals when the predictor 
## at pred.idx is not accounted for
## Inputs: 
##   - resp.mat: matrices of the responses of dimension num.pred -by- 
##               num.resp
##   - add.fun.arr: 3D-array of dimension num.obs -by- num.pred -by- num.resp
##                  containing all the additive functions  
##   - pred.idx: index of the predictor of interest
## Output: returns the partial residual for all the responses Y^(k)
##         when the predictors [x_j^(1), ..., x_j^(K)] are left out 
##         with j being the predictor of interest
   num.obs <- dim(resp.mat)[1]
   num.resp <- dim(resp.mat)[2]
   num.pred <- dim(add.fun.arr)[3]
   part.res <- resp.mat - apply(add.fun.arr[ ,-pred.idx, ], MARGIN = c(1, 3), 
                                FUN = sum)
   return(part.res)
}

BackfittingMTSpAM <- function(resp.mat, pred.arr, init.fun.arr, 
                            lambda, FUN = GaussianKernel,
                            max.it = 100, tol = 10^(-3), sm.type = 1, deg = 1, bw =.8){
## This function returns the estimated additive functions in the 
## the multi-task (or multivariate) regression setting
## Input:
##   - resp.mat: matrix of dimension num.obs -by- num.resp that 
##               contains all the responses
##   - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
##   - init.fun.arr: array of dimension num.obs - by- num.pred -by- num.resp
##                   that contains the additive function values
##   - lambda: regularization parameter value
##   - max.it: maximum number of iteration for the outerloop 
##   - tol: tolerance level for the convergence criterion
## Outputs: 
##   - fun.arr: returns an array of dimension  num.obs -by- num.pred -by-
##              num.resp with the estimated additive functions 
##   - fun.norm: the l2 \ l1 norm of the functions f_j
  if((dim(resp.mat)[1] != dim(pred.arr)[1])|| 
     (dim(resp.mat)[1] != dim(init.fun.arr)[1])){
    stop("check dimensions (rows) of resp.mat and pred.arr and init.fun.arr")
  }
  if((dim(resp.mat)[2] != dim(pred.arr)[3]) || 
     (dim(resp.mat)[2] != dim(init.fun.arr)[3])){
    stop("check dimensions (cols) of resp.mat and  3rd dim of pred.arr and init.fun.arr")
  }
  
  num.obs <- dim(resp.mat)[1]
  num.pred <- dim(pred.arr)[2]
  num.resp <- dim(resp.mat)[2]
  # initialize the array of additive functions
  fun.arr <- init.fun.arr
  fun.norm <- matrix(0, num.resp, num.pred)
  # Compute the matrix of Smoothers labelled as S in the paper
  if (sm.type == 1){
    Smoothers <- ArrayToKernel(pred.arr, FUN, bw = bw)
  } else {
    Smoothers <- ArrayToSmoothing(pred.arr, deg = deg, bw = bw)
  }
  it <- 1
  conv.crit <- 1
  part.res.mat <- resp.mat - matrix(rep(colMeans(resp.mat), num.obs), nrow = num.obs, 
                                    ncol = num.resp, byrow = TRUE)
  while ((it < max.it) & (conv.crit > tol)){
    old.fun.arr <- fun.arr
    old.fun.norm <- fun.norm
    for (pred.idx in seq_len(num.pred)){
      # Compute the partial residuals 
      # part.res.mat labelled as R in the paper 
      part.res.mat <- UpdateMTPartialResiduals(part.res.mat, -fun.arr, pred.idx)     
      #part.res.mat <- ComputePartialResiduals(resp.mat, fun.arr, pred.idx)     

      # Use the Algo1 of the paper to update the additive functions 
      fun.arr.norm <- UpdateMTAdditiveComponent(Smoothers[,,pred.idx], part.res.mat, lambda)
      fun.arr[,pred.idx,] <- fun.arr.norm[[1]]
      fun.norm[ ,pred.idx] <- fun.arr.norm[[3]]
      part.res.mat <- UpdateMTPartialResiduals(part.res.mat, fun.arr, pred.idx)     
    }
    conv.crit <- sqrt(sum(sum((old.fun.norm - fun.norm)^2)))
  }
  fun.group.norm <- fun.norm 
  fun.norms <- fun.arr.norm[[3]]
  return(list(fun.arr = fun.arr, fun.norm = fun.norm, part.res.mat = part.res.mat))
}

ComputeGCV <- function(resp.mat, fun.arr, fun.norms, smoothers){
## This function computes the Generalized Cross Validation for the Multitask 
## L2\L1 additive model
## Inputs:
##   - resp.mat: matrix of dimension num.obs -by- num.resp 
##   - fun.arr: an array of dimension num.obs -by- num.pred -by- num.resp
##   - fun.norms: a matrix of dimension num.resp -by- num.pred
##   - smoothers: an array of dimension num.obs -by- (num.obs * num.resp)
##                                              -by- num.resp
##  Outputs: The function returns the GCV for the model 
   num.obs <- dim(resp.mat)[1]
   num.resp <- dim(resp.mat)[2]
   num.pred <- dim(fun.arr)[2]
   residuals <- (resp.mat - apply(fun.arr, MARGIN = c(1, 3), FUN = sum))
   least.square <- sum(sum(residuals^2)) 
   smooth.traces <- NULL
   ind.fun       <- NULL
   for (pred.idx in seq_len(num.pred)){
       mat.idx <- 1
       for (resp.idx in seq_len(num.resp)){
          cur.trace <- sum(diag(smoothers[ ,mat.idx:(mat.idx + num.obs - 1), 
                                           pred.idx]))
          smooth.traces <- c(smooth.traces, cur.trace)
          ind.fun <- c(ind.fun, as.numeric(fun.norms[resp.idx, pred.idx] != 0))
          mat.idx <- mat.idx + num.obs
       }
   }
   df <- sum(smooth.traces * ind.fun)
   gcv <- least.square / (num.obs^2*num.resp^2 - (num.obs * num.resp)*df)^2
   return(gcv)
}


MTLambdaMax <- function(resp.mat, pred.arr, FUN = GaussianKernel,
                        max.lam = 4, max.it = 100, tol = 10^(-3), 
                        max.lam.it = 10, max.lam.tol = 10^(-2), sm.type = 1, deg = 1, bw = .8){
## This function returns the lambda value at which the first value 
## enter in the model 
  upper.lam <- max.lam 
  lower.lam <- 0
  est.lam <- upper.lam 
  gap <- upper.lam  - lower.lam
  dim.arr <- dim(pred.arr)
  init.fun.arr <- array(0, dim = c(dim.arr[1], dim.arr[2], dim.arr[3]))
  not.found <- TRUE
  it <- 1
  while ((gap > max.lam.tol) && (it < max.lam.it)){
    model <- BackfittingMTSpAM(resp.mat, pred.arr, init.fun.arr, 
                            lambda = est.lam, FUN = GaussianKernel,
                            max.it = max.it, tol = tol, sm.type = sm.type, deg = deg, bw = bw)
    max.fun.norm <- max(colMeans(model[[2]]))
    if (max.fun.norm == 0){
      upper.lam <- est.lam
      est.lam <- upper.lam - (upper.lam - lower.lam) / 2
    } else {
       lower.lam <- est.lam 
       est.lam <- lower.lam + (upper.lam - lower.lam)/2
    }
    gap <- upper.lam - lower.lam
    it <- it + 1
  }
  return(list(upper.lam, max.fun.norm))
}

MTSpAM <- function(resp.mat, pred.arr, lambda.seq = NULL, 
                   lam.scale = (1 - seq(0, 0.95, by = 0.05)), 
                   FUN = GaussianKernel, 
                  max.it = 1000, tol = 10^(-4), sm.type = 1, deg = 1, bw = .8){
## This function is used to fit the additive functions for the whole set 
## of values of the parameter lambda
##   - resp.mat: matrix of dimension num.obs -by- num.resp that 
##               contains all the responses
##   - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
##   - lambda.seq: regularization values vector default is null
##   - max.it: maximum number of iteration for the outerloop 
##   - tol: tolerance level for the convergence criterion
## Outputs:
##   - fun.arr.list: list of size num.lam that contains the arrays of 
##                   dimension num.obs -by- num.pred -by- num.resp
##                   with values of additive functions for each value of 
##                   lambda.
##   - fun.norm.mat: matrix of size num.pred -by- num.lam that contains
##                   the l2/l1 norm of functions f_j for each value of the
##                   the regularization parameter lambda

  fun.arr.list <- list()
  part.res.mat.list <- list()
  if (is.null(lambda.seq)){
   max.lam <- MTLambdaMax(resp.mat, pred.arr = pred.arr, FUN = FUN, sm.type = sm.type, 
                          deg = deg)
   lambda.seq <- lam.scale * max.lam[[1]]
  }  
  num.lam <- length(lambda.seq)
  num.obs <- dim(resp.mat)[1]
  num.pred <- dim(pred.arr)[2]
  num.resp <- dim(resp.mat)[2]
  fun.norm.mat <- NULL
  fun.norm.arr <- array(0, dim= c(num.resp, num.pred, num.lam))
  init.fun.arr <- array(0, dim= c(num.obs, num.pred, num.resp))
  if (sm.type == 1){
    smoothers <- ArrayToKernel(pred.arr, FUN, bw = bw)
  } else {
    smoothers <- ArrayToSmoothing(pred.arr, deg = deg, bw = bw)  
  }
  gcv.vec <- rep(0, num.lam)
  #intercept <- matrix(rep(colMeans(resp.mat), num.obs),
  #                    nrow = num.obs, ncol = num.resp)
  #response.mat <- resp.mat - intercept
  for (lam.idx in seq_len(num.lam)){
     lam <- lambda.seq[lam.idx]
     fun.arr.norm <- BackfittingMTSpAM(resp.mat, pred.arr = pred.arr,
                     init.fun.arr = init.fun.arr, 
                            lam, FUN = FUN, max.it = max.it,
                            tol = tol, sm.type = sm.type, deg = deg, bw = bw)
     init.fun.arr <- fun.arr.norm[[1]]
     fun.arr.list[[lam.idx]] <-  init.fun.arr
     fun.norm.arr[,,lam.idx] <- fun.arr.norm[[2]]
     part.res.mat.list[[lam.idx]] <- fun.arr.norm[[3]]
     # Computation of the GCV fun.arr.norm[[1]] (all the additive functions)
     # fun.arr.norm[[2]] (their norms)
     gcv.vec[lam.idx] <- ComputeGCV(resp.mat, fun.arr.norm[[1]], fun.arr.norm[[2]],
                                    smoothers)
  }
  return(list(fun.arr.list, fun.norm.arr, part.res.mat.list, 
              lambda = lambda.seq, gcv.vec = gcv.vec))
}	



plotMTSpAM <- function(spam.model, num.rows , num.cols, 
                       plot.gcv.curve = FALSE, all.resp = TRUE, spec.resp = 1){
## This function plots the path of the all the additive functions 
## for all the responses
## Inputs:
##   - spam.model: output from method MTSpAM
##   - num.rows: number of rows in the plot frame
##   - num.cols: number of columns in the plot frame
## Output: Void, as a side effect plot the regularization paths 
##         of all the predictors for all the responses
  all.dim <- dim(spam.model[[2]])
  num.pred <- all.dim[2]
  num.resp <- all.dim[1]
  num.lam  <- all.dim[3]
  lambda.seq <- spam.model[[4]]
  max.lambda <- max(lambda.seq)
  lam.x.axis <- (max.lambda - lambda.seq) / max.lambda
  lam.scale <- lam.x.axis
  # determination of the minimum value for the gcv 
  all.gcv <- spam.model[[5]]
  min.lam.idx <- order(all.gcv)[1]
  lam.scal.gcv <- (max.lambda - lambda.seq[min.lam.idx]) / max.lambda
  #[num.lam:1]
  par(mfrow = c(num.rows, num.cols))
  if (all.resp){
    for (resp.idx in seq_len(num.resp)){
       matplot(lam.scale, t(spam.model[[2]][resp.idx,,]), type ="l",
              ylab = "Functional Norms",
       xlab = expression(paste(" Regularization Parameter  ", frac(lambda, lambda[max]))),
       main = paste("Regularization path for Response", resp.idx))
       #abline(v = lam.scal.gcv, lty = 2, col = "red")
       axis(4, at =(spam.model[[2]][resp.idx,,num.lam]), 
            cex.axis = 0.75, label = seq_len(num.pred))
    }
  } else {
       matplot(lam.scale, t(spam.model[[2]][spec.resp,,]), type ="l",
              ylab = "Functional Norms",
       xlab = expression(paste(" Regularization Parameter  ", frac(lambda, lambda[max]))),
       main = paste("Regularization path for Response", spec.resp))
       abline(v = lam.scal.gcv, lty = 2, col = "red")
       axis(4, at =(spam.model[[2]][spec.resp,,num.lam]), 
            cex.axis = 0.75, label = seq_len(num.pred))    
   
  }
  if (plot.gcv.curve){
     X11()
     plot(lambda.seq, all.gcv, xlab = expression(lambda), ylab = "GCV", 
          main = "Generalized Cross Validation")
     abline(v = lambda.scale[min.lam.idx], lty = 2, col = "red")
  }
}

CreateUnivariateInputs <- function(resp.mat, init.fun.arr, 
                                        cov.arr, resp.idx){
## This function is used to create a univariate set of 
## inputs that will be used to fit the model using l2\l1 regularization 
## Inputs:
##   - resp.mat: response matrix of dimension num.obs by num.resp
##   - init.fun.arr: initialization functions array of dimension 
##                   num.obs -by- num.pred -by- num.resp
##   - cov.arr: an array that contains the covariates for 
##              the multitask model 
##              (num.obs - by- num.pred -by- num.resp)
   
     cov.dim <- dim(cov.arr)
     uni.cov.arr <- array(cov.arr[,,resp.idx], 
                          dim =c(cov.dim[1], cov.dim[2], 1))

     fun.dim <- dim(init.fun.arr)
     uni.fun.arr <- array(init.fun.arr[,,resp.idx], 
                     dim =c(fun.dim[1], fun.dim[2], 1))
     uni.resp <- as.matrix(resp.mat[,resp.idx], ncol = 1)

     return(list(resp.mat = uni.resp, init.fun.arr = uni.fun.arr, 
                 cov.arr  = uni.cov.arr))
}


plotFittedAdditiveFun <- function(cov.arr, add.models, responses_addi, weights, plot.lam.idx,
                                  n.plots = 9, p.n.r = 3, p.n.c = 3, 
                                  y.lim = c(-5, 5), resp.plot.id = NULL, max.pred = 4,  
                                  FUN = GenerateTrueFunctions, x.v.min = 0, x.v.max = 1){
## This function plots some of the true additive functions and
## the estimated additive functions
## Inputs:
##   - cov.arr: an array with the covariates of dimension 
##              num.obs -by- num.pred -by- num.resp 
##   - add.models: an object returned by function MTSpAM_C (refer to function)
##   - plot.lam.idx: an integer for the lambda of interest
## Output: void, side effect plot all the additive functions and the 
##         estimated additive functions.
   add_fun <- add.models[[1]][[plot.lam.idx]]
   add.dim <- dim(add_fun)
   n.resp <- add.dim[3]
   n.obs <- add.dim[1]
   n.pred <- add.dim[2]
   x.val <- seq(x.v.min, x.v.max, length = n.obs)
   responses <- responses_addi[[1]]
   s.x.val <- (x.val - min(x.val)) / (max(x.val) - min(x.val))
   for (resp.idx in seq_len(n.resp)){
     if (is.null(resp.plot.id)){
        true.matrix  <- FUN(x.val,  weights, resp.idx,
                                               n.pred)
     } else {
       true.matrix  <- FUN(x.val,  weights, resp.plot.id,
                                           n.pred)
     }
     X11()
     add_par <- responses_addi[[resp.idx + 1]]
     par(mfrow = c(p.n.r, p.n.c))
     for (pred.idx in seq_len(n.plots)){
        cur.cov <- (cov.arr[, pred.idx,resp.idx] - min(cov.arr[, pred.idx,resp.idx])) /
                   (max(cov.arr[, pred.idx,resp.idx])- min(cov.arr[, pred.idx,resp.idx]))
        true.fun <- true.matrix[, pred.idx]
        cur.add.fun <- add_fun[order(cur.cov), pred.idx, resp.idx]
        if(pred.idx <= max.pred){
          true.add <- functional(responses, add_par, idx.fun = pred.idx, 
                                 idx.resp = resp.idx)          
          y.max <- max(c(true.fun, cur.add.fun, true.add))
          y.min <- min(c(true.fun, cur.add.fun, true.add))
          plot(cur.cov, true.add, ylim = c(y.min - 0.5, y.max + 0.5), cex = 0.5)
          lines(s.x.val, true.fun, 
             col = "red", type = "l", lty = 2)
        } else{
          y.max <- max(c(true.fun, cur.add.fun))
          y.min <- min(c(true.fun, cur.add.fun))
          plot(s.x.val, true.fun, 
               col = "red", type = "l", lty = 2,  
               ylim = c(y.min - 0.5, y.max + 0.5))
        } 
        lines(sort(cur.cov), cur.add.fun, type = "l", col = "blue")
      }
   }
}

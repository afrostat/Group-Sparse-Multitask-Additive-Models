source("Utilies.R")
require(Rcpp)
require(RcppArmadillo)
sourceCpp("multitaskspam.cpp")
MTSpAM_C <- function(resp.mat, pred.arr, lambda.seq = NULL, 
                     lam.scale = (1 - seq(0, 0.95, by = 0.05)), 
                     FUN = GaussianKernel, 
                     max.it = 1000, tol = 10^(-4), max.lam = 4,
                     max.lam.it = 50, max.lam.tol = 10^(-6), bth.plugin = NULL){
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
  smoothers <- ArrayToKernel(pred.arr, FUN, bth.plugin = bth.plugin)

  if (is.null(lambda.seq)){
   max.lam <- MTLambdaMaxR(resp.mat, pred.arr, smoothers, max.lam, max.it, tol,
                          max.lam.it, max.lam.tol)
   lambda.seq <- lam.scale * max.lam
  }  
  num.lam <- length(lambda.seq)
  num.obs <- dim(resp.mat)[1]
  num.pred <- dim(pred.arr)[2]
  num.resp <- dim(resp.mat)[2]
  fun.norm.mat <- NULL
  fun.norm.arr <- array(0, dim= c(num.resp, num.pred, num.lam))
  init.fun.arr <- array(0, dim= c(num.obs, num.pred, num.resp))
  gcv.vec <- rep(0, num.lam)
  #intercept <- matrix(rep(colMeans(resp.mat), num.obs),
  #                    nrow = num.obs, ncol = num.resp)
  #response.mat <- resp.mat - intercept
  for (lam.idx in seq_len(num.lam)){
     lam <- lambda.seq[lam.idx]
     fun.arr.norm <- BackfittingMTSpAMR(resp.mat, pred.arr, init.fun.arr, 
                                        smoothers, lam, max.it, tol)
     init.fun.arr <- fun.arr.norm[[1]]
     fun.arr.list[[lam.idx]] <-  init.fun.arr
     fun.norm.arr[,,lam.idx] <- fun.arr.norm[[2]]
     part.res.mat.list[[lam.idx]] <- fun.arr.norm[[3]]
     # Computation of the GCV fun.arr.norm[[1]] (all the additive functions)
     # fun.arr.norm[[2]] (their norms)
     gcv.vec[lam.idx] <- ComputeGCV(resp.mat, fun.arr.norm[[1]], fun.arr.norm[[2]],
                                    smoothers)
  }
  return(list(fun_arr = fun.arr.list, fun_norm = fun.norm.arr, part_res = part.res.mat.list, 
              lambda = lambda.seq, gcv.vec = gcv.vec))
}	


SMTSpAM_C <- function(resp.mat, pred.arr, alpha, lambda.seq = NULL, 
                     lam.scale = (1 - seq(0, 0.95, by = 0.05)), 
                     FUN = GaussianKernel, final.max.lam = NULL, 
                     max.it = 1000, tol = 10^(-4), max.lam = 4,
                     max.lam.it = 50, max.lam.tol = 10^(-6), init.val = 0.01, 
                     bth.plugin = NULL){
## This function is used to fit the additive functions for the whole set 
## of values of the parameters alpha and  lambda
##   - resp.mat: matrix of dimension num.obs -by- num.resp that 
##               contains all the responses
##   - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
##   - alpha.seq: sequence of values for the weigthing parameter alpha
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
  smoothers <- ArrayToKernel(pred.arr, FUN, bth.plugin = bth.plugin)


  # We use the parameter max.lam determined by the method without sparsity 
  # within groups by default to save some time
  if (is.null(lambda.seq) & is.null(final.max.lam)){
   max.lam <- MTLambdaMaxR(resp.mat, pred.arr, smoothers, max.lam, max.it, tol,
                          max.lam.it, max.lam.tol)
   # just to be sure that the max will be 0 for all alphas.
   final.max.lam <- max.lam + 0.5 
   lambda.seq <- lam.scale * final.max.lam
  }
  if (!is.null(final.max.lam) & is.null(lambda.seq)){
      lambda.seq <- lam.scale * final.max.lam
  }  
  num.lam <- length(lambda.seq)
  num.obs <- dim(resp.mat)[1]
  num.pred <- dim(pred.arr)[2]
  num.resp <- dim(resp.mat)[2]
  fun.norm.mat <- NULL
  total <- num.pred * num.obs * num.resp
  fun.norm.arr <- array(0, dim= c(num.resp, num.pred, num.lam))
  # initialize with random values 
  init.fun.arr <- array(runif(total, -init.val, init.val), dim= c(num.obs, num.pred, num.resp))
  gcv.vec <- rep(0, num.lam)
  #intercept <- matrix(rep(colMeans(resp.mat), num.obs),
  #                    nrow = num.obs, ncol = num.resp)
  #response.mat <- resp.mat - intercept
  for (lam.idx in seq_len(num.lam)){
     lam <- lambda.seq[lam.idx]
     fun.arr.norm <- BackfittingSMTSpAMR(resp.mat, pred.arr, init.fun.arr, 
                                         smoothers, alpha, lam, max.it, tol)
     #init.fun.arr <- fun.arr.norm[[1]]
     fun.arr.list[[lam.idx]] <-  fun.arr.norm[[1]]
     fun.norm.arr[,,lam.idx] <-  fun.arr.norm[[2]]
     part.res.mat.list[[lam.idx]] <- fun.arr.norm[[3]]
     # Computation of the GCV fun.arr.norm[[1]] (all the additive functions)
     # fun.arr.norm[[2]] (their norms)
     gcv.vec[lam.idx] <- ComputeGCV(resp.mat, fun.arr.norm[[1]], fun.arr.norm[[2]],
                                    smoothers)
  }
  return(list(fun_arr = fun.arr.list, fun_norm = fun.norm.arr, part_res = part.res.mat.list, 
              lambda = lambda.seq, gcv.vec = gcv.vec, final.max.lam = final.max.lam))
}	

GSMALR <- function(resp.mat, pred.arr, alpha, lambda, 
                     FUN = GaussianKernel, final.max.lam = NULL, 
                     max.it = 1000, tol = 10^(-4), max.lam = 4,
                     max.lam.it = 50, max.lam.tol = 10^(-6), init.val = 0.01){

  fun.arr.list <- list()
  part.res.mat.list <- list()
  smoothers <- ArrayToKernel(pred.arr, FUN)

  # We use the parameter max.lam determined by the method without sparsity 
  # within groups by default to save some time
  num.lam <- length(lambda.seq)
  num.obs <- dim(resp.mat)[1]
  num.pred <- dim(pred.arr)[2]
  num.resp <- dim(resp.mat)[2]
  fun.norm.mat <- NULL
  total <- num.pred * num.obs * num.resp
  fun.norm.arr <- array(0, dim= c(num.resp, num.pred, num.lam))
  # initialize with random values 
  init.fun.arr <- array(runif(total, -init.val, init.val), dim= c(num.obs, num.pred, num.resp))
  gcv.vec <- rep(0, num.lam)
  intercepts <- matrix(colSums(resp.mat) / (num.obs - sum(sum(resp.mat))), 	
                       nrow = 1)
  all.intercepts <- matrix(rep(intercepts, num.obs), num.obs, num.resp,
                           byrow = TRUE)
  num.it <- 0 
  conv.crit <- 1 
  fun.arr <- array(0, dim = c(num.obs, num.pred ,num.resp))
  while(((conv.crit > tol) & (num.it < max.it))){
     old.fun.arr <- fun.arr
     sum.add.fun <- apply(fun.arr, MARGIN = c(1, 3), FUN = sum)
     prob.mat.num <- exp(sum.add.fun + all.intercepts) 
     prob.mat.den <- 1 + rowSums(prob.mat.num)
     prob.mat <- prob.mat.num / matrix(rep(prob.mat.den,num.resp), 
                 nrow = num.obs, ncol = num.resp)  
     # Calculate the transformed responses 
     trans.resp <- sum.add.fun + all.intercepts + 4 *(resp.mat - prob.mat)  
     fun.arr.obj <- BackfittingSMTSpAMR(trans.resp, pred.arr, init.fun.arr, 
                                        smoothers, alpha, sqrt(2) * lambda, 
                                        max.it, tol)
    fun.arr <- fun.arr.obj[['add_fun']]
    intercepts <- colMeans(trans.resp)
    all.intercepts <- matrix(rep(intercepts, num.obs), num.obs, num.resp,
                           byrow = TRUE)  
     conv.crit <- sqrt(sum((old.fun.arr - fun.arr)^2))
     num.it <- num.it + 1
   }  
   return(list(fun_arr = fun_arr, fun_norms = fun.arr.obj[["fun_norms"]], 
               intercepts = intercepts))
}


GSMALR_C <- function(resp.mat, pred.arr, alpha, lambda.seq, 
                     FUN = GaussianKernel, final.max.lam = NULL, 
                     in.max.it = 1000, in.tol = 10^(-3), out.max.it = 2000, 
                     out.tol = 10^(-4),max.lam = 4,
                     max.lam.it = 50, max.lam.tol = 10^(-6), init.val = 10^(-5), 
                     bth.plugin = NULL){
## This function is used to fit the additive functions for the whole set 
## of values of the parameters alpha and  lambda
##   - resp.mat: matrix of dimension num.obs -by- num.resp that 
##               contains all the responses
##   - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
##   - alpha.seq: sequence of values for the weigthing parameter alpha
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
  smoothers <- ArrayToKernel(pred.arr, FUN, bth.plugin = bth.plugin)

  # We use the parameter max.lam determined by the method without sparsity 
  # within groups by default to save some time
  num.lam <- length(lambda.seq)
  num.obs <- dim(resp.mat)[1]
  num.pred <- dim(pred.arr)[2]
  num.resp <- dim(resp.mat)[2]
  fun.norm.mat <- NULL
  total <- num.pred * num.obs * num.resp
  fun.norm.arr <- array(0, dim= c(num.resp, num.pred, num.lam))
  # initialize with random values 
  #init.fun.arr <- array(runif(total, -init.val, init.val), dim= c(num.obs, num.pred, num.resp))
  init.fun.arr <- array(init.val, dim= c(num.obs, num.pred, num.resp))
  gcv.vec <- rep(0, num.lam)
  intercepts <- matrix(colSums(resp.mat) / (num.obs - sum(sum(resp.mat))), 	
                       nrow = 1)
  all.intercepts <- matrix(0, num.lam, num.resp)
  for (lam.idx in seq_len(num.lam)){
     lam <- lambda.seq[lam.idx]
     fun.arr.norm <- BackfittingGSMALR(resp.mat, pred.arr, init.fun.arr, 
                                       smoothers, intercepts, alpha, lam, in.max.it, in.tol, 
                                       out.max.it, out.tol)

     fun.arr.list[[lam.idx]] <-  fun.arr.norm[[1]]
     fun.norm.arr[,,lam.idx] <-  fun.arr.norm[['fun_norms']]
     #part.res.mat.list[[lam.idx]] <- fun.arr.norm[[3]]
     trans_resp <- fun.arr.norm[['trans_resp']]
     gcv.vec[lam.idx] <- ComputeGCV(trans_resp, fun.arr.norm[[1]], fun.arr.norm[["fun_norms"]],
                                    smoothers)
     all.intercepts[lam.idx,] <- fun.arr.norm[['intercepts']]
     
  }
  return(list(fun_arr = fun.arr.list, fun_norm = fun.norm.arr, intercepts = all.intercepts, 
              lambda = lambda.seq,
              gcv.vec = gcv.vec, num_it = fun.arr.norm[['num_iter']], 
              trans_resp = trans_resp, prob_mat = fun.arr.norm[['prob_mat']]))
}	


SMTSpAM_All_Alphas <- function(resp.mat, pred.arr, alpha.seq, lambda.seq = NULL, 
                               lam.scale = (1 - seq(0, 0.95, by = 0.05)), 
                               FUN = GaussianKernel, final.max.lam = NULL, 
                               max.it = 1000, tol = 10^(-4), max.lam = 4,
                               max.lam.it = 50, max.lam.tol = 10^(-6), 
                               init.val = 0.01, bth.plugin = NULL){
## This function is used to fit the additive functions for the whole set 
## of values of the parameters alpha and  lambda
##   - resp.mat: matrix of dimension num.obs -by- num.resp that 
##               contains all the responses
##   - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
##   - alpha.seq: sequence of values for the weigthing parameter alpha
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

   all.outputs <- list()
   num.alpha <- length(alpha.seq)
   if (is.null(lambda.seq)){ 
      num.lam   <- length(lam.scale)
   } else {
      num.lam <- length(lambda.seq)
   }
   gcv.mat   <- matrix(0, nrow = num.alpha, ncol = num.lam)
   for (alpha.idx in seq_len(num.alpha)){
       alpha <- alpha.seq[alpha.idx]
       all.outputs[[alpha.idx]] <- SMTSpAM_C(resp.mat, pred.arr, alpha, lambda.seq, 
                                   lam.scale, FUN, final.max.lam, max.it, tol, max.lam, 
                                   max.lam.it, max.lam.tol, init.val = init.val, bth.plugin = bth.plugin)
       final.max.lam <- all.outputs[[alpha.idx]][['final.max.lam']]
       gcv.mat[alpha.idx,] <- all.outputs[[alpha.idx]][['gcv.vec']]
   }
   lambda.seq <- all.outputs[[num.alpha]][['lambda']]
   return(list(all.models = all.outputs, gcv.mat = gcv.mat, lambda = lambda.seq))
}



#require(mvtnorm)
require(splines)
#require(fda)
#setwd("Code")
#source("Spam")
#setwd("..")


#GenerateCovariates <- function(num.obs, num.cov, l.b = -2.5, u.b = 2.5, 
#                               cor.par = 0.5){
## This function is used to generate a set of covariates
## via the uniform distribution
## Inputs:
##   - num.obs: number of observations for each covariate
##   - num.cov: number of covariates generated
##   - l.b: lower bound of the uniform distribution
##   - u.b: upper bound of the uniform distribution
##   - cor.par: parameter used to change the degree of correlation
##              between the predictors
#  cov.basis <- matrix(runif(num.obs * num.cov, l.b, u.b) , 
#                    nrow = num.obs, ncol = num.cov)
#  com.basis <- runif(num.obs, l.b, u.b)
#  cov.mat <- (cov.basis + cor.par * com.basis) / ( 1 + cor.par)
  #max.vec <- apply(cov.mat, MARGIN = 2, FUN = max)
  #min.vec <- apply(cov.mat, MARGIN = 2, FUN = min)
  #range.vec   <- max.vec - min.vec
  #range.mat <- matrix(rep(range.vec, num.obs), nrow = num.obs, ncol = num.cov, 
  #                byrow = TRUE)
  #min.mat <- matrix(rep(min.vec, num.obs), nrow = num.obs, ncol = num.cov, 
  #                byrow = TRUE)
  #pred.mat <- (cov.mat - min.mat) / (range.mat)
#  return(cov.mat)
#}

# Generate Multi-task models

#GenerateModel <- function(covariates, num.resp, variances, 
#                          idx){
# f1 <- exp(-5 * (covariates[,1])) - exp(-1)
# f1 <- f1 - mean(f1)
# f2 <- 3 * sin(exp((covariates[,2]))) - 2
# f2 <- f2 - mean(f2)
# f3 <- 3 * ((covariates[,3])- 0.2)^3 -((covariates[,3]) - 0.3)^2 
# f3 <- f3 - mean(f3)
# f4 <-  -(covariates[,4])
# f4 <- f4 - mean(f4)
# f5 <- (covariates[,5])
# f5 <- f5 - mean(f5)
# f6 <- ( 2 * (covariates[, 6]) - 1)^2 
# f6 <- f6 - mean(f6)
# f7 <- sin(2 * pi * (covariates[,7])) / 
#       (2 - sin(2 * pi * (covariates[, 7])))
# f7 <- f7 - mean(f7)
# x8 <- (covariates[,8])
# f8 <- 0.1 * sin(2 * pi * x8) + 0.2 * cos(2 * pi * x8) + 
#          0.3 * sin(2 * pi * x8)^2 + 0.4 * cos(2 * pi * x8) + 
#          0.5 * sin(2 * pi * x8)^3
# f8 <- f8 - mean(f8)
# errors <- rmvnorm(dim(covariates)[1], rep(0, num.resp), 
#                       diag(variances))
# responses <- matrix(0, nrow = dim(covariates)[1], ncol = num.resp)
# responses <- idx[1] * f1 + idx[2] * f2 + idx[3] * f3 + idx[4] * f4 + 
#              idx[5] * f5 + idx[6] * f6 + idx[7] * f7 + idx[8] * f8 + errors
#   all_additive <- cbind(idx[1] * f1, idx[2] * f2, idx[3] * f3, idx[4] * f4, 
#                         idx[5] * f5, idx[6] * f6, idx[7] * f7, idx[8] * f8)
#   return(list(responses = responses, additive = all_additive))
#}

#GenerateMultiResponses <- function(covariates, variances, weights){
## This function generates the responses used 
## for the simulations
## Inputs:
##   - covariates: a matrix of dimension num.obs -by- num.cov
##   - variances: a vector with the standard deviations of the 
##                errors associated with the responses
##   - weights: a matrix of dimension 3 -by- 4
#  num.obs <- dim(covariates)[1]
#  num.resp <- 3
#  resp.matrix <- matrix(0, nrow = num.obs, ncol = num.resp)
#  x1 <- (covariates[,1] + 2.5) / 5
#  x2 <- (covariates[,2] + 2.5) / 5
#  x3 <- covariates[,3]
#  x4 <- covariates[,4]
#  f1 <- 0.1 * sin(2 * pi * x1) + 0.2 * cos(2 * pi * x1) + 
#          0.3 * sin(2 * pi * x1)^2 + 0.4 * cos(2 * pi * x1) + 
#          0.5 * sin(2 * pi * x1)^3
#  f1 <- f1 - mean(f1)
  
#  f21 <- sin(2 * pi * x2) / (2 - sin(2 * pi * x2))
#  f21 <- f21 - mean(f21)
#  f22 <- - f21
#  f23 <- (2 * x2 - 1)^2 
#  f23 <- f23 - mean(f23)
  
#  f31 <- - x3
#  f31 <- f31 - mean(f31)

#  f32 <- exp(-5 * x3) - exp(-1)
#  f32 <- f32 - mean(f32)

#  f33 <- exp(-3 * x3) - exp(-1)
#  f33 <- f33 - mean(f33)

#  f41 <- x4 
#  f41 <- f41 - mean(f41)
  
#  f42 <- -(2 * x4 - 1)^2
#  f42 <- f42 - mean(f42)
 
#  f43 <- -(x4 - 1)^3 - 1.5 * (x4 -1 )^2 + 0.25
#  f43 <- f43 - mean(f43)

#  resp.matrix[,1] <- weights[1, 1] * f1 + weights[1, 2] * f21 +
#                     weights[1, 3] * f31 + weights[1, 4] * f41 

#  resp.matrix[,2] <- weights[2, 1] * f1 + weights[2, 2] * f22 +
#                     weights[2, 3] * f32 + weights[2, 4] * f42 

#  resp.matrix[,3] <- weights[3, 1] * f1 + weights[3, 2] * f23 +
#                     weights[3, 3] * f33 + weights[3, 4] * f43 
   
#  errors <- rmvnorm(dim(covariates)[1], rep(0, num.resp), 
#                       diag(variances))
#  add_fun1 <- cbind(weights[1, 1] * f1 , weights[1, 2] * f21,
#                     weights[1, 3] * f31, weights[1, 4] * f41 )
#  add_fun2 <- cbind(weights[2, 1] * f1 , weights[2, 2] * f22 ,
#                     weights[2, 3] * f32 , weights[2, 4] * f42)
#  add_fun3 <- cbind(weights[3, 1] * f1 , weights[3, 2] * f23 ,
#                    weights[3, 3] * f33 , weights[3, 4] * f43 )

#  resp.matrix <- resp.matrix + errors
#  return(list(resp.matrix, add_fun1, add_fun2, add_fun3))
#}

#GenerateMultiResponses.2 <- function(covariates, variances, weights){
## This function generates the responses used 
## for the simulations
## Inputs:
##   - covariates: a matrix of dimension num.obs -by- num.cov
##   - variances: a vector with the standard deviations of the 
##                errors associated with the responses
##   - weights: a matrix of dimension 3 -by- 4
#  num.obs <- dim(covariates)[1]
#  num.resp <- 4
#  resp.matrix <- matrix(0, nrow = num.obs, ncol = num.resp)
#  x1 <- covariates[,1] 
#  x2 <- covariates[,2]
#  x3 <- covariates[,3]

#  f11 <- -5 * sin(2 *x1)
#  f11 <- f11 - mean(f11)


#  f12 <- x1^2
#  f12 <- f12 - mean(f12)


#  f13 <- (2 * sin(x1)) / (2 - sin(x1))
#  f13 <- f13 - mean(f13)


#  f14 <- exp(-x1)
#  f14 <- f14 - mean(f14)


#  f21 <- x2^3 + 1.5 * (x2 - 1)^2
#  f21 <- f21 - mean(f21)

  
#  f22 <- x2
#  f22 <- f22 - mean(f22)

   
#  f23 <- 3 * sin(exp(-.5 * x2))
#  f23 <- f23 - mean(f23)


#  f24 <- -5 * pnorm(x2, 0.5, 0.8)
#  f24 <- f24 - mean(f24)
  
#  x4 <- (covariates[,3] + 2.5) / 5
#  f31 <- 0.5 * sin(2 * pi * x4) + 1 * cos(2 * pi * x4) + 
#         1.5 * sin(2 * pi * x4)^2 + 2 * cos(2 * pi * x4)^3 + 
#         2.5 * sin(2 * pi * x4)^3
#  f31 <- f31 - mean(f31)
  
#  f32 <- 4 *(sin(2 * pi * x4) / (2 - sin(2 * pi * x4)))
#  f32 <- f32 - mean(f32)

#  f33 <- -x3
#  f33 <- f33 - mean(f33)

#  f34 <- - x3^2
#  f34 <- f34 - mean(f34)

#  resp.matrix[,1] <- weights[1, 1] * f11 + weights[1, 2] * f21 +
#                     weights[1, 3] * f31 
#  resp.matrix[,2] <- weights[2, 1] * f12 + weights[2, 2] * f22 +
#                     weights[2, 3] * f32
 
#  resp.matrix[,3] <- weights[3, 1] * f13 + weights[3, 2] * f23 +
#                     weights[3, 3] * f33

#  resp.matrix[,4] <- weights[4, 1] * f14 + weights[4, 2] * f24 +
#                     weights[4, 3] * f34


   
#  errors <- rmvnorm(dim(covariates)[1], rep(0, num.resp), 
#                       diag(variances))
#  add_fun1 <- cbind(weights[1, 1] * f11 , weights[1, 2] * f21, 
#                    weights[1, 3] * f31)
#  add_fun2 <- cbind(weights[2, 1] * f12 , weights[2, 2] * f22, 
#                    weights[2, 3] * f32)
#  add_fun3 <- cbind(weights[3, 1] * f13 , weights[3, 2] * f23,
#                    weights[3, 3] * f33)
#  add_fun4 <- cbind(weights[4, 1] * f14 , weights[4, 2] * f24, 
#                    weights[4, 3] * f34)

#  resp.matrix <- resp.matrix + errors
#  return(list(resp.matrix, add_fun1, add_fun2, add_fun3, add_fun4))
#}

functional <- function(responses, add_par, idx.fun , idx.resp = 1){
  all_add_idx <- rowSums(add_par) - add_par[, idx.fun]
  functional_idx <- responses[, idx.resp] - all_add_idx
  return(functional_idx)
}

CovariatesArray <- function(num.resp, covariates){
 n.obs <- dim(covariates)[1]
 n.pred <- dim(covariates)[2]
 cov.array <- array(0, dim = c(n.obs, n.pred, num.resp))
 for (idx in seq_len(num.resp)){
    cov.array[ , ,idx] <- covariates
 }
 return(cov.array)
}

GenerateTrueFunctions <- function(input.vec, weights, resp.idx = 1, 
                                  num.cov = 16){
## This function returns the functions that are used in the 
## simulations
  num.obs <- length(input.vec)
  x  <- input.vec
  f1 <- 0.1 * sin(2 * pi * x) + 0.2 * cos(2 * pi * x) + 
          0.3 * sin(2 * pi * x)^2 + 0.4 * cos(2 * pi * x) + 
          0.5 * sin(2 * pi * x)^3
  f1 <- f1 - mean(f1)
  f11 <- weights[1, 1]  * f1
  f12 <- weights[2, 1] * f1
  f13 <- weights[3, 1] * f1

  f2 <- sin(2 * pi * x) / (2 - sin(2 * pi * x))
  f21 <- f2 - mean(f2)
  f21 <- weights[1, 2] * f21

  f22 <- - f2
  f22 <- weights[2,2] * f22
  f23 <- (2 * x - 1)^2 
  f23 <- f23 - mean(f23)
  f23 <- weights[3, 2] * f23

  f31 <- - x
  f31 <- f31 - mean(f31)
  f31 <- weights[1, 3] * f31

  f32 <- exp(-5 * x) - exp(-1)
  f32 <- f32 - mean(f32)
  f32 <- weights[2, 3] * f32

  f33 <- exp(-3 * x) - exp(-1)
  f33 <- f33 - mean(f33)
  f33 <- weights[3, 3]  * f33

  f41 <- x 
  f41 <- f41 - mean(f41)
  f41 <- weights[1, 4] * f41  

  f42 <- -(2 * x - 1)^2
  f42 <- f42 - mean(f42)
  f42 <- weights[2, 4]  * f42

  f43 <- -(x - 1)^3 - 1.5 * (x -1 )^2 + 0.25
  f43 <- f43 - mean(f43)
  f43 <- weights[3, 4] * f43

  fun.mat <- matrix(0, nrow = num.obs, ncol = num.cov)
  if (resp.idx == 1){
    fun.mat[,1:4] <- cbind(f11, f21, f31, f41)
  } 
  if (resp.idx == 2){
    fun.mat[,1:4] <- cbind(f12, f22, f32, f42)
  }
  if (resp.idx == 3) {
    fun.mat[,1:4] <- cbind(f13, f23, f33, f43)
  }
  return(fun.mat)
}

GenerateTrueFunctions.2 <- function(input.vec, weights, resp.idx = 1, 
                                  num.cov = 16){
## This function returns the functions that are used in the 
## simulations
  num.obs <- length(input.vec)
  x  <- input.vec
  f11 <- - 5 * sin(2 *x)
  f11 <- f11 - mean(f11)
  f11 <- weights[1, 1] * f11

  f12 <- x^2
  f12 <- f12 - mean(f12)
  f12 <- weights[2,1] * f12

  f13 <- (2 * sin(x)) / (2 - sin(x))
  f13 <- f13 - mean(f13)
  f13 <- weights[3, 1] * f13

  f14 <- exp(-x)
  f14 <- f14 - mean(f14)
  f14 <- weights[4, 1] * f14

  f21 <- x^3 + 1.5 * (x - 1)^2
  f21 <- f21 - mean(f21)
  f21 <- weights[1, 2] * f21
  
  f22 <- x
  f22 <- f22 - mean(f22)
  f22 <- weights[2, 2] * f22 
   
  f23 <- 3 * sin(exp(-.5 * x))
  f23 <- f23 - mean(f23)
  f23 <- weights[3, 2] * f23

  f24 <- -5 * pnorm(x, 0.5, 0.8)
  f24 <- f24 - mean(f24)
  f24 <- weights[4, 2] * f24

  y <- (x + 2.5) / 5
  f31 <- 0.5 * sin(2 * pi * y) + 1 * cos(2 * pi * y) + 
         1.5 * sin(2 * pi * y)^2 + 2 * cos(2 * pi * y)^3 + 
         2.5 * sin(2 * pi * y)^3
  f31 <- f31 - mean(f31)
  f31 <- weights[1, 3] * f31

  f32 <- 4 * (sin(2 * pi * y) / (2 - sin(2 * pi * y)))
  f32 <- f32 - mean(f32)
  f32 <- weights[2, 3] * f32

  f33 <- -x
  f33 <- f33 - mean(f33)
  f33 <- weights[3, 3] * f33

  f34 <- - x^2
  f34 <- f34 - mean(f34)
  f34 <- weights[4, 3] * f34

  fun.mat <- matrix(0, nrow = num.obs, ncol = num.cov)
  if (resp.idx == 1){
    fun.mat[,1:3] <- cbind(f11, f21, f31)
  } 
  if (resp.idx == 2){
    fun.mat[,1:3] <- cbind(f12, f22, f32)
  }
  if (resp.idx == 3){
    fun.mat[,1:3] <- cbind(f13, f23, f33)
  } 
  if (resp.idx == 4){
    fun.mat[,1:3] <- cbind(f14, f24, f34)
  }
  return(fun.mat)
}


#ConstructSplineBasis <- function(covariate, num.knots, norder){
# knots <- quantile(covariate, seq(0, 1, length = num.knots)) 
# basis <- create.bspline.basis(rangeval = range(knots), breaks = knots, 
#                              norder = norder)
# design.cov <- eval.basis(covariate, basis)
# rough.penalty <- getbasispenalty(basis, Lfdobj = 2)
# return(list(design.cov = design.cov, penalty = rough.penalty))
#}

#ConstructAllBasis <- function(covariates, num.knots = 10, norder = 4){
#  num.cov <-  dim(covariates)[2]
#  nbasis <- num.knots + norder - 2
#  B <- matrix(0, nrow = dim(covariates)[1],
#                 ncol = nbasis * num.cov)
#  Omega <- matrix(0, nrow = nbasis,
#                 ncol = nbasis*num.cov)
# 
#  beg.idx <- 0
#  end.idx <- 0
#  for (cov.idx in seq_len(num.cov)){
#    beg.idx <- end.idx + 1
#    end.idx <- end.idx + nbasis 
#     b_pen<- ConstructSplineBasis(covariates[, cov.idx],
#                            num.knots = num.knots, norder = norder)
 #    B[,beg.idx: end.idx] <- b_pen[[1]]
#     Omega[,beg.idx: end.idx] <- b_pen[[2]]

#  } 
#  return(list(B = B, Omega = Omega))
#}


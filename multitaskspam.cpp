// [[Rcpp::depends(RcppArmadillo)]]
#include <algorithm>
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

mat updateMTPartialResidual(mat part_res, cube add_fun_arr, int pred_idx){
// This function removes the function f_j^{(k)} from the response y^{(k)} 
// Inputs:
//   - part_res: a matrix of partial residuals (n_obs -by- n_resp)
//   - add_fun_arr: 3D array containing all the additive functions 
//                  n_obs -by- n_pred -by- n_resp
//   - pred_idx: predictor index that starts at 0
// Output: returns an update partial residual matrix 
  mat add_fun_j = add_fun_arr.subcube(span(), span(pred_idx), span());
  mat new_part_res = part_res - add_fun_j;
  return new_part_res;
}

// [[Rcpp::export]]
mat updateMTPartialResidualR(NumericMatrix part_res_r, NumericVector add_fun_arr_r, int pred_idx){
  IntegerVector arrayDims = add_fun_arr_r.attr("dim");
  cube add_fun_arr(add_fun_arr_r.begin(), arrayDims[0], arrayDims[1], arrayDims[2], arrayDims[3], false);
  mat part_res(part_res_r.begin(), part_res_r.nrow(), part_res_r.ncol() ,false);
  //mat add_fun_mat(add_fun_matr.begin(), add_fun_matr.nrow(), add_fun_matr.ncol() ,false)  ;
  return updateMTPartialResidual(part_res, add_fun_arr, (pred_idx - 1));
}


mat updateMTAdditiveComponent( mat kernel_mat, mat part_res, double lam){
 /* 
  * This function updates all the additive functions associated  with a 
  * the set of predictors X_j^{(1)}, ...., X_j^{(K)}
  * Inputs:
  *   - kernel_mat: a matrix of dimension n_obs -by- n_obs * n_resp
  *                 all the entries are related to the same type of predictors
  *   - part_res: a matrix of dimension n_obs -by- n_resp that contains 
  *               the partial residuals for all the responses
  *   - lambda: regularization parameter
  * Outputs: returns the matrix of updated additive functions 
  *          [f_j(1), ..., f_j(K)] of dimension n_obs -by- n_resp
  *          
  */
  double f_n, W_j;
  int  n_obs = part_res.n_rows, n_resp = part_res.n_cols, beg_idx = 0;
  mat add_fun_mat, add_fun_mean;
  add_fun_mat.zeros(n_obs, n_resp);
  for (int resp_idx = 0; resp_idx < n_resp; resp_idx++){
	  add_fun_mat.col(resp_idx) = kernel_mat.cols(beg_idx, beg_idx + n_obs - 1) *
		                           part_res.col(resp_idx) ;
      beg_idx += n_obs;
  }

  mat sq_add_fun_mat = add_fun_mat % add_fun_mat;
  rowvec fun_norms = sum(sq_add_fun_mat, 0);
  f_n = sqrt(sum(fun_norms) / n_obs);
  vec y = zeros<vec>(2);
  y(0) = 1 - ((lam * sqrt(n_resp)) / f_n);
  W_j = max(y);
  add_fun_mat *= W_j;
  //rowvec fun_mean_vec = mean(add_fun_mat, 0);
  //add_fun_mat -= repmat(fun_mean_vec, n_obs, 1);
  return add_fun_mat;
}

// [[Rcpp::export]]
mat updateMTAdditiveComponent(NumericMatrix kernel_mat_r, NumericMatrix part_res_r, double lam){
  mat part_res(part_res_r.begin(), part_res_r.nrow(), part_res_r.ncol() ,false);
  mat kernel_mat(kernel_mat_r.begin(), kernel_mat_r.nrow(), kernel_mat_r.ncol() ,false);
 return updateMTAdditiveComponent(kernel_mat, part_res, lam);
}

mat updateMTSAdditiveComponent(mat kernel_mat, mat part_res, mat current_fun, double alpha, double lambda){
/* 
  * This function updates all the additive functions associated  with a 
  * the set of predictors X_j^{(1)}, ...., X_j^{(K)}
  * Inputs:
  *   - kernel_mat: a matrix of dimension n_obs -by- n_obs * n_resp
  *                 all the entries are related to the same type of predictors
  *   - part_res: a matrix of dimension n_obs -by- n_resp that contains 
  *               the partial residuals for all the responses
  *   - current_fun: a matrix of dimension n_obs -by- n_resp 
  *   - lambda: regularization parameter
  *   - alpha: trade off regularization 
  *   - current_fun: matrix of dimension n_obs -by- n_resp with the 
  *                  [f_j(1), ... ,f_j(K)] at time t
  * Outputs: returns the matrix of updated additive functions 
  *          [f_j(1), ..., f_j(K)] of dimension n_obs -by- n_resp
  *          at time t+1
  */
	int n_obs = part_res.n_rows, n_resp = part_res.n_cols, beg_idx = 0;
	mat add_fun_mat;
	add_fun_mat.zeros(n_obs, n_resp);
	for (int resp_idx = 0; resp_idx < n_resp; resp_idx++){
	   add_fun_mat.col(resp_idx) = kernel_mat.cols(beg_idx, beg_idx + n_obs - 1) * 
		                           part_res.col(resp_idx);
	   beg_idx += n_obs;
	}
	mat sq_add_fun = add_fun_mat % add_fun_mat;
	rowvec fun_norms = sum(sq_add_fun, 0);
	//rowvec zero_vec = zeros<rowvec>(n_resp);
	//rowvec g = max(zero_vec, (1 - alpha / fun_norms));
	mat mat_fun_n; 
	mat_fun_n.zeros(2, n_resp);
	mat_fun_n.row(1) = 1 - (alpha * n_obs) / fun_norms;
	rowvec g = max(mat_fun_n, 0);
	rowvec h_vec = square(g % sqrt(fun_norms));
	double h = sqrt(sum(h_vec) / n_obs);
	mat output_fun; 
	output_fun.zeros(n_obs, n_resp);
	double group_threshold = (1 - alpha) * lambda * sqrt(n_resp) ;
    
	if (h > group_threshold) {
	   rowvec curr_fun_norm = sum(current_fun % current_fun, 0) / n_obs ;
	   double curr_fun_gr_norm = sqrt(sum(curr_fun_norm)), al_lam = alpha * lambda ;

	   // Update each function if necessary
	   for (int resp_idx = 0; resp_idx < n_resp; resp_idx++) {
		   if (sqrt(fun_norms(resp_idx) / n_obs) > al_lam) {
		      output_fun.col(resp_idx) = add_fun_mat.col(resp_idx) / 
	             (1 + (group_threshold / curr_fun_gr_norm) + (al_lam / sqrt(curr_fun_norm(resp_idx))));
		   }
	   }
	} 

	//rowvec fun_mean_vec = mean(output_fun, 0);
	//output_fun -= repmat(fun_mean_vec, n_obs, 1);
	return output_fun; 
	//return h_vec;
}



// [[Rcpp::export]]
mat updateMTSAdditiveComponent(NumericMatrix kernel_mat_r, NumericMatrix part_res_r, NumericMatrix current_fun_r, double alpha, double lambda){
  mat part_res(part_res_r.begin(), part_res_r.nrow(), part_res_r.ncol() ,false);
  mat kernel_mat(kernel_mat_r.begin(), kernel_mat_r.nrow(), kernel_mat_r.ncol() ,false);
  mat current_fun(current_fun_r.begin(), current_fun_r.nrow(), current_fun_r.ncol() ,false);
 return updateMTSAdditiveComponent(kernel_mat, part_res, current_fun, alpha, lambda);
}

mat computeFunNorm(cube fun_norm_arr){
	int n_pred = fun_norm_arr.n_cols,
		n_resp = fun_norm_arr.n_slices;
	mat temp_fun_mat, fun_norm_mat = zeros<mat>(n_resp, n_pred);
	
	for (int pred_idx = 0; pred_idx < n_pred; pred_idx++){
		temp_fun_mat = fun_norm_arr.subcube( span(), span(pred_idx), span());
		fun_norm_mat.col(pred_idx) = sqrt(sum(pow(temp_fun_mat, 2), 0).t());
	}
	return fun_norm_mat;
	//return temp_fun_mat;
}

cube BackfittingMTSpAM(mat resp_mat, cube pred_arr, cube init_fun_arr, 
	                   cube smoother_arr, double lam, int max_it, double tol){
/* This function returns the estimated additive functions in the 
   the multi-task (or multivariate) regression setting
   Input:
      - resp.mat: matrix of dimension num.obs -by- num.resp that 
                contains all the responses
     - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
     - init.fun.arr: array of dimension num.obs - by- num.pred -by- num.resp
                     that contains the additive function values
    - smoother_arr: array of dim n-obs -by- n_obs * n_pred -by- 
    - lambda: regularization parameter value
    - max.it: maximum number of iteration for the outerloop 
    - tol: tolerance level for the convergence criterion
  Outputs: 
*   - fun.arr: returns an array of dimension  num.obs -by- num.pred -by-
*              num.resp with the estimated additive functions 
*/
  cube fun_arr = init_fun_arr;
  mat part_res, smoother_mat;
  cube old_fun_arr; 
  //vec it ;
  double conv_crit = 1;
  int n_pred = pred_arr.n_cols, n_resp = resp_mat.n_cols, num_it = 1;
  part_res = resp_mat - repmat(mean(resp_mat, 0), resp_mat.n_rows, 1);
  while ((num_it < max_it) && (conv_crit > tol)){
     old_fun_arr = fun_arr;
	 //= randi<vec>(100,distr_parameter(0, 99));
	 for (int pred_idx = 0; pred_idx < n_pred; pred_idx++){
	    part_res = updateMTPartialResidual(part_res, -fun_arr, pred_idx);
		smoother_mat = smoother_arr.slice(pred_idx);
		fun_arr.subcube(span(), span(pred_idx), span())  = updateMTAdditiveComponent(smoother_mat, part_res, lam);
		part_res = updateMTPartialResidual(part_res, fun_arr, pred_idx);
	 }
	 conv_crit = 0;
	 for (int resp_idx = 0; resp_idx < n_resp; resp_idx ++){
		 conv_crit += sqrt(pow(accu(old_fun_arr.slice(resp_idx) - fun_arr.slice(resp_idx)), 2)) ;
	 }
	 
  }
  return fun_arr;
}

// [[Rcpp::export]]
List BackfittingMTSpAMR(NumericMatrix resp_mat_r, NumericVector pred_arr_r, NumericVector init_fun_arr_r, 
                    	NumericVector smoother_arr_r, double lam, int max_it, double tol){

 // Conversion of R inputs into C++ inputs
 IntegerVector pred_dim = pred_arr_r.attr("dim"), fun_arr_dim = init_fun_arr_r.attr("dim");
 IntegerVector smoother_dim = smoother_arr_r.attr("dim");
 mat resp_mat(resp_mat_r.begin(), resp_mat_r.nrow(), resp_mat_r.ncol(), false);
 //cube add_fun_arr; mat add_fun_norms;
 cube pred_arr(pred_arr_r.begin(), pred_dim[0], pred_dim[1], pred_dim[2], false), 
	  init_fun_arr(init_fun_arr_r.begin(), fun_arr_dim[0], fun_arr_dim[1], fun_arr_dim[2]),
	  smoother_arr(smoother_arr_r.begin(), smoother_dim[0], smoother_dim[1], smoother_dim[2]);

  cube fun_arr = init_fun_arr;
  mat part_res, smoother_mat, add_fun_norms;
  cube old_fun_arr;
  double conv_crit = 1;
  int n_pred = pred_arr.n_cols, n_resp = resp_mat.n_cols, num_it = 1;
  part_res = resp_mat - repmat(mean(resp_mat, 0), resp_mat.n_rows, 1);
  while ((num_it < max_it) && (conv_crit > tol)){
     old_fun_arr = fun_arr;
	 for (int pred_idx = 0; pred_idx < n_pred; pred_idx++){
	    part_res = updateMTPartialResidual(part_res, -fun_arr, pred_idx);
		smoother_mat = smoother_arr.slice(pred_idx);
		fun_arr.subcube(span(), span(pred_idx), span())  = updateMTAdditiveComponent(smoother_mat, part_res, lam);
		part_res = updateMTPartialResidual(part_res, fun_arr, pred_idx);
	 }
	 conv_crit = 0;
	 for (int resp_idx = 0; resp_idx < n_resp; resp_idx ++){
		 conv_crit += sqrt(pow(accu(old_fun_arr.slice(resp_idx) - fun_arr.slice(resp_idx)), 2)) ;
	 }
	 
  }
 //add_fun_arr = BackfittingMTSpAM(resp_mat, pred_arr, init_fun_arr, smoother_arr, lam, max_it, tol);
 add_fun_norms = computeFunNorm(fun_arr);
 return List::create(Named("add_fun") = fun_arr, Named("fun_norms") = add_fun_norms,
	                 Named("part_residuals") = part_res);
}


cube BackfittingSMTSpAM(mat resp_mat, cube pred_arr, cube init_fun_arr, 
	cube smoother_arr, double alpha, double lambda, int max_it, double tol){

/* This function returns the estimated additive functions in the 
   the multi-task (or multivariate) regression setting
   Input:
      - resp.mat: matrix of dimension num.obs -by- num.resp that 
                contains all the responses
     - pred.arr: array of dimension num.obs -by- num.pred -by- num.resp
     - init.fun.arr: array of dimension num.obs - by- num.pred -by- num.resp
                     that contains the additive function values
    - smoother_arr: array of dim n-obs -by- n_obs * n_pred -by- 
	- alpha: weighting parameter between group penalty and single penalty
    - lambda: regularization parameter value
    - max.it: maximum number of iteration for the outerloop 
    - tol: tolerance level for the convergence criterion
  Outputs: 
*   - fun.arr: returns an array of dimension  num.obs -by- num.pred -by-
*              num.resp with the estimated additive functions 
*/
  cube fun_arr = init_fun_arr;
  mat part_res, smoother_mat, add_fun_norms;
  cube old_fun_arr;
  double conv_crit = 1;
  int n_obs = resp_mat.n_rows,n_pred = pred_arr.n_cols, n_resp = resp_mat.n_cols, num_it = 1;
  part_res = resp_mat - repmat(mean(resp_mat, 0), resp_mat.n_rows, 1);
  mat old_fun_mat; 
  old_fun_mat.zeros(n_obs, n_resp);
  while ((num_it < max_it) && (conv_crit > tol)){
     old_fun_arr = fun_arr;
	 for (int pred_idx = 0; pred_idx < n_pred; pred_idx++){
	    part_res = updateMTPartialResidual(part_res, -fun_arr, pred_idx);
		smoother_mat = smoother_arr.slice(pred_idx);
		old_fun_mat = old_fun_arr.subcube(span(), span(pred_idx), span());
		fun_arr.subcube(span(), span(pred_idx), span())  = updateMTSAdditiveComponent(smoother_mat, part_res, old_fun_mat, alpha, lambda);
		part_res = updateMTPartialResidual(part_res, fun_arr, pred_idx);
	 }
	 conv_crit = 0;
	 for (int resp_idx = 0; resp_idx < n_resp; resp_idx ++){
		 conv_crit += sqrt(pow(accu(old_fun_arr.slice(resp_idx) - fun_arr.slice(resp_idx)), 2)) ;
	 }
	 num_it++;
  }
  return fun_arr;
 //add_fun_arr = BackfittingMTSpAM(resp_mat, pred_arr, init_fun_arr, smoother_arr, lam, max_it, tol);
 //add_fun_norms = computeFunNorm(fun_arr);
 //return List::create(Named("add_fun") = fun_arr, Named("fun_norms") = add_fun_norms,
	//                 Named("part_residuals") = part_res);
 }


// [[Rcpp::export]]
List BackfittingSMTSpAMR(NumericMatrix resp_mat_r, NumericVector pred_arr_r, NumericVector init_fun_arr_r, 
	NumericVector smoother_arr_r, double alpha, double lambda, int max_it, double tol){
// Conversion of R inputs into C++ inputs
 IntegerVector pred_dim = pred_arr_r.attr("dim"), fun_arr_dim = init_fun_arr_r.attr("dim");
 IntegerVector smoother_dim = smoother_arr_r.attr("dim");
 mat resp_mat(resp_mat_r.begin(), resp_mat_r.nrow(), resp_mat_r.ncol(), false);
 //cube add_fun_arr; mat add_fun_norms;
 cube pred_arr(pred_arr_r.begin(), pred_dim[0], pred_dim[1], pred_dim[2], false), 
	  init_fun_arr(init_fun_arr_r.begin(), fun_arr_dim[0], fun_arr_dim[1], fun_arr_dim[2]),
	  smoother_arr(smoother_arr_r.begin(), smoother_dim[0], smoother_dim[1], smoother_dim[2]);

  cube fun_arr = init_fun_arr;
  mat part_res, smoother_mat, add_fun_norms;
  cube old_fun_arr;
  double conv_crit = 1;
  int n_obs = resp_mat.n_rows,n_pred = pred_arr.n_cols, n_resp = resp_mat.n_cols, num_it = 1;
  part_res = resp_mat - repmat(mean(resp_mat, 0), resp_mat.n_rows, 1);
  mat old_fun_mat; 
  old_fun_mat.zeros(n_obs, n_resp);
  while ((num_it < max_it) && (conv_crit > tol)){
     old_fun_arr = fun_arr;
	 for (int pred_idx = 0; pred_idx < n_pred; pred_idx++){
	    part_res = updateMTPartialResidual(part_res, -fun_arr, pred_idx);
		smoother_mat = smoother_arr.slice(pred_idx);
		old_fun_mat = old_fun_arr.subcube(span(), span(pred_idx), span());
		fun_arr.subcube(span(), span(pred_idx), span())  = updateMTSAdditiveComponent(smoother_mat, part_res, old_fun_mat, alpha, lambda);
		part_res = updateMTPartialResidual(part_res, fun_arr, pred_idx);
	 }
	 conv_crit = 0;
	 for (int resp_idx = 0; resp_idx < n_resp; resp_idx ++){
		 conv_crit += sqrt(pow(accu(old_fun_arr.slice(resp_idx) - fun_arr.slice(resp_idx)), 2)) ;
	 }
         num_it++;	 
  }
 //add_fun_arr = BackfittingMTSpAM(resp_mat, pred_arr, init_fun_arr, smoother_arr, lam, max_it, tol);
 add_fun_norms = computeFunNorm(fun_arr);
 return List::create(Named("add_fun") = fun_arr, Named("fun_norms") = add_fun_norms,
	                 Named("part_residuals") = part_res);
 }


// [[Rcpp::export]]
List BackfittingGSMALR(NumericMatrix resp_mat_r, NumericVector pred_arr_r, NumericVector init_fun_arr_r, 
	               NumericVector smoother_arr_r, NumericMatrix intercepts_r, 
	               double alpha, double lambda, int in_max_it, double in_tol, int out_max_it , double out_tol){
// Conversion of R inputs into C++ inputs
 IntegerVector pred_dim = pred_arr_r.attr("dim"), fun_arr_dim = init_fun_arr_r.attr("dim");
 IntegerVector smoother_dim = smoother_arr_r.attr("dim");
 mat resp_mat(resp_mat_r.begin(), resp_mat_r.nrow(), resp_mat_r.ncol(), false);
 mat intercepts(intercepts_r.begin(), intercepts_r.nrow(), intercepts_r.ncol(), false);
 //cube add_fun_arr; mat add_fun_norms;
 cube pred_arr(pred_arr_r.begin(), pred_dim[0], pred_dim[1], pred_dim[2], false), 
	  init_fun_arr(init_fun_arr_r.begin(), fun_arr_dim[0], fun_arr_dim[1], fun_arr_dim[2]),
	  smoother_arr(smoother_arr_r.begin(), smoother_dim[0], smoother_dim[1], smoother_dim[2]);

  
  mat part_res, smoother_mat, add_fun_norms;
  cube old_fun_arr;
  double conv_crit = 1;
  int n_obs = resp_mat.n_rows,n_pred = pred_arr.n_cols, n_resp = resp_mat.n_cols, num_it = 1;
  part_res = resp_mat - repmat(mean(resp_mat, 0), resp_mat.n_rows, 1);
  //mat old_fun_mat; 
  //old_fun_mat.zeros(n_obs, n_resp);
  cube fun_arr = zeros<cube>(n_obs, n_pred, n_resp);
  mat exp_additive_sums = zeros<mat>(n_obs, n_resp);
  mat trans_resp = zeros<mat>(n_obs, n_resp);
  mat prob_mat = zeros<mat>(n_obs, n_resp);
  colvec inte = ones <colvec>(n_obs);
  while ((num_it < out_max_it) && (conv_crit > out_tol)){
	old_fun_arr = fun_arr;
    for (int resp_idx = 0; resp_idx < n_resp; resp_idx++){
       exp_additive_sums.col(resp_idx) = sum(fun_arr.slice(resp_idx), 1);  
	   exp_additive_sums.col(resp_idx).each_row() += intercepts.col(resp_idx);
	}
	exp_additive_sums = exp(exp_additive_sums);
    colvec denominator = inte + sum(exp_additive_sums, 1);
	for (int resp_idx = 0; resp_idx < n_resp; resp_idx++){
		prob_mat.col(resp_idx) = exp_additive_sums.col(resp_idx) / denominator;
	}
	for (int resp_idx = 0; resp_idx < n_resp; resp_idx++){
	  //trans_resp.col(resp_idx) = intercepts.col(resp_idx) * inte + 
	  //	                         4 * (resp_mat.col(resp_idx) - prob_mat.col(resp_idx)) +
      //                           sum(fun_arr.slice(resp_idx), 1); 
	  trans_resp.col(resp_idx) = 4 * (resp_mat.col(resp_idx) - prob_mat.col(resp_idx)) +
                                 sum(fun_arr.slice(resp_idx), 1); 

	  trans_resp.col(resp_idx).each_row() += intercepts.col(resp_idx);
	}
	fun_arr = BackfittingSMTSpAM(trans_resp, pred_arr, init_fun_arr, 
		                         smoother_arr, alpha, sqrt(2) * lambda, in_max_it, in_tol);
	intercepts = mean(trans_resp, 0);
	conv_crit = 0;
	 for (int resp_idx = 0; resp_idx < n_resp; resp_idx ++){
		 conv_crit += sqrt(pow(accu(old_fun_arr.slice(resp_idx) - fun_arr.slice(resp_idx)), 2)) ;
	 }
	num_it++;
  }
  add_fun_norms = computeFunNorm(fun_arr);
  return List::create(Named("add_fun") = fun_arr, Named("intercepts") = intercepts,
	                  Named("fun_norms") = add_fun_norms, Named("trans_resp") = trans_resp, 
					  Named("prob_mat") = prob_mat,
					  Named("conv_crit") = conv_crit, Named("num_iter") = num_it);
}
	               
// [[Rcpp::export]]
double maxofmatrix(NumericVector arr_r){
	IntegerVector arr_dim = arr_r.attr("dim");
	cube arr(arr_r.begin(), arr_dim[0], arr_dim[1], arr_dim[2], false);
	mat matrix =  computeFunNorm(arr);
	double max_mean;
	max_mean = max(max(mean(matrix, 0)));
	return max_mean; 
}

//[[Rcpp::export]]
double MTLambdaMaxR(NumericMatrix resp_mat_r, NumericVector pred_arr_r, NumericVector smoother_arr_r, 
	double max_lam, int max_it, double tol, int max_lam_it, double max_lam_tol){
  
  IntegerVector pred_dim = pred_arr_r.attr("dim");
  IntegerVector smoother_dim = smoother_arr_r.attr("dim");
  mat resp_mat(resp_mat_r.begin(), resp_mat_r.nrow(), resp_mat_r.ncol(), false);
  cube pred_arr(pred_arr_r.begin(), pred_dim[0], pred_dim[1], pred_dim[2], false), 
  smoother_arr(smoother_arr_r.begin(), smoother_dim[0], smoother_dim[1], smoother_dim[2]);

  double upper_lam = max_lam, lower_lam = 0, est_lam = upper_lam;
  double gap = upper_lam - lower_lam, max_fun_norm;
  cube init_fun_arr = zeros<cube>(pred_dim[0], pred_dim[1], pred_dim[2]);
  int it = 1;
  cube add_fun_arr; mat fun_norm_mat;
  while ((gap > max_lam_tol) & (it < max_lam_it)){
     add_fun_arr  = BackfittingMTSpAM(resp_mat, pred_arr, init_fun_arr, smoother_arr, 
		                          est_lam, max_it, tol) ;
	 fun_norm_mat = computeFunNorm(add_fun_arr);
	 max_fun_norm = max(max(mean(fun_norm_mat, 0)));
	 if (max_fun_norm == 0){
		 upper_lam = est_lam;
		 est_lam = upper_lam - (upper_lam - lower_lam) / 2;
	 } else {
	     lower_lam = est_lam;
		 est_lam = lower_lam + (upper_lam - lower_lam)/2;
	 }
	 gap = upper_lam - lower_lam; 
	 it++;
  }
  return upper_lam;
} 

// Check this function
//[[Rcpp::export]]
double SMTLambdaMaxR(NumericMatrix resp_mat_r, NumericVector pred_arr_r, NumericVector smoother_arr_r, 
	NumericVector init_fun_arr_r, double max_lam, int max_it, double tol, int max_lam_it, double max_lam_tol, double alpha){
  
  IntegerVector pred_dim = pred_arr_r.attr("dim");
  IntegerVector smoother_dim = smoother_arr_r.attr("dim");
  IntegerVector init_dim = init_fun_arr_r.attr("dim");
  mat resp_mat(resp_mat_r.begin(), resp_mat_r.nrow(), resp_mat_r.ncol(), false);
  cube pred_arr(pred_arr_r.begin(), pred_dim[0], pred_dim[1], pred_dim[2], false), 
  smoother_arr(smoother_arr_r.begin(), smoother_dim[0], smoother_dim[1], smoother_dim[2]),
  init_fun_arr(init_fun_arr_r.begin(), init_dim[0], init_dim[1], init_dim[2]);
  double upper_lam = max_lam, lower_lam = 0, est_lam = upper_lam;
  double gap = upper_lam - lower_lam, max_fun_norm;
  // bad initialization !!!! randomize this 
  // cube init_fun_arr = (randu<cube>(pred_dim[0], pred_dim[1], pred_dim[2]) - 0.5) / 5 ;
  int it = 1;
  cube add_fun_arr; mat fun_norm_mat;
  while ((gap > max_lam_tol) & (it < max_lam_it)){
     add_fun_arr  = BackfittingSMTSpAM(resp_mat, pred_arr, init_fun_arr, smoother_arr, 
		                          alpha, est_lam, max_it, tol) ;
	 fun_norm_mat = computeFunNorm(add_fun_arr);
	 max_fun_norm = max(max(mean(fun_norm_mat, 0)));
	 if (max_fun_norm == 0){
		 upper_lam = est_lam;
		 est_lam = upper_lam - (upper_lam - lower_lam) / 2;
	 } else {
	     lower_lam = est_lam;
		 est_lam = lower_lam + (upper_lam - lower_lam)/2;
	 }
	 gap = upper_lam - lower_lam; 
	 it++;
  }
  return upper_lam;
  //return max_fun_norm;
} 

This file contains a brief description of the R files containing the implementation of the methods 
presented in the paper High-Dimensional Multivariate Additive Regression for Uncovering Contributing 
Factors to Healthcare Expenditures. 

1- Utilies.R: contains multiple utilities functions that compute the linear smoother matrices associated 
              with each predictor / functions for plotting regularization paths and others. 
              For each function the inputs and outputs are described (types)

2- MTSpAM.R: contains implementations of the outer layers of algorithm 1 and 3 for different values of 
             the regularization parameters.

3- multitaskspam.cpp: contains C++ routines used to compute the inner looops of the algorithm 1 and 2
                      (updating each additive function)

4- prediction_functions.R: contains predictor functions and cross validation functions 

5- generatesimulations.R: contains the R code used to generate simulation functions presented in 
                          Appendix of the paper. 
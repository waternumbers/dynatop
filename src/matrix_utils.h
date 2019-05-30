// header file for the matrix utility functions
#ifndef MATUTILS
#define MATUTILS
// prevents multiple definitions
#include <Rcpp.h>

Rcpp::NumericMatrix mmult(Rcpp::NumericMatrix a, Rcpp::NumericMatrix b);
Rcpp::NumericVector matrix_vect_mult_cpp(Rcpp::NumericMatrix A, Rcpp::NumericVector b);
Rcpp::NumericVector vect_matrix_mult_cpp(Rcpp::NumericVector b, Rcpp::NumericMatrix A);


#endif

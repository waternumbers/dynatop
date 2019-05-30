// Matrix utilities utilising Rcpp

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mmult(NumericMatrix a, NumericMatrix b)
{
 int ncol_a = a.ncol();
 int ncol_b = b.ncol();
 int nrow_b = b.nrow();
 int nrow_a = a.nrow();
 if(ncol_a != nrow_b)
 {
   Rcout << "LH matrix is " << nrow_a << " x " << ncol_a << ", RH matrix is " << nrow_b << " x " << ncol_b;
   stop("incompatible matrix dimensions");
 }
 NumericMatrix c = NumericMatrix(a.nrow(), b.ncol());

 for(int i=0; i < nrow_a; i++)
 {
   for(int j=0; j < b.ncol(); j++)
   {
     c(i,j) = 0;
     for(int k=0; k < nrow_b; k++)
     {
       c(i,j) = c(i,j) + b(k, j)*a(i,k);
     }
   }
 }

 return c;
}

// Matrix x vector
// [[Rcpp::export]]
NumericVector matrix_vect_mult_cpp(NumericMatrix A, NumericVector b)
{
  // Vector multipled by matrix
  if(A.ncol() != b.size())
  {
    Rcout << "Vector size is " << b.size();
    stop("incompatible matrix dimensions");
  }
  NumericVector c = NumericVector(A.nrow());

  for(int i=0; i < A.nrow(); i++)
  {
    c(i) = 0;
    for(int j=0; j < b.size(); j++)
    {
      c(i) = c(i) + b(j)*A(i,j);
    }
  }

  return c;
}

// test them
// vect <- 1:1000
// A <- matrix(ncol=1000, nrow=1000)
// A[1:1000,]<- 1:1000
// mult_matrix_vect_cpp(A, vect)
// mult_vect_matrix_cpp(vect, mt)

// vector x matrix
// [[Rcpp::export]]
NumericVector vect_matrix_mult_cpp(NumericVector b, NumericMatrix A)
{
// returns b %*% A (b implicitly is transposed )
  if(A.nrow() != b.size())
  {
    stop("incompatible matrix dimensions");
  }
  NumericVector c = NumericVector(A.ncol());

  for(int j=0; j < A.ncol(); j++)
  {
    c(j) = 0;
    for(int i=0; i < b.size(); i++)
    {
      c(j) = c(j) + b(i)*A(i,j);
    }
  }

  return c;
}


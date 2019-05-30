#include <Rcpp.h>
// bring the useful matrix utilities
#include "matrix_utils.h"

using namespace Rcpp;

// Unsaturated drainage as function of storage deficit and unsaturated zone moisture content
// [[Rcpp::export]]
NumericVector funcpp_uz(NumericVector suz, NumericVector sd, NumericVector td, double dt=0.0001)
{
  NumericVector uz = NumericVector(suz.size());
// limit the amount that can be drained over the given time step (defaults to very lareg value)
  NumericVector uz_max = suz/dt;
  for(int i=0; i < suz.size(); i++)
  {
    if(sd(i)>0 && suz(i)>0)
    {
      // only drainage if there's spare moisture content
      // unsaturated drainage is ratio of unsaturated zone storage to defict / residence time td
      uz(i) = std::min(suz(i)/(sd(i)*td(i)), uz_max(i));
    }
    else
    {
      // prevent divide by zero - no drainage when saturated
      uz(i) = 0;
    }
  }
  return uz;
}


// Derivative function for ODE solver.
// Rate of change of subsurface discharges
// [[Rcpp::export]]
List funcpp_dqdt(NumericVector t, 
              NumericVector y, 
              List parms) 
{
  // unpack paramters
  NumericVector m = as<NumericVector>(parms["m"]);
  NumericMatrix Wdash = as<NumericMatrix>(parms["Wdash"]);
  NumericVector uz = as<NumericVector>(parms["uz"]);
  NumericVector qbf_max = as<NumericVector>(parms["qbf_max"]);

  NumericVector res = (y/m) + matrix_vect_mult_cpp(Wdash,y) + uz;
  
  for(int i=0; i < res.length(); i++){
    if(res[i]>qbf_max[i] && res[i]>0){
      res[i] = 0;
    }
  }
					    
  return as<List>(res);
}

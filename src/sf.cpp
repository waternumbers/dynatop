#include "sf.h"

// solve 
sfc::sfc(){ }
double sfc::fq(double &s, double &qin){return(-999.9);}

// constant celerity and diffusivity
sfc_cnstCD::sfc_cnstCD(std::vector<double> const &param, double const &Dx){
  kappa = param[0]/Dx; // celerity divided by length to get q from storage
  eta = 0.5 - (param[1] / (param[0]*Dx)); // could retrun negative eta....
}
double sfc_cnstCD::fq(double &s, double &qin){
  //Rcpp::Rcout << "in flow call " << kappa << " " << eta << " " << s << " " << qin << std::endl;
  return( std::max(0.0, ( (kappa*s) - (eta*qin)) / (1-eta) ) );
}

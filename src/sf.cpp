#include "sf.h"

// solve 
sfc::sfc(){ eta = -999.9; }
double sfc::fq(double &s, double &ain){ return(-999.9); }
double sfc::fa(double &q){ return(-999.9); }

// constant celerity
sfc_cnstC::sfc_cnstC(std::vector<double> const &param, double const &Dx_):
  Dx(Dx_), celerity(param[0])
{
  //Rcpp::Rcout << "sf celerity param[0] " << param[0] << std::endl;
  //Rcpp::Rcout << "sf celerity " << celerity << std::endl;
  eta = 0.5;
}
double sfc_cnstC::fq(double &s, double &ain){
  //Rcpp::Rcout << "sf celerity " << celerity << std::endl;
  double a = std::max(0.0, (s - (Dx*eta*ain))/(Dx*(1-eta))); // outflow area
  //Rcpp::Rcout << "sf computing q " << a << " " << Dx << " " << eta << " " << ain << " " << celerity << std::endl;
  return( a*celerity ); }
double sfc_cnstC::fa(double &q){
  //Rcpp::Rcout << "sf celerity " << celerity << std::endl;
  //Rcpp::Rcout << "sf computing a " << q << " " << celerity << std::endl;
  return( q/celerity ); }

// constant celerity and difusivity
sfc_cnstCD::sfc_cnstCD(std::vector<double> const &param, double const &Dx_):
  Dx(Dx_), celerity(param[0])
{
  double D = param[1];
  eta = std::max(0.0, 0.5 - D/(celerity*Dx));
}
double sfc_cnstCD::fq(double &s, double &ain){
  double a = std::max(0.0, (s - (Dx*eta*ain))/(Dx*(1-eta))); // outflow area
  return( a*celerity );
}
double sfc_cnstCD::fa(double &q){ return( q/celerity ); }

#include "sf.h"

// solve 
sfc::sfc(){};
double sfc::fq(double &s, double &qin){ return(-999.9); }
double sfc::finit(double &q, double &qin){ return(-999.9); }

// constant celerity
sfc_cnstC::sfc_cnstC(std::vector<double> const &param, double const &A, double const &w){
  Dx = A/w;
  celerity = param[0];
  eta = 0.5;
};
double sfc_cnstC::fq(double &s, double &qin){
  double a = qin/celerity; // inflow area
  a = std::max(0.0, (s - (Dx*eta*a))/(Dx*(1-eta))); // outflow area
  return( a*celerity );
}
double sfc_cnstC::finit(double &q, double &qin){
  if( q==0.0 ){  return(0.0); }
  return( (Dx/(2.0*celerity))*(q+qin) );
  //return( (1.0/(2.0*celerity))*(q+qin) );
}

// constant celerity and difusivity
sfc_cnstCD::sfc_cnstCD(std::vector<double> const &param, double const &A, double const &w){
  Dx = A/w;
  celerity = param[0];
  double D = param[1];
  eta = std::max(0.0, 0.5 - D/(celerity*Dx));
}
double sfc_cnstCD::fq(double &s, double &qin){
  double a = qin/celerity; // inflow area
  a = std::max(0.0, (s - (Dx*eta*a))/(Dx*(1-eta))); // outflow area
  return( a*celerity );
}
double sfc_cnstCD::finit(double &q, double &qin){
  if( q==0.0 ){  return(0.0); }
  return( (Dx/celerity)*( eta*qin + (1-eta)*q ) );
}

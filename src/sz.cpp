#include "sz.h"

szc::szc(){
  Dx = -999.9;
  D = -999.9;
  q_szmax = -999.9;
};

void szc::fke(double &kappa, double &eta, double const &s){};
double szc::fc(double const &s){ return(-9999.9); };


// bounded exponential
szc_bexp::szc_bexp(std::vector<double> const &param, double const &grd, double const &A, double const &w){
  szc();
  const double &t_0(param[0]), &m(param[1]);
  D = param[2];
  Dx = A/w;
  // Rcpp::Rcout << "D is " << D << std::endl;
  double beta = std::atan(grd);
  psi = std::cos(beta) / m ;
  omega = t_0*std::sin(beta);
  q_szmax = ( w * omega * ( 1 -  std::exp(-psi*D) ) ) / A ;
};
void szc_bexp::fke(double &kappa, double &eta, double const &s){
  kappa  = Dx / fc(s);
  eta = 0.5;
};
double szc_bexp::fc(double const& s){
  if(s > D){ return( 0.0 ); };
  return( psi * omega *  std::exp(-s*psi)  );
};


// constant celerity/velocity
szc_cnst::szc_cnst(std::vector<double> const &param, double const &A, double const &w){
  szc();
  //const double &vsz(param[0]), &maxH(param[1]);
  celerity = param[0];
  D = param[1];
  kappa0 = A / (w*param[0]);
  q_szmax = (w*D*param[0]) / A;
};
void szc_cnst::fke(double &kappa, double &eta, double const &s){
  kappa = kappa0;
  eta = 0.5;
};
double szc_cnst::fc(double const& s){
  if(s > D){ return( 0.0 ); };
  return( celerity );
};

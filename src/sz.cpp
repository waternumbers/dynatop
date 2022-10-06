#include "sz.h"

szc::szc(){}
double szc::fq(double &s, double &qin){ return(-999.9); } // flux for crossectional area

// exponential bounded below
szc_exp::szc_exp(std::vector<double> const param, double const grd,
		 double const A, double const w, double const Dx){
  szc();
  double const &D(param[0]), t0(param[1]), &m(param[2]);
  double beta = std::atan(grd);
  
  sz_max = D*A; // maximum storage D time area
  kappa =  w*t0*std::sin(beta);
  psi = std::cos(beta) / (m*w*Dx); // scaling to get crosssectional depth from storage
  eta = std::exp( - std::cos(beta)*D/m ); // offset due to minimum depth
  //Rcpp::Rcout << D << " " << t0 << " " << m << " " << grd << " " << beta << std::endl;
  //Rcpp::Rcout << kappa << " " << psi << " " << eta << std::endl;
  q_szmax = kappa * (1.0 - eta);
}
double szc_exp::fq(double &s, double &qin){
  double q = kappa * ( std::exp(-psi*s) - eta );
  q = std::min( q_szmax, std::max( 2*q - qin, 0.0)); // bracket for unreasonable values of s
  return( q );
}


  
// // bounded exponential
// szc_bexp::szc_bexp(std::vector<double> const &param, double const &grd, double const &A, double const &w)
// {
//   szc();
//   double const &t_0(param[0]), &m(param[1]), &D(param[2]);
//   Dx = A/w;
//   // Rcpp::Rcout << "D is " << D << std::endl;
//   double beta = std::atan(grd);
//   psi = std::cos(beta) / m ;
//   omega = w*t_0*std::sin(beta)/A;
//   kappa = std::exp(-psi*D);
//   q_szmax = omega * ( 1 -  kappa );
// };
// double szc_bexp::fa(double &q){
//   return( -std::log((q/omega)+kappa)/psi );
// };
// // constant celerity/velocity
// szc_cnst::szc_cnst(std::vector<double> const &param, double const &A, double const &w){
//   szc();
//   //const double &vsz(param[0]), &maxH(param[1]);
//   celerity = param[0];
//   D = param[1];
//   kappa0 = A / (w*param[0]);
//   q_szmax = (w*D*param[0]) / A;
// };
// void szc_cnst::fke(double &kappa, double &eta, double const &s){
//   kappa = kappa0;
//   eta = 0.5;
// };
// double szc_cnst::fc(double const& s){
//   if(s > D){ return( 0.0 ); };
//   return( celerity );
// };

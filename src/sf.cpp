#include "sf.h"

// solve 
sfc::sfc(){
  kappa0 = -999.9;
  eta0 = -999.9;
   };

// the following are valid for constant kappa and eta
void sfc::fke(double &kappa, double &eta, double const &s){
  // compute 1-eta allowing for positive outflow condition
  kappa = kappa0;
  eta = eta0;
};

// constant celerity
sfc_cnstC::sfc_cnstC(std::vector<double> const &param, double const &A, double const &w){
  // param[0] is celerity
  eta0 = 0.5;
  kappa0 = A / (w* param[0]); // Dx / celerity
};
// constant celerity
sfc_cnstCD::sfc_cnstCD(std::vector<double> const &param, double const &A, double const &w){
  // param[0] is celerity
  double Dx = A/w;
  kappa0 = Dx / param[0]; // Dx / celerity
  eta0 = 0.5 - ( param[1] / (param[0]*Dx) );
};

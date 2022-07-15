#include "sf.h"

sfc::sfc(){};
void sfc::init(double &s, double &q, double &qin, double &r){ };
double sfc::max_down(double &s, double &qin, double const &Dt){ return (-9999.9); };
void sfc::solve(double &s, double &q, double &qin, double &r, double const &Dt){ };

// constant celerity
sfc_cnst::sfc_cnst(std::vector<double> const &param, double const &A, double const &w){
  // param[0] is celerity
  kappa = A / (w * param[0]); // Dx / celerity
};

void sfc_cnst::init(double &s, double &q, double &qin, double &r){
  // revise s and q on basis of known inflwo qin and outflow r assuming steady state
  q = qin - r;
  s = kappa*q;
};

double sfc_cnst::max_down(double &s, double &qin, double const &Dt){
  eta=0.5; // since kinematic
  if( Dt < kappa*eta ){ eta = Dt/kappa; }; // introduce dispersion to stop negative flows...
  return( (s + qin*(Dt-(kappa*eta)))/Dt );
};

void sfc_cnst::solve(double &s, double &q, double &qin, double &r, double const &Dt){
  eta=0.5; // since kinematic
  if( Dt < kappa*eta ){ eta = Dt/kappa; };
  
  // revise s and q on basis of known inflow qin and outflow r
  double z = ( (kappa*(1-eta)) / (kappa*(1-eta) + Dt) ) * (s + Dt * ( (qin/(1-eta)) - r) );
  q = (s-z)/Dt + qin - r ;
  s = z;

};

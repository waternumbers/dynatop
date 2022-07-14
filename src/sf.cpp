#include "sf.h"

sfc::sfc(){};
void sfc::init(double &s, double &q, double &qin, double &r){ };
double sfc::max_down(double &s, double &qin, double const &Dt){ return (-9999.9); };
void sfc::solve(double &s, double &q, double &qin, double &r, double const &Dt){ };

// constant celerity
sfc_cnst::sfc_cnst(std::vector<double> const &param, double const &A, double const &w){
  celerity = param[0];
  Dx = A/w;
};

void sfc_cnst::init(double &s, double &q, double &qin, double &r){
  // revise s and q on basis of known inflwo qin and outflow r assuming steady state
  q = qin - r;
  s = Dx*q/celerity;
};

double sfc_cnst::max_down(double &s, double &qin, double const &Dt){ return( (s/Dt) + qin ); };

void sfc_cnst::solve(double &s, double &q, double &qin, double &r, double const &Dt){
  // revise s and q on basis of known inflow qin and outflow r
  double z = ( Dx / (Dx + (2*Dt*celerity)) )* (s + Dt * ( 2*qin - r) );
  q = (s -z)/Dt + qin - r ;
  s = z;
};

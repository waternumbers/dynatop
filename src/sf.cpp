#include "sf.h"

// solve 
sfc::sfc(){ }
// fq compute the multiplier of storage that gives outflow
double sfc::fk(double const &s){
  double vbar = (v_1 * std::min(s,s_1) + v_2 * std::max(0.0,s-s_1))/s;
  return( vbar*rho );
}

// two section constant velocity
sfc_cnst::sfc_cnst(std::vector<double> const &param, double const &width, double const &area){
  v_1 = param[0]; // lower section
  v_2 = param[2]; // upper section
  s_1 = param[1]; // threshold volumne
  rho = width/area;
}


// Kinematic with raf
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, double const &width, double const &grad, double const &area){
  v_1 = param[0]; // lower section
  s_1 = param[1]; // threshold volume
  v_2 = std::pow(grad,1.0/2.0) / param[2]; // gradient and mannings n part of velocity equation
  rho = width/area;
}
// based on wide channel approximation
double sfc_kin::fk(double const &s){
  double ss = std::max(0.0,s-s_1);
  double vv = v_2 * std::pow(rho*ss,2.0/3.0);
  double vbar = (v_1 * std::min(s,s_1) + vv * ss)/s;
  return( vbar*rho );
}

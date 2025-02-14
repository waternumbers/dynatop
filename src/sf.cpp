#include "sf.h"

// solve 
sfc::sfc(){ }
void szc::update(double const &q){}


// two section constant velocity
sfc_cnst::sfc_cnst(std::vector<double> const &param, double const &_Dx){
  v_1 = param[0]; // lower section
  v_2 = param[2]; // upper section
  q_1 = param[0] * param[1] / _Dx; // threshold volume
  Dx = _Dx;
  eta = 0.0;
}
// fq compute the multiplier of storage that gives outflow
double sfc::fk(double const &q){
  double a = ( std::min(q,q_1) /v_1 ) + std::max(0.0, q-q_1)/v_2 ;
  kappa = Dx *a /q;
  return( vbar*rho );
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

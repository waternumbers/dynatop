#include "sf.h"

// solve 
sfc::sfc(){ }
void szc::update(double const &q){}


// two section constant velocity
sfc_cnst::sfc_cnst(std::vector<double> const &param, double const &_Dx):
  v_1(param[0]), v_2(param[2]), q_1(param[0] * param[1] / _Dx),
  Dx(_Dx), eta(0.0)
{}

//   v_1 = param[0]; // velocity in upper section
//   v_2 = param[2]; // upper section
//   q_1 = param[0] * param[1] / _Dx; // threshold volume
//   Dx = _Dx;
//   eta = 0.0;
// }
// fq compute the multiplier of storage that gives outflow
double sfc::update(double const &q){
  double a = ( std::min(q,q_1) /v_1 ) + std::max(0.0, q-q_1)/v_2 ;
  kappa = Dx*a/q;
  return( vbar*rho );
}


// Kinematic with raf
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, double const &_width, double const &grad, double const &_Dx){
  width = _width;
  Dx = _Dx;
  rho = param[0] / (width * std::pow(grad,1.0/2.0)); // needed in computing h
  v_1 = _Dx / param[2]; // lower section velocity
  q_1 = v_1 * param[1]/_Dx; // threshold flow
  eta = 0.5;
}
// based on wide channel approximation
double sfc_kin::update(double const &q){
  double a_2 = width * std::power( std::max(0.0, q-q_1) * rho, 3.0/5.0 );
  kappa = Dx*( (std::min(q,q_1) /v_1) + a_2 )/q;
}

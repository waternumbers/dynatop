#include "sf.h"

// solve 
sfc::sfc(){ }
void sfc::update(double const &q){}

/// <TODO> check order fo properties
// two section constant velocity
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &prop){
  v_1=param[0];
  v_2=param[2];
  q_1=param[0] * param[1] / prop[2];
  Dx=prop[2];
  eta=0.0;
  kappa=0.0;
};

void sfc_cnst::update(double const &q){
  double a = ( std::min(q,q_1) /v_1 ) + std::max(0.0, q-q_1)/v_2 ;
  kappa = Dx*a/q;
}


// Kinematic with raf
// TODO check order of properties
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &prop){ //double const &_width, double const &grad, double const &_Dx){
  width = prop[1];
  Dx = prop[2];
  rho = param[0] / (width * std::pow(prop[3],1.0/2.0)); // needed in computing h
  v_1 = Dx / param[2]; // lower section velocity
  q_1 = v_1 * param[1]/Dx; // threshold flow
  eta = 0.5;
}
// based on wide channel approximation
void sfc_kin::update(double const &q){
  double a_2 = width * std::pow( std::max(0.0, q-q_1) * rho, 3.0/5.0 );
  kappa = Dx*( (std::min(q,q_1) /v_1) + a_2 )/q;
}

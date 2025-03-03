#include "sf.h"

// solve 
sfc::sfc(){ }
void sfc::update(double const &q){}

/// <TODO> check order fo properties
// two section constant velocity
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &prop){
  // param is v_sf, s_raf, t_raf
  Dx = prop[1];
  v_raf = prop[1] / param[2];
  q_raf = param[1] / param[2];
  v_sf = param[0];
  eta=0.0;
};

void sfc_cnst::update(double const &q){
  if( q==0.0 ){
    kappa = 0.0;
  }else{
    double a = ( std::min(q,q_raf) /v_raf ) + std::max(0.0, q-q_raf)/v_sf ;
    kappa = Dx*a/q;
  }
}

// Kinematic with raf
// TODO check order of properties
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &prop){
  // param is n, s_raf, t_raf
  Dx = prop[1];
  v_raf = prop[1] / param[2];
  q_raf = param[1] / param[2];
  
  width = prop[0]/Dx;
  rho = param[0] / (width * std::pow(prop[2],1.0/2.0)); // needed in computing h
  eta = 0; //.5;
}
// based on wide channel approximation
// TODO <check>
void sfc_kin::update(double const &q){
  if( q==0 ){
    kappa = 0.0;
  }else{
    double a_2 = width * std::pow( std::max(0.0, q-q_raf) * rho, 3.0/5.0 );
    kappa = Dx*( (std::min(q,q_raf) /v_raf) + a_2 )/q;
  }
}

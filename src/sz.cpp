#include "sz.h"

szc::szc(){}
double szc::fq(double &s){ return(-999.9); } // flux for crossectional area
double szc::fs(double &q){ return(-999.9); } // flux for crossectional area

// exponential
szc_exp::szc_exp(std::vector<double> const &param, std::vector<double> const &prop){
  szc();
  
  double const &t0(param[0]), &m(param[1]);
  double const &area(prop[0]), &width(prop[1]), &grd(prop[3]);

  double beta = std::atan(grd);
  q_szmax =  width*t0*std::sin(beta);
  psi = std::cos(beta) / (m*area); // scaling to get crosssectional depth from storage
}
double szc_exp::fq(double &s){ // get flow from storage
  double q = q_szmax * std::exp(-psi*s);
  return( q );
}
double szc_exp::fs(double &q){ // get storage from flow
  double s = -std::log(q/q_szmax) / psi;
  return( s );
}

  
// bounded exponential
szc_bexp::szc_bexp(std::vector<double> const &param, std::vector<double> const &prop){
  szc();
  double const &t_0(param[0]), &m(param[1]), &h_sz_max(param[2]);
  double const &area(prop[0]), &width(prop[1]), &grd(prop[3]);
  double beta = std::atan(grd);
  
  psi = std::cos(beta) / (m*area) ;
  omega = width*t_0*std::sin(beta);
  kappa = std::exp(-psi*h_sz_max);
  q_szmax = omega * ( 1 -  kappa );
}
double szc_bexp::fq(double &s){ // get flow from storage
  double q = std::max(0.0, omega*( std::exp(-psi*s) - kappa ) );
  return( q );
}
double szc_bexp::fs(double &q){ // get storage from flow
  return( -std::log((q/omega)+kappa)/psi );
};


// constant celerity/velocity
szc_cnst::szc_cnst(std::vector<double> const &param,  std::vector<double> const &prop){
  szc();
  //const double &vsz(param[0]), &maxH(param[1]);
  double const &v_sz(param[0]), &h_sz_max(param[1]);
  double const &area(prop[0]), &width(prop[1]);
  
  omega = width*v_sz;
  psi=1.0/area;
  kappa = h_sz_max;
  q_szmax = omega*h_sz_max;
};
double szc_cnst::fq(double &s){
  return( std::max(0.0, omega*(kappa - (s*psi))) );
};
double szc_cnst::fs(double &q){
  return( -psi*((q/omega)-kappa) );
};

// double exponential
szc_dexp::szc_dexp(std::vector<double> const &param, std::vector<double> const &prop){
  szc();
  
  double const &t0(param[0]), &m(param[1]), &m2(param[2]);
  double const &area(prop[0]), &width(prop[1]), &grd(prop[3]);

  double beta = std::atan(grd);
  omega = param[3]; // weight
  q_szmax =  width*t0*std::sin(beta);
  psi = std::cos(beta) / (m*area); // scaling to get crosssectional depth from storage
  kappa = std::cos(beta) / (m2*area); // scaling to get crosssectional depth from storage
}
double szc_dexp::fq(double &s){ // get flow from storage
  double q = q_szmax * ( omega*std::exp(-psi*s) + (1.0-omega)*std::exp(-kappa*s) );
  return( q );
}
double szc_dexp::fs(double &q){ // get storage from flow
  double lwr = -std::log(q/q_szmax) / psi;
  double upr = -std::log(q/q_szmax) / kappa;
  if(upr < lwr){
    double tmp(upr);
    upr=lwr;
    lwr=tmp;
  }
  //bisection to find solution
  int it(0), max_it(1000);
  double z, qq;
  while( (it <= max_it) and ( (upr-lwr)>1e-10 ) ){
    z = (lwr+upr)/2.0;
    qq = fq(z);
    if( qq <= q ){ upr = z; } else { lwr = z; }
  }
  z = (lwr+upr)/2.0;
  return( z );
}

#include "sf.h"

// solve 
sfc::sfc(){ }
// fq compute the outflow given the other variables
double sfc::fq(double const &s){ //, double const &qin, double const &r){
  double qq = std::max(0.0, std::min(s,s_1)*kappa_1 )
    + std::max(0.0, std::max(s-s_1,0.0)*kappa_2) ;
  return( qq );
}
// fs computes  steady state storage given the inflows
double sfc::fs(double const &q){ //in, double const &r){
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  double q_1_max = std::max(0.0, s_1*kappa_1);
  double q_1 = std::min(q,q_1_max);
  double q_2 = std::max(0.0,q-q_1_max);
  double s = std::max(0.0,q_1*kappa_1) + std::max(0.0,q_2/kappa_2);
  return(s);
}
void sfc::update(double &s, double &q, double const &qin, double const &vin,
		 double const &Dt, double const &vtol, int const &max_it){
  
  double sfmax = s + Dt*qin - vin;
  //double rin = vin / Dt;
  std::pair<double,double> lbnd(0.0, 999.9);
  double qq = fq(lbnd.first); //,qin,rin);
  lbnd.second = sfmax - Dt*qq - lbnd.first;
  
  std::pair<double,double> ubnd(sfmax, 999.9);
  qq = fq(ubnd.first); //,qin,rin);
  ubnd.second = sfmax - Dt*qq - ubnd.first;
  
  int it = 0;
  while( (it <= max_it) and ( lbnd.second > vtol ) ){ //( (bnd.second - bnd.first)>vtol ) ){
    double iW = ubnd.second / (ubnd.second-lbnd.second);
    iW = std::max(0.001,std::min(iW,0.999));
    double z = (iW*lbnd.first) + (1.0-iW)*ubnd.first;
    qq = fq(z); //,qin,rin);
    double Sw = sfmax - Dt*qq - z;
    if( Sw <= 0 ){ //bnd.second= z; } else { bnd.first=z; }
      ubnd.first = z;
      ubnd.second = Sw;
    }else{
      lbnd.first = z;
      lbnd.second = Sw;
    }
    it += 1;
  }
  double z = lbnd.first;
  q = qin + (s - vin - z)/Dt;
  s = z;
};

void sfc::iter_update(double &s, double &q, double const &qin, double const &vin,
		      double const &Dt, double const &vtol, int const &max_it){
  
  double sfmax = s + Dt*qin - vin;
  double z = s;
  double v(0.0);
  int it(0);
  while( it<= 3 ){ //max_it ){
    if(z==0.0){ v = 0.0; }
    else{ v = fq(z)/z; }
    z = sfmax / ( 1.0 + Dt*v);
    it += 1;
  }
  s = z;
  q = (sfmax - s)/Dt;
}
    

// constant celerity, diffusivity with raf
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  // lower path is linear tank
  kappa_1 = 1.0/param[2]; // param[2] is t_raf
  s_1 = param[1]; // param[2] is raf storage
  kappa_2 = param[0]/ properties[1]; // celerity divided by length to get q from storage
}

// Kinematic with raf
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const &Dx(properties[1]), &area(properties[0]), &grd(properties[2]);
  double const width = area/Dx;
  double const &n(param[0]);
  kappa_1 = 1.0 / param[2]; // param[2] is raf time constant
  eta_1 = 0.0;
  s_1 = param[1]; // raf storage
  kappa_2 = Dx;
  eta_2 = std::pow(grd,0.5) / (n * std::pow(width,(2.0/3.0)));
}
double sfc_kin::fq(double const &s){ //, double const &qin, double const &r){
  double q_1 = kappa_1*std::min(s_1,s);
  double q_2 = eta_2 * std::pow( std::max(0.0,(s-s_1)/kappa_2), (5.0/3.0) );
  return( q_1 + q_2 );
}
double sfc_kin::fs(double const &qin){ //, double const &r){
  double q = qin;
  if( q<= 0.0 ){ return(0.0); } // handl case of no outflow

  double q_1_max = s_1*kappa_1;
  double q_1 = std::min(q,q_1_max);
  double s_1 = q_1 / kappa_1;
  double q_2 = std::max(0.0,q-q_1_max);
  double s_2 = kappa_2 * std::pow( q_2/eta_2, 3.0/5.0 );
  return( s_1 + s_2 );
}

// compound channel
sfc_comp::sfc_comp(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[1]);
  kappa_1 = param[0]/Dx; // velocity divided by length to get q from storage for lower part of channel
  s_1 = param[1]; // max stoage in lower part of channel
  kappa_2 = param[2]/Dx; // velocity divided by length to get q from storage for upper part of channel
}

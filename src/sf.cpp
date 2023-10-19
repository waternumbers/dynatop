#include "sf.h"

// solve 
sfc::sfc(){ }
double sfc::fq(double &s, double &qin){return(-999.9);}
double sfc::fs(double &q, double &qin){return(-999.9);}


// constant celerity, diffusivity with raf
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[2]);
  kappa = param[0]/Dx; // celerity divided by length to get q from storage
  eta = 0.5 - (param[1] / (param[0]*Dx)); // could retrun negative eta....
  s_raf = param[2]; // raf storage
  t_raf = param[3]; // raf time constant
  q_raf = s_raf/t_raf; // raf max flow
}
double sfc_cnst::fq(double &s, double &qin){
  double qq(-999.9);
  if( s <= s_raf ){
    qq = s / t_raf;
  }else{
    qq = ( kappa*(s-s_raf) - (eta*std::max(qin-q_raf,0.0)) ) / (1-eta); // flow from outside raf
    qq = std::max(0.0,qq); // ensure positive
    qq += q_raf; // add raf outfloe
  }
  return( qq );
}
double sfc_cnst::fs(double &q, double &qin){
  double ss(-999.9);
  if( q > q_raf ){
    ss = ( (1-eta)*(q-q_raf) + eta*(std::max(0.0,qin-q_raf)) ) / kappa;
    ss += s_raf;
  }else{
    ss = q*t_raf;
  }
  return(ss);
}

// Kinematic with raf
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const &Dx(properties[2]), &width(properties[1]), &grd(properties[3]);
  double const &n(param[0]);
  kappa = Dx;
  omega = std::pow(grd,0.5) / (n * std::pow(width,(2.0/3.0)));
  s_raf = param[1]; // raf storage
  t_raf = param[2]; // raf time constant
  q_raf = s_raf/t_raf; // raf max flow
  
}
double sfc_kin::fq(double &s, double &qin){
  double qq(-999.9);
  if( s <= s_raf ){
    qq = s / t_raf;
  }else{
    qq = omega * std::pow( (s-s_raf)/kappa, (5.0/3.0) );
    qq = std::max( 0.0, 2*qq - std::max(qin-q_raf,0.0) );
    qq += q_raf; // add raf outfloe
  }
  return( qq );
}

double sfc_kin::fs(double &q, double &qin){
  double ss(-999.9);
  if( q > q_raf ){
    ss = kappa * std::pow( (q+qin)/(2*omega), 3.0/5.0 ) ;
    ss += s_raf;
  }else{
    ss = q*t_raf;
  }
  return( ss );
}


// compound channel
sfc_comp::sfc_comp(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[2]);
  kappa_1 = param[0]/Dx; // velocity divided by length to get q from storage for lower part of channel
  eta_1 = 0.5 - (param[1] / (param[0]*Dx));  // Dispersion property for lower part of channel
  s_1 = param[2]; // max stoage in lower part of channel
  kappa_2 = param[3]/Dx; // velocity divided by length to get q from storage for upper part of channel
  eta_2 = 0.5 - (param[4] / (param[3]*Dx));  // Dispersion property for upper part of channel
  q_1_max = s_1*kappa_1; // max inflow to lower part of channel
 
  //Rcpp::Rcout << eta_1 << " " << kappa_1 << " " << s_1 << " " << q_1_max << std::endl;
  //Rcpp::Rcout << eta_2 << " " << kappa_2 << std::endl;
  
}
double sfc_comp::fq(double &s, double &qin){
  double q_1 = std::min(qin,q_1_max);
  double q_2 = std::max(0.0,qin-q_1_max);
  
  double qq = std::max(0.0, (std::min(s,s_1)*kappa_1 - eta_1*q_1) / (1-eta_1) )
    + std::max(0.0, (std::max(s-s_1,0.0)*kappa_2 - eta_2*q_2) / (1-eta_2)) ;
  return( qq );
 
}
double sfc_comp::fs(double &q, double &qin){
  // bool flg(false);
  //   if( (q==0.0) and (qin > 0.0) ){
  //   flg = true;
  // }
  
  double q_1 = std::min(qin,q_1_max) ; //std::max(0.0,std::min(qin,q_1_max));
  double q_2 = std::max(0.0,qin-q_1_max);
  double qq = std::max(0.0, (q_1_max - eta_1*q_1)/(1-eta_1)); // outflow at stroage at s_1
  // if( flg ){
  //   Rcpp::Rcout << q << " " << qin << " " << qq << " " << q_1 << " " << q_2 << " " << s_1 << std::endl;
  // }
  double ss(-999.9);
  if( qq >= q ){ // then storage between 0 & s_1
    //    if( flg ){ Rcpp::Rcout <<  "In stage 1" << " " << eta_1 << std::endl; }
    ss = ( (1-eta_1)*q + eta_1*q_1 ) / kappa_1;
  }else{
    qq = q-qq; // flow from upper part of channel
    // if( flg ){ Rcpp::Rcout <<  "In stage 2" << " " << qq << std::endl; }
    ss = ( (1-eta_2)*qq + eta_2*q_2 ) / kappa_2;
    ss += s_1;

  }
  return(ss);
  
}

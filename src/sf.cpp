#include "sf.h"

// solve 
sfc::sfc(){ }
// fq compute the outflow given the other variables
double sfc::fq(double const &s, double const &qin, double const &r){
  double q_1_max = std::max(0.0, s_1*kappa_1 + (1-eta_1)*r);
  if( s_1*kappa_1 == 0.0 ){ q_1_max = 0.0; }

  double q_1 = std::min(qin,q_1_max);
  double q_2 = std::max(0.0,qin-q_1_max);
  
  double qq = std::max(0.0, (std::min(s,s_1)*kappa_1 - eta_1*q_1) / (1-eta_1) )
    + std::max(0.0, (std::max(s-s_1,0.0)*kappa_2 - eta_2*q_2) / (1-eta_2)) ;
  return( qq );
}
// fs computes  steady state storage given the inflows
double sfc::fs(double const &qin, double const &r){
  double q = qin - r;
  if( q<= 0.0 ){ return(0.0); } // handl case of no outflow

  double q_1_max = std::max(0.0, s_1*kappa_1 + (1-eta_1)*r);
  if( s_1*kappa_1 == 0.0 ){ q_1_max = 0.0; }

  double q_1 = std::min(qin,q_1_max);
  double q_2 = std::max(0.0,qin-q_1_max);
  double qq = std::max(0.0, (s_1*kappa_1 - eta_1*q_1) / (1-eta_1) );// flow at s_1
  double s(-999.9);
  if( qq <= q ){// then in upper part of the storage
    s = s_1 + ( (1-eta_2)*(q-qq) + eta_2*q_2 )/kappa_2;
  }else{ //in lower part of the storage
    s = ( (1-eta_1)*q + eta_1*q_1 )/kappa_1;
  }
  return(s);
}
void sfc::update(double &s, double &q, double const &qin, double const &vin,
		 double const &Dt, double const &vtol, int const &max_it){
  
  double sfmax = s + Dt*qin - vin;
  double rin = vin / Dt;
  std::pair<double,double> lbnd(0.0, 999.9);
  double qq = fq(lbnd.first,qin,rin);
  lbnd.second = sfmax - Dt*qq - lbnd.first;
  
  std::pair<double,double> ubnd(sfmax, 999.9);
  qq = fq(ubnd.first,qin,rin);
  ubnd.second = sfmax - Dt*qq - ubnd.first;
  
  int it = 0;
  while( (it <= max_it) and ( lbnd.second > vtol ) ){ //( (bnd.second - bnd.first)>vtol ) ){
    //z = (bnd.first+bnd.second)/2.0;
    //z = (lbnd.first+ubnd.first)/2.0;
    double iW = ubnd.second / (ubnd.second-lbnd.second);
    iW = std::max(0.001,std::min(iW,0.999));
    double z = (iW*lbnd.first) + (1.0-iW)*ubnd.first;
    qq = fq(z,qin,rin);
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
  //z = bnd.first;
  q = qin + (s - vin - z)/Dt;
  s = z;
};
    

// constant celerity, diffusivity with raf
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[2]);
  // lower path is linear tank
  kappa_1 = 1.0/param[3]; // param[3] is t_raf
  eta_1 = 0.0;
  s_1 = param[2]; // param[2] is raf storage
  kappa_2 = param[0]/Dx; // celerity divided by length to get q from storage
  eta_2 = 0.5 - (param[1] / (param[0]*Dx)); // could retrun negative eta....
  //Rcpp::Rcout << eta_1 << " " << kappa_1 << " " << s_1 << " " << std::endl;
  //Rcpp::Rcout << eta_2 << " " << kappa_2 << std::endl;  
}

// Kinematic with raf
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const &Dx(properties[2]), &width(properties[1]), &grd(properties[3]);
  double const &n(param[0]);
  kappa_1 = 1.0 / param[2]; // param[2] is raf time constant
  eta_1 = 0.0;
  s_1 = param[1]; // raf storage
  kappa_2 = Dx;
  eta_2 = std::pow(grd,0.5) / (n * std::pow(width,(2.0/3.0)));
}
double sfc_kin::fq(double const &s, double const &qin, double const &r){
  double q_1_max = std::max(0.0, s_1*kappa_1 + (1-eta_1)*r);
  double q_1 = std::min(qin,q_1_max);
  double q_2 = std::max(0.0,qin-q_1_max);

  double qc = eta_2 * std::pow( std::max(0.0,(s-s_1)/kappa_2), (5.0/3.0) );
  double qq = std::max(0.0, (std::min(s,s_1)*kappa_1 - eta_1*q_1) / (1-eta_1) )
    + std::max(0.0, (2*qc - q_2) / (1-eta_2)) ;
  return( qq );
}

double sfc_kin::fs(double const &qin, double const &r){
  double q = qin - r;
  if( q<= 0.0 ){ return(0.0); } // handl case of no outflow
  
  double q_1_max = std::max(0.0, s_1*kappa_1 + (1-eta_1)*r);
  double q_1 = std::min(qin,q_1_max);
  double q_2 = std::max(0.0,qin-q_1_max);
  double qq = std::max(0.0, (s_1*kappa_1 - eta_1*q_1) / (1-eta_1) );// flow at s_1
  double s(-999.9);
  if( qq < q ){// then in upper part of the storage
    double qc = 0.5*(q-qq+q_2);
    s = s_1 + kappa_2 * std::pow( qc/eta_2, 3.0/5.0 );
  }else{ //in lower part of the storage
    s = ( (1-eta_1)*q + eta_1*q_1 )/kappa_1;
  }
  return(s);
}


// compound channel
sfc_comp::sfc_comp(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[2]);
  kappa_1 = param[0]/Dx; // velocity divided by length to get q from storage for lower part of channel
  eta_1 = 0.5 - (param[1] / (param[0]*Dx));  // Dispersion property for lower part of channel
  s_1 = param[2]; // max stoage in lower part of channel
  kappa_2 = param[3]/Dx; // velocity divided by length to get q from storage for upper part of channel
  eta_2 = 0.5 - (param[4] / (param[3]*Dx));  // Dispersion property for upper part of channel
 
  //Rcpp::Rcout << eta_1 << " " << kappa_1 << " " << s_1 << " " << std::endl;
  //Rcpp::Rcout << eta_2 << " " << kappa_2 << std::endl;  
}





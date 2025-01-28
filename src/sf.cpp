#include "sf.h"

// solve 
sfc::sfc(){ }
double sfc::fq(double const &s){
  if( s <= 0.0 ){ return(0.0); } // handle case of no storage
  double q = eta_1 * std::min(s,s_1) +
    eta_2 * std::min(s-s_1,0.0);  
  return( q );
}
double sfc::fs(double const &q){
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  double q_1 = fq(s_1); // flow at change storage
  double s = ( std::min(q,q_1)/eta_1 ) +
    ( std::min(q-q_1,0.0)/eta_2 );
  return( s );
}
void sfc::update(double &s, double &q,
		 double const &Dt, double const &vtol, int const &max_it){
  // presume s is at the maximum value (starting volume + inflows when passed in
  double s0 = s;
  
  // lower bound of search - start at 0.0
  std::pair<double,double> lbnd(0.0, 999.9);
  double qq = fq(lbnd.first);
  lbnd.second = s0 - Dt*qq - lbnd.first;
  
  std::pair<double,double> ubnd(s, 999.9);
  qq = fq(ubnd.first);
  ubnd.second = s0 - Dt*qq - ubnd.first;
  
  int it = 0;
  while( (it <= max_it) and ( lbnd.second > vtol ) ){
    double iW = ubnd.second / (ubnd.second-lbnd.second);
    iW = std::max(0.001,std::min(iW,0.999));
    double z = (iW*lbnd.first) + (1.0-iW)*ubnd.first;
    qq = fq(z);
    double Sw = s0 - Dt*qq - z;
    if( Sw <= 0 ){ //bnd.second= z; } else { bnd.first=z; }
      ubnd.first = z;
      ubnd.second = Sw;
    }else{
      lbnd.first = z;
      lbnd.second = Sw;
    }
    it += 1;
  }
  s = lbnd.first;
  q = (s0-s)/Dt;
  //z = bnd.first;
};


// constant celerity in each block - presume flow area of s/Dx
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[2]);
  eta_1 = param[0] / Dx;
  kappa_1 = 0.0; // not used
  s_1 = param[1];
  eta_2 = param[2] / Dx;
  kappa_2 = 0.0; // not used
}

// generic power law with constant parameters
sfc_power_law::sfc_power_law(std::vector<double> const &param, std::vector<double> const &properties){
  //double const &Dx(properties[2]), &width(properties[1]), &grd(properties[3]);
  kappa_1 = param[0];
  eta_1 = param[1];
  s_1 = param[2];
  kappa_2 = param[3];
  eta_2 = param[4];
}
// fq compute the outflow given the storage volume
double sfc_power_law::fq(double const &s){
  if( s <= 0.0 ){ return(0.0); } // handle case of no storage
  double q = eta_1 * std::pow( std::min(s,s_1), kappa_1 ) +
    eta_2 * std::pow( std::min(s-s_1,0.0), kappa_2 );  
  return( q );
}
// fs computes storage given the outflow
double sfc_power_law::fs(double const &q){
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  double q_1 = fq(s_1); // flow at change storage
  double s = 0.0;
  if ( q > q_1 ){
    s += s_1;
    double ss = (q-q_1)/eta_2;
    if(kappa_2 > 0.0){
      ss = std::pow( ss , 1.0 /kappa_2);
    }
    s += ss;
  }else{
    double ss = q/eta_1;
    if(kappa_1 > 0.0 ){
      ss = std::pow( ss, 1.0 /kappa_1);
    }
    s = ss;
  }
  return(s);
}


// Two stage Kinematic with raf
// Assumes shallow water so wetted perimeter ~ height or s/(Dx*w)
// give v = (1/n) * grd^0.5 * s^0.666 / (Dx*w)*0.666
// or q = (1/n) * grd^0.5 * s*(5/3) / ( Dx^(5/3) * w^(2/3) )
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const &Dx(properties[2]), &width(properties[1]), &grd(properties[3]);
  double const &n(param[0]);
  eta_1 = 1.0 / param[2]; // param[2] is raf time constant
  kappa_1 = 0.0; // not used
  s_1 = param[1]; // raf storage
  kappa_2 = 5.0/3.0;
  eta_2 = std::pow(grd,0.5) / (n * std::pow(width,(2.0/3.0)) * std::pow(Dx, 5.0/3.0) );
}
double sfc_kin::fq(double const &s){
  double q = eta_1*std::min(s_1,s);
  q += eta_2 * std::pow( std::max(0.0,(s-s_1)), kappa_2 );
  return( q );
}
double sfc_kin::fs(double const &q){
  if( q<= 0.0 ){ return(0.0); } // handl case of no outflow
  double q_1 = fq(s_1); // flow at change storage
  double s = 0.0;
  if ( q > q_1 ){
    s += s_1;
    double ss = (q-q_1)/eta_2;
    s += std::pow( ss , 3.0 / 5.0 );
  }else{
    s += q/eta_1;
  }
  return(s);
}



#include "sz.h"

szc::szc(){}
double szc::fs(double const &q){ return(9999.9); }// compute storage deficit for a given flow
std::pair<double,double> szc::fq(double const &s){
  std::pair<double,double> out(9999.9,9999.9);
  return(out);
}

// exponential
szc_exp::szc_exp(std::vector<double> const &param, std::vector<double> const &prop){
  szc();
  
  double const &t0(param[0]), &m(param[1]);
  //double const &area(prop[0]), &width(prop[1]), &grd(prop[3]);
  double const &Dx(prop[1]), &grd(prop[2]);
  double area = prop[0];
  double width = area/Dx;
  double beta = std::atan(grd);
  q_szmax =  width*t0*std::sin(beta);
  psi = std::cos(beta) / (m*area); // scaling to get crosssectional depth from storage
}
double szc_exp::fs(double const &q){ // get storage from flow
  if( q_szmax==0.0 ){ return(0.0); } // since there can be no flow or storage
  if( q >= q_szmax ){ return(0.0); } // since saturated
  double s = -std::log( q/q_szmax ) / psi;
  return( s );
}
std::pair<double,double> szc_exp::fq(double const &s){ // get storage from flow
  std::pair<double,double> out(q_szmax,0.0);
  // if( q_szmax<=0.0 ){ return(out); } // since there can be no flow or storage
  if(s>=0){
    out.first = q_szmax * std::exp(-psi*std::max(0.0,s));
    out.second = -psi*out.first;
  }
  return(out) ; //q_szmax * std::exp(-psi*std::max(0.0,s)) );
}
  
// bounded exponential
szc_bexp::szc_bexp(std::vector<double> const &param, std::vector<double> const &prop){
  szc();
  double const &t_0(param[0]), &m(param[1]), &h_sz_max(param[2]);
  double const &width(prop[0]/prop[1]), &grd(prop[2]);
  double area = prop[0];
  double beta = std::atan(grd);
  
  psi = std::cos(beta) / (m*area) ;
  omega = width*t_0*std::sin(beta);
  kappa = std::exp(-psi*h_sz_max);
  q_szmax = omega * ( 1 -  kappa );
}
double szc_bexp::fs(double const &q){ // get storage from flow
  if( omega ==0.0 ){ return( 0.0 ); }  // since there can be no flow or storage
  if( q >= q_szmax ){ return(0.0); } // since saturated
  return( -std::log((q/omega)+kappa)/psi );
};
std::pair<double,double> szc_bexp::fq(double const &s){ // get flow from storage
  std::pair<double,double> out(q_szmax,0.0);
  if(s>0){
    out.first = std::max(0.0, omega*( std::exp(-psi*s) - kappa ) );
    out.second = -psi*out.first;
  }
  return( out );
}

// constant celerity/velocity
szc_cnst::szc_cnst(std::vector<double> const &param,  std::vector<double> const &prop){
  szc();
  //const double &vsz(param[0]), &maxH(param[1]);
  double const &v_sz(param[0]), &h_sz_max(param[1]);
  double const &width(prop[0]/prop[1]);
  double area = prop[0];
  
  omega = width*v_sz;
  psi=1.0/area;
  kappa = h_sz_max;
  q_szmax = omega*h_sz_max;
};
double szc_cnst::fs(double const &q){
  if( q_szmax==0.0 ){ return(0.0); } // since there can be no flow or storage
  return( std::max(0.0,-psi*((q/omega)-kappa)) );
};
std::pair<double,double> szc_cnst::fq(double const &s){
  std::pair<double,double> out(omega*kappa,0.0);
  if( s>0 ){
    out.first = std::max(0.0, omega*(kappa - (s*psi)));
    if( out.first > 0 ){
      out.second = -omega*psi;
    }
  }
  return( out ) ;
};


// double exponential
szc_dexp::szc_dexp(std::vector<double> const &param, std::vector<double> const &prop){
  szc();
  
  double const &t0(param[0]), &m(param[1]), &m2(param[2]);
  double const &area(prop[0]), &Dx(prop[1]), &grd(prop[2]);
  double width = area/Dx;

  double beta = std::atan(grd);
  omega = param[3]; // weight
  q_szmax =  width*t0*std::sin(beta);
  psi = std::cos(beta) / (m*area); // scaling to get crosssectional depth from storage
  kappa = std::cos(beta) / (m2*area); // scaling to get crosssectional depth from storage
}

std::pair<double,double> szc_dexp::fq(double const &s){ // get flow from storage
  std::pair<double,double> out(q_szmax,0.0);
  if( s>0 ){
    out.first = q_szmax * ( omega*std::exp(-psi*s) + (1.0-omega)*std::exp(-kappa*s) );
    out.second = -q_szmax * ( psi*omega*std::exp(-psi*s) + kappa*(1.0-omega)*std::exp(-kappa*s) );
  }
  return(out);
}

double szc_dexp::fs(double const &q){ // get storage from flow
  if( q_szmax==0.0 ){ return(0.0); } // since there can be no flow or storage
  double z;
  if( q > q_szmax ){
    Rcpp::Rcout << "q > qmax " << q << " " << q_szmax << " " << q - q_szmax << std::endl;
    z = 0.0;
  }
  if( q >= q_szmax ){ z = 0.0; }
  else{
    
    double lwr = -std::log(q/q_szmax) / psi;
    double upr = -std::log(q/q_szmax) / kappa;
    if(upr < lwr){
      double tmp(upr);
      upr=lwr;
      lwr=tmp;
    }
    //bisection to find solution
    int it(0), max_it(100);
    double qq; //z, qq;
    while( (it <= max_it) and ( (upr-lwr)>1e-10 ) ){
      z = (lwr+upr)/2.0;
      qq = q_szmax * ( omega*std::exp(-psi*z) + (1.0-omega)*std::exp(-kappa*z) );
      if( qq <= q ){ upr = z; } else { lwr = z; }
      it += 1;
    }
    z = (lwr+upr)/2.0;
    if(it == max_it){ Rcpp::Rcout << "max_it reached " <<lwr << " " << z << " " << upr << std::endl; }
  }
  return( z );
}

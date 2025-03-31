#include "sz.h"

szc::szc(){}
double szc::fs(double const &q){ return(9999.9); }// compute storage deficit for a given flow
double szc::fq(double const &s){ return(9999.9); }

// double szc::update(double &s, double &q, double const &qin, double const &vin,
// 		   double const &Dt, int const &max_it){

//   double s0 = s - (Dt*qin) - vin;
//   q = qin;
//   for(int it=0; it<max_it; ++it){
//     double Qref = (q + qin)/2.0;
//     double Sref = fs( Qref );
//     q = std::min(q_szmax, std::max(0.0, (Sref - s0)/Dt ));
//   }
//   s = std::max(0.0, s0 + Dt*q);
// };

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
  if( q_szmax<=0.0 ){ return(0.0); } // since there can be no flow or storage
  if (q == 0.0){ return( 1e+66 ); }
  //  double s = std::max(0.0,q); // temp covert to computing q
  //s = q_szmax * std::exp(-psi*s);
  double qq = std::min(q,q_szmax);
  double s = -std::log( qq/q_szmax ) / psi;
  return( s );
}
double szc_exp::fq(double const &s){ // get storage from flow
  if( q_szmax<=0.0 ){ return(0.0); } // since there can be no flow or storage
  return( q_szmax * std::exp(-psi*std::max(0.0,s)) );
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
  double qq = std::min(q,q_szmax);
  return( -std::log((qq/omega)+kappa)/psi );
};


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
  double qq = std::min(q,q_szmax);
  return( -psi*((qq/omega)-kappa) );
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
double szc_dexp::fs(double const &q){ // get storage from flow
  if( q_szmax==0.0 ){ return(0.0); } // since there can be no flow or storage
  double z;
  if( q > q_szmax ){
    Rcpp::Rcout << "q > qmax " << q << " " << q_szmax << " " << q - q_szmax << std::endl;
    z = 0.0;
  }
  if( q == q_szmax ){ z = 0.0; }
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
      Rcpp::Rcout << "TO FIX in dexp" << std::endl;
      qq = 0.0; //ftq(z);
      if( qq <= q ){ upr = z; } else { lwr = z; }
      it += 1;
    }
    z = (lwr+upr)/2.0;
    if(it == max_it){ Rcpp::Rcout << "max_it reached " <<lwr << " " << z << " " << upr << std::endl; }
  }
  return( z );
}

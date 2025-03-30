#include "sf.h"


// solve 
sfc::sfc(){ }
void sfc::update( double const &q ){ // set for two partions of constant velocity with storage
  if( q<= 0.0 ){ kappa = 0.0; return; } // handle case of no outflow
  kappa = ( (std::min(q,q_1))*kappa_1 + std::max(0.0,q-q_1)*kappa_2 )/ q;
}
    

// constant celerity, diffusivity with raf
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  // lower path is linear tank
  kappa_1 = param[2]; // param[2] is t_raf
  q_1 = param[1] / kappa_1; // param[1] is raf storage
  kappa_2 = properties[1] / param[0]; // celerity divided by length to get q from storage
  eta = 0.0;
}

// Kinematic with raf
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const &Dx(properties[1]), &area(properties[0]), &grd(properties[2]);
  // double const width = area/Dx;
  double const &n(param[0]);
  kappa_1 = param[2]; // param[2] is raf time constant
  q_1 = param[1] / kappa_1; // param[1] is raf storage
  kappa_2 = area;
  eta_1 = (area * std::sqrt(grd))/(Dx*n);
  eta = 0.5;
}
void sfc_kin::update(double const &q ){
  if( q<= 0.0 ){ kappa = 0.0; return; } // handle case of no outflow
  double h = std::pow( std::max(0.0,q - q_1)/eta_1, 2.0/5.0 );
  kappa = ( std::min(q,q_1)*kappa_1 + kappa_2*h )/q;
}

// compound channel
sfc_comp::sfc_comp(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[1]);
  kappa_1 = Dx / param[0]; // velocity divided by length to get q from storage for lower part of channel
  q_1 = param[1] / kappa_1; // max flow from the lower store
  kappa_2 = Dx /param[2]; // velocity divided by length to get q from storage for upper part of channel
  eta = 0.0;
}

// arbitary area, discharge relationship
sfc_arb_kin::sfc_arb_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[1]);
  unsigned int n = param.size()/2;
  // Rcpp::Rcout << "param is length " << param.size() << std::endl;
  // Rcpp::Rcout << "n is" << n << std::endl;
  for(unsigned int ii = 0; ii<n; ++ii){
    // Rcpp::Rcout << "ii is " << ii << std::endl;
    // Rcpp::Rcout << "storage is " << param[ii]*Dx << std::endl;
    // Rcpp::Rcout << "Flow is " << param[ii+n] << std::endl;
    s_val.push_back( param[ii]*Dx );
    q_val.push_back( param[ii+n] );
  }
  eta = 0.5;
}
void sfc_arb_kin::update(double const &q ){
  if( q<= 0.0 ){ kappa = 0.0; return; } // handle case of no outflow
  unsigned int n = s_val.size();
  unsigned int ii = 1;
  while( (q_val[ii] < q) & (ii < (n-1)) ){
    ii += 1;
  };
  double s = s_val[ii-1] + ( (s_val[ii] - s_val[ii-1])/(q_val[ii]-q_val[ii-1]) )* (q - q_val[ii-1]);
  kappa = s/q;
}

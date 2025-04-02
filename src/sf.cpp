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
  for(unsigned int ii = 0; ii<n; ++ii){
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

// generic power law with constant parameters
sfc_power_law::sfc_power_law(std::vector<double> const &param, std::vector<double> const &properties){
  Dx = properties[1];
  kappa_1 = param[3]; // param[3] is raf time constant
  eta_1 = param[2] / kappa_1; // flow when raf is full param[2] is raf storage
  kappa_2 = param[1]; // power
  eta_2 = param[0]; // scale
  eta = 0;
}
void sfc_power_law::update(double const &q ){
  if( q<= 0.0 ){ kappa = 0.0; return; } // handle case of no outflow
  double s = ( kappa_1 * std::min(q,eta_1) ) + std::pow( std::max(0.0,q-eta_1)/eta_2 , 1.0/kappa_2 );
  kappa = Dx*q/s;
}


// Muskingham Cunge after Todini
sfc_mct::sfc_mct(std::vector<double> const &param, std::vector<double> const &properties){
  Dx = properties[1];
  grd = properties[2];
  n = param[0];
  ca = param[1];
  sa = std::sin( std::atan(param[1]) );
  B0 = param[2];
}
// internal update
void sfc_mct::update(double const&Q){
  auto Ay = [&](double y){ return( (B0 + y*ca)*y ); };
  auto By = [&](double y){ return( B0 + 2*y*ca ); };
  auto Py = [&](double y){ return( B0 + 2*(y/sa) ); };
  auto Qy = [&](double y){ return( (std::sqrt(grd)/n) * std::pow(Ay(y),(5/3)) / std::pow(Py(y),(2/3)) ); };
  auto cy = [&](double y){ return( 
				   (5/3)* (std::sqrt(grd)/n) * std::pow(Ay(y),(2/3)) / std::pow(Py(y),(2/3)) *
				   ( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*sa) ) ) );
  };
  auto vy = [&](double y){ return( (std::sqrt(grd)/n) * std::pow(Ay(y)/Py(y), 2/3) ); };
  //auto betay = [&](double y){ return( (5/3)*( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*sa) ) ) ); };
  
  // solve for y
  // lower bound of search - start at 0.0 - this flow should be less then Q
  std::pair<double,double> lbnd(0.0, 999.9);
  lbnd.second = Qy(lbnd.first);
  // pick a high value - this flow should be greater then Q
  std::pair<double,double> ubnd(1000.0, 999.9);
  ubnd.second = Qy(ubnd.first);
  if( ubnd.second < Q ){
    Rcpp::Rcout <<"need adaptive mC range" << std::endl;
    Rcpp::Rcout << grd << " " << B0 << " " << ca << " " << sa << " " << n << " " << Dx << std::endl;
    Rcpp::Rcout << Ay(ubnd.first) << " " << Py(ubnd.first) << std::endl;
    Rcpp::Rcout <<"Lower bound:" << lbnd.first << " " << lbnd.second << std::endl;
    Rcpp::Rcout <<"Upper bound:" << ubnd.first << " " << ubnd.second << std::endl;
  }
  int it = 0;
  double y(0.0);
  while( (it <= 1000) and ( ubnd.second - lbnd.second > 1e-3 ) ){
    //double iW = (Q - lbnd.second) / (ubnd.second-lbnd.second);
    //iW = std::max(0.001,std::min(iW,0.999));
    y = (ubnd.first + lbnd.first)/2.0; //(iW*ubnd.first) + (1.0-iW)*lbnd.first;
    double qq = Qy(y);
    if( qq <= Q ){ //bnd.second= z; } else { bnd.first=z; }
      lbnd.first = y;
      lbnd.second = qq;
    }else{
      ubnd.first = y;
      ubnd.second = qq;
    }
    it += 1;
  }
  if( (lbnd.second > Q ) or (ubnd.second < Q) ){
    Rcpp::Rcout << "error in solving for height" << std::endl;
    Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
    Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
  }

  double vel = vy(y);
  double cel = cy(y);
  double tw = By(y);
  double D = Q / (2*tw*grd);
  kappa = Dx / vy(y);
  eta = 0.5*( 1.0 -  ((2*D*vel)/(Dx*cel*cel)) );
  
  // double beta = betay(y);
  // double cel = cy(y);
  // Cs = cel/(beta*Dx); // removed Dt compared to paper
  // if(Q == 0){
  //   Ds = 0.0;
  // }else{
  //   Ds = Q/(beta*By(y)*grd*cel*Dx);
  // }
  
};
// double sfc_mct::fq(double const &s){ return(-999.9); }
// //   double q = eta_1*std::min(s_1,s);
// //   q += eta_2 * std::pow( std::max(0.0,(s-s_1)), kappa_2 );
// //   return( q );
// // }
// double sfc_mct::fs(double const &q, double const &qin ){
//   if( qin < 0.0 ){ return( std::nan("") ); } // flag negative inflows
//   if( q < 0.0 ){ return( std::nan("") ); } // handl case of no outflow
//   double Qref = (q+qin)/2;
//   if( Qref<= 0.0 ){ return(0.0); }
//   internal_update(Qref);
//   return( (1/(2*Cs))*( (1-Ds)*qin + (1+Ds)*q ) ); // can be simplified
// }
// void sfc_mct::update(double &s, double &q, double const &qin, double const &vout,
//   double const &Dt, double const &vtol, int const &max_it){
    
//   double s0 = s + Dt*qin - vout;
//   if( s0 == 0.0 ){ // no stroage so no outflow
//     s = s0;
//     q = 0.0;
//     return;
//   }

//   double q_hat = qin;
//   for(int it=0; it<3; ++it){
//     double Qref = (q_hat + qin)/2.0;
//     internal_update(Qref);
//     q_hat = std::max(0.0, (2*Cs*s0 - (1-Ds)*qin) / (1+Ds+2*Cs*Dt) );
//   }
//   q = q_hat;
//   s = s0 - Dt*q;

// };

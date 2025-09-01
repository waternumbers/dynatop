// each function takes a storage and 
#include "sf.h"

// solve 
sfc::sfc(){ }

std::pair<double,double> sfc::fq( double const &s ){ // set for two partions with different constants
  std::pair<double,double> out(0.0,0.0);
  if (s>0) {
    out.first = k_1*std::min(s,S_1) + k_2*std::max(0.0,s-S_1) ;
    out.second = k_1;
    if( s > S_1 ){ out.second = k_2; }
  }
  return( out );
}
double sfc::fs( double const &q ){ // set for two partions with different time constants
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  return( std::min( q/k_1,S_1 ) + std::max(q - k_1*S_1, 0.0)/k_2 );
}
    
// constant celerity, diffusivity with raf
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  // lower path is linear tank
  k_1 = 1.0/param[2]; // param[2] is t_raf
  S_1 = param[1]; // param[1] is raf storage
  k_2 = param[0] / properties[1]; // velocity divided by Dx
}

// Mannings with raf - assume shallow water for hydraulic radius ~ S/area
// so outside of raf
// v = (S^(2/3) * sqrt(gradient) ) / (n * area^(2/3))
// q = (S^(5/3) * sqrt(gradient) ) / (Dx * n * area^(2/3))
// q = eta * s^(5/3)
// where
// eta = sqrt(gradient) / (Dx * n * area^(2/3))
// Assumes shallow water so wetted perimeter ~ width
sfc_kin::sfc_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const &Dx(properties[1]), &area(properties[0]), &grd(properties[2]);
  double const &n(param[0]);
  k_1 = 1.0/param[2]; // param[2] is raf time constant
  S_1 = param[1]; // param[1] is raf storage
  eta = std::sqrt(grd) / (Dx * n * std::pow(area, 2.0/3.0));
}
std::pair<double,double> sfc_kin::fq(double const &s ){
  std::pair<double,double> out(0.0,0.0);
  if( s>0 ){
    out.first = (k_1*std::min(s,S_1)) + (eta*std::pow(std::max(s-S_1,0.0), 5.0/3.0));
    out.second = k_1;
    if( s > S_1 ){
      out.second = (5.0/3.0)*eta*std::pow( s-S_1, 2.0/3.0 );
    }
  }
  return( out );
}
double sfc_kin::fs(double const &q ){
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  return( std::min( q/k_1,S_1 ) + std::pow( std::max(q - k_1*S_1, 0.0) / eta, 3.0/5.0 ) );
}


// compound channel
sfc_comp::sfc_comp(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[1]);
  k_1 = param[0] / Dx; // velocity divided by length to get q from storage for lower part of channel
  S_1 = param[1] ; // max value of lower store
  k_2 = param[2] / Dx; // velocity divided by length to get q from storage for upper part of channel
}

// arbitary area, discharge relationship
sfc_arb_kin::sfc_arb_kin(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[1]);
  unsigned int n = param.size()/2;
  for(unsigned int ii = 0; ii<n; ++ii){
    s_val.push_back( param[ii]*Dx );
    q_val.push_back( param[ii+n] );
  }
}
std::pair<double,double> sfc_arb_kin::fq(double const &s ){
  std::pair<double,double> out(0.0,0.0);
  if( s > 0 ){
    unsigned int n = s_val.size();
    unsigned int ii = 1;
    while( (s_val[ii] < s) & (ii < (n-1)) ){
      ii += 1;
    };
    out.first = q_val[ii-1] + ( (q_val[ii] - q_val[ii-1])/(s_val[ii]-s_val[ii-1]) )* (s - s_val[ii-1]);
    out.second = ( (q_val[ii] - q_val[ii-1])/(s_val[ii]-s_val[ii-1]) );
  }
  return(out);
}
double sfc_arb_kin::fs(double const &q ){
  if( q<= 0.0 ){ return( 0.0 ); } // handle case of no outflow
  unsigned int n = s_val.size();
  unsigned int ii = 1;
  while( (q_val[ii] < q) & (ii < (n-1)) ){
    ii += 1;
  };
  double s = s_val[ii-1] + ( (s_val[ii] - s_val[ii-1])/(q_val[ii]-q_val[ii-1]) )* (q - q_val[ii-1]);
  return(s);
}

// generic power law with constant parameters
sfc_power_law::sfc_power_law(std::vector<double> const &param, std::vector<double> const &properties){
  k_1 = 1.0 / param[3]; // param[3] is raf time constant
  S_1 = param[2]; // param[2] is raf storage
  kappa = param[1]; // power
  eta = param[0]; // scale
}
std::pair<double,double> sfc_power_law::fq(double const &s ){
  std::pair<double,double> out(0.0,0.0);
  if( s > 0 ){
    out.first = (k_1*std::min(s,S_1)) + (eta * std::pow(std::max(s-S_1,0.0), kappa));
    out.second = k_1;
    if( s > S_1 ){
      out.second = eta*kappa*std::pow(s-S_1, kappa-1.0);
    }
  }
  return( out );
}
double sfc_power_law::fs(double const &q ){
  if( q<= 0.0 ){ return(0); } // handle case of no outflow
  return( std::min( q/k_1,S_1 ) + std::pow( std::max(q - k_1*S_1, 0.0) / eta, 1.0/kappa ) );
}


// Trapezoid channel with Mannings after Todini
sfc_mct::sfc_mct(std::vector<double> const &param, std::vector<double> const &properties){
  Dx = properties[1];
  grd = properties[2];
  n = param[0];
  ca = 1.0 / param[1] ; // cotangent (cot) of side slop angle
  sa = std::sin( std::atan(param[1]) ); // sin of side slope angle
  B0 = param[2]; // bed width
  // auto Ay = [&](double y){ return( (B0 + y*ca)*y ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
  // auto Py = [&](double y){ return( B0 + 2.0*(y/sa) ); };
  // auto Qy = [&](double y){ return( (std::sqrt(grd)/n) * std::pow(Ay(y),(5.0/3.0)) / std::pow(Py(y),(2.0/3.0)) ); };
  // auto dQ_dy = [&](double y){ return( Qy(y)*( (5.0/3.0)*(B0+2.0*ca*y)/Ay(y) - (4.0/3.0)/(sa*Py(y)) ) ); };
}
//helper function
double sfc_mct::Ay(double const&y){ return( (B0 + y*ca)*y ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
double sfc_mct::Py(double const&y){ return( B0 + 2.0*(y/sa) ); };
double sfc_mct::Qy(double const&y){ return( (std::sqrt(grd)/n) * std::pow(Ay(y),(5.0/3.0)) / std::pow(Py(y),(2.0/3.0)) ); };
double sfc_mct::dQ_dy(double const&y){
  if( y == 0.0 ){ return(0.0); } // PJJS test
  return( Qy(y)*( (5.0/3.0)*(B0+2.0*ca*y)/Ay(y) - (4.0/3.0)/(sa*Py(y)) ) ); };
// internal update
std::pair<double,double> sfc_mct::fq(double const&s){
  std::pair<double,double> out(0.0,0.0);
  if( s > 0 ){
    double A = s/Dx;
    double h = ( -B0 + std::sqrt( std::pow(B0,2.0) + 4.0*A*ca ) ) / (2.0*ca) ;
    double dh_ds = 1.0 / (Dx * std::sqrt( std::pow(B0,2.0) + 4.0*A*ca ));
    out.first = Qy(h);
    out.second = dQ_dy(h) * dh_ds;
  }
  return(out);
}

double sfc_mct::fs(double const&q){
  if( q<= 0.0 ){ return(0); } // handle case of no outflow    
  double y = 1.0; //initial estimate
  double e = q - Qy(y);
  double y_old = 100; // previous guess
  int it = 0;
  
  while( (it<100) and (std::abs(y_old - y) > 1e-6) and (std::abs(e)>1e-6) ){
    y_old = y;
    y = std::max(y + (e/dQ_dy(y)) , 0.0);
    if(y == 0 ){
      return(0); //break;
    }
    e = q - Qy(y);
    it +=1;
  }
  
  // Rcpp::Rcout << " y " << y << " q " << q << " e "<< e<< std::endl;
  return( Dx*Ay(y) );

}



// //////////////////////////
// Muskingham Cunge after Todini with two level rectangular channel
sfc_mct_rect::sfc_mct_rect(std::vector<double> const &param, std::vector<double> const &properties){
  //Rcpp:Rcout << "in initialisation" << std::endl;
  // store inputs
  Dx = properties[1];
  double const& grd(properties[2]);
  // parameters
  double const& n = param[0];
  B0 = param[1]; // bed width for rectangular segment
  ca = 1.0 / param[2] ; // cotangent (cot) of side slop angle
  sa = std::sin( std::atan(param[2]) ); // sin of side slope angle
  double const& q_crit = param[3]; // threshold flow
  // functions
 
  // computed values
  beta = std::sqrt(grd) / n;
  y_crit = 1e300; // set large so next part stays within the rectangular channel part
  y_crit = solve_depth(q_crit);
  A_crit = B0 * y_crit;
  //Rcpp:Rcout << "A_crit: " << A_crit << "y_crit: " << y_crit << " Dx: " << Dx << std::endl;
  //Rcpp::Rcout << "y_crit: " << y_crit << " A_crit: " << A_crit << " q_crit: " << q_crit << std::endl;
}
// helper function
double sfc_mct_rect::Ay(double const&y){ return( (B0*y) + std::max(0.0,y-y_crit)*ca*std::max(0.0,y-y_crit) ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
double sfc_mct_rect::Py(double const&y){ return( B0 + 2.0*std::min(y,y_crit) + 2.0*(std::max(0.0,y-y_crit)/sa) ); };
double sfc_mct_rect::Qy(double const&y){ return( (beta * std::pow(Ay(y),(5.0/3.0))) / std::pow(Py(y),(2.0/3.0)) ); };
double sfc_mct_rect::dQ_dy(double const&y){
  double dP_dy = 2;
  if( y> y_crit){ dP_dy = 2.0/sa; }
  return( Qy(y)*( (5.0/3.0)*(B0+2.0*ca*std::max(0.0,y-y_crit))/Ay(y) - (2.0/3.0)*dP_dy/Py(y) ) );
};
// q = A*sqrt(s)*(R^2/3)/n;
// y<yc
// A=b_lower*y;
// wp = b_lower+ 2*y;
// y>yc
// alpha = atan(b_upper);
// A = b_lower*y + (y-yc)*(y-yc)/b_upper;
// wp = b_lower + 2*yc + (y-yc)/sin(alpha);
// celerity
// beta = sqrt(s)/n;
// celerity dq/dA = beta*(R^2/3) + A*beta*(2/3)*(R^-1/3)/wp
// dq/dy = dq/dA dA/dy
// dA/dy = b_lower + 2(y-yc)/b_upper
 
// solve depth for a given outflow
double sfc_mct_rect::solve_depth(double const&Q){
  if(Q==0){
    return(0.0);
  }
  double y = 1.0; //initial estimate
  double e = Q - Qy(y);
  double y_old = 100; // previous guess
  int it = 0;
  //Rcpp:Rcout << "it: " << it << " y " << y << " q " << Qy(y) << " e "<< e << std::endl;
  while( (it<100) and (std::abs(y_old - y) > 1e-6) and (std::abs(e)>1e-6) ){
    y_old = y;
    y = std::max(y + (e/dQ_dy(y)) , 0.0);
    if(y == 0 ){
      return(0); //break;
    }
    e = Q - Qy(y);
    it +=1;
    //Rcpp:Rcout << "it: " << it << " y " << y << " q " << Qy(y) << " e "<< e << std::endl;
  }
  
  //Rcpp:Rcout << " y " << y << " q " << Qy(y) << " e "<< e << std::endl;
  return( y );
}

// internal update
std::pair<double,double> sfc_mct_rect::fq(double const&s){
  std::pair<double,double> out(0.0,0.0);
  if( s > 0 ){
    // solve for height given the cross sectional area
    double A = s/Dx;
    double y(A/B0), dh_ds(1.0/B0*Dx);
    if(A > A_crit){ // then some trapezoid part
      y = y_crit + ( (-B0 + std::sqrt( std::pow(B0,2.0) + 4*(A-A_crit)*ca )) / (2.0*ca) );
      if( y < y_crit){
	Rcpp::Rcout <<"Negative h_2 " << y << " "<< y_crit << std::endl;
	y = y_crit;
      }
      dh_ds = 1.0 / (Dx * std::sqrt( std::pow(B0,2.0) + 4.0*(A-A_crit)*ca ));
    }
    
    // compute flow
    out.first = Qy(y);
    out.second = dQ_dy(y)*dh_ds;
  }  
  return( out );
};

// internal update
double sfc_mct_rect::fs(double const&q){
  //Rcpp:Rcout << "in fS" << std::endl;
  if( q<= 0.0 ){ return(0.0); }
  double y = solve_depth(q);
  return( Ay(y)*Dx );
};


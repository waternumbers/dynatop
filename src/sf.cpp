// each function takes a storage and 
#include "sf.h"

// fT returns a such that q=aS => a = q/S
// fS returns S such that q=aS

// solve 
sfc::sfc(){ }
double sfc::fT( double const &s ){ // set for two partions with different time constants
  if( s<= 0.0 ){ return(0.0); } // handle case of no outflow
  return( ( std::min(s,S_1)/T_1 + std::max(0.0,s-S_1)/T_2 ) / s );
}
double sfc::fS( double const &q ){ // set for two partions with different time constants
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  return( std::min( q*T_1,S_1 ) + T_2*std::max(q - S_1/T_1, 0.0) );
}
    

// constant celerity, diffusivity with raf
sfc_cnst::sfc_cnst(std::vector<double> const &param, std::vector<double> const &properties){
  // lower path is linear tank
  T_1 = param[2]; // param[2] is t_raf
  S_1 = param[1]; // param[1] is raf storage
  T_2 = param[0] / properties[1]; // constant velocity divided by Dx
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
  T_1 = param[2]; // param[2] is raf time constant
  S_1 = param[1]; // param[1] is raf storage
  eta = std::sqrt(grd) / (Dx * n * std::pow(area, 2.0/3.0));
}
double sfc_kin::fT(double const &s ){
  if( s<= 0.0 ){ return(0.0); } // handle case of no outflow
  double q = (std::min(s,S_1)/T_1) + (std::pow(std::max(s-S_1,0.0), 5.0/3.0)*eta);
  return( q/s );
}
double sfc_kin::fS(double const &q ){
  if( q<= 0.0 ){ return(0.0); } // handle case of no outflow
  return( std::min( q*T_1,S_1 ) + std::pow( std::max(q - S_1/T_1, 0.0) / eta, 3/5 ) );
}


// compound channel
sfc_comp::sfc_comp(std::vector<double> const &param, std::vector<double> const &properties){
  double const& Dx(properties[1]);
  T_1 = param[0] / Dx; // velocity divided by length to get q from storage for lower part of channel
  S_1 = param[1] ; // max value of lower store
  T_2 = param[2] / Dx; // velocity divided by length to get q from storage for upper part of channel
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
double sfc_arb_kin::fT(double const &s ){
  if( s<= 0.0 ){ return(0.0); } // handle case of no outflow
  unsigned int n = s_val.size();
  unsigned int ii = 1;
  while( (s_val[ii] < s) & (ii < (n-1)) ){
    ii += 1;
  };
  double q = q_val[ii-1] + ( (q_val[ii] - q_val[ii-1])/(s_val[ii]-s_val[ii-1]) )* (s - s_val[ii-1]);
  return( q/s );
}
double sfc_arb_kin::fS(double const &q ){
  if( q<= 0.0 ){ return( 0 ); } // handle case of no outflow
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
  T_1 = param[3]; // param[3] is raf time constant
  S_1 = param[2]; // param[2] is raf storage
  kappa = param[1]; // power
  eta = param[0]; // scale
}
double sfc_power_law::fT(double const &s ){
  if( s<= 0.0 ){ return(0.0); } // handle case of no outflow
  double q = (std::min(s,S_1)/T_1) + (eta * std::pow(std::max(s-S_1,0.0), kappa));
  return( q/s );
}
double sfc_power_law::fS(double const &q ){
  if( q<= 0.0 ){ return(0); } // handle case of no outflow
  return( std::min( q*T_1,S_1 ) + std::pow( std::max(q - S_1/T_1, 0.0) / eta, 1/kappa ) );
}

// Trapezoid channel with Mannings after Todini
sfc_mct::sfc_mct(std::vector<double> const &param, std::vector<double> const &properties){
  Dx = properties[1];
  grd = properties[2];
  n = param[0];
  ca = 1.0 / param[1] ; // cotangent (cot) of side slop angle
  sa = std::sin( std::atan(param[1]) ); // sin of side slope angle
  B0 = param[2]; // bed width
}
// internal update
double sfc_mct::fT(double const&s){
  if( s<= 0.0 ){ return(0.0); } // handle case of no outflow
  
  auto Ay = [&](double y){ return( (B0 + y*ca)*y ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
  auto Py = [&](double y){ return( B0 + 2*(y/sa) ); };
  auto Qy = [&](double y){ return( (std::sqrt(grd)/n) * std::pow(Ay(y),(5/3)) / std::pow(Py(y),(2/3)) ); };
  
  // solve for height given the cross sectional area
  double A = s/Dx;
  double h = ( -B0 + std::sqrt( std::pow(B0,2.0) + 4*A*ca ) ) / (2*ca) ;
  if( h<0.0){
    Rcpp::Rcout <<"Negative h " << h << std::endl;
    h = 0.0;
  }
  // compute flow
  double q = Qy(h);
  return( q/s );
};
double sfc_mct::fS(double const&q){
  if( q<= 0.0 ){ return(0); } // handle case of no outflow
  // find area giving outflow
  auto Ay = [&](double y){ return( (B0 + y*ca)*y ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
  auto Py = [&](double y){ return( B0 + 2*(y/sa) ); };
  auto Qy = [&](double y){ return( (std::sqrt(grd)/n) * std::pow(Ay(y),(5/3)) / std::pow(Py(y),(2/3)) ); };
  auto dQ_dy = [&](double y){ return( Qy(y)*( (5/3)*(B0+2*ca*y)/Ay(y) - (4/3)/(sa*Py(y)) ) ); };
    
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
  // store inputs
  Dx = properties[1];
  double const& grd(properties[2]);
  // parameters
  double const& n = param[0];
  B0 = param[1]; // bed width for rectangular segment
  ca = 1.0 / param[2] ; // cotangent (cot) of side slop angle
  sa = std::sin( std::atan(param[2]) ); // sin of side slope angle
  double const& q_crit = param[3]; // threshold flow
  // computed values
  beta = std::sqrt(grd) / n;
  y_crit = 1e300; // set large so next part stays within the rectangular channel part
  y_crit = solve_depth(q_crit);
  A_crit = B0 * y_crit;
}

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
  auto Ay = [&](double y){ return( (B0*y) + std::max(0.0,y-y_crit)*ca*std::max(0.0,y-y_crit) ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
  auto Py = [&](double y){ return( B0 + 2*std::min(y,y_crit) + 2*(std::max(0.0,y-y_crit)/sa) ); };
  auto Qy = [&](double y){ return( (beta * std::pow(Ay(y),(5/3))) / std::pow(Py(y),(2/3)) ); };
  auto dQ_dy = [&](double y){
    double dP_dy = 2;
    if( y> y_crit){ dP_dy = 2/sa; }
    return( Qy(y)*( (5/3)*(B0+2*ca*std::max(0.0,y-y_crit))/Ay(y) - (2/3)*dP_dy/Py(y) ) );
  };

  double y = 1.0; //initial estimate
  double e = Q - Qy(y);
  double y_old = 100; // previous guess
  int it = 0;
  
  while( (it<100) and (std::abs(y_old - y) > 1e-6) and (std::abs(e)>1e-6) ){
    y_old = y;
    y = std::max(y + (e/dQ_dy(y)) , 0.0);
    if(y == 0 ){
      return(0); //break;
    }
    e = Q - Qy(y);

    it +=1;
    
  }
  
  //Rcpp::Rcout << " y " << y << " q " << Qy(y) << " e "<< e<< std::endl;
  return( y );
}

// internal update
double sfc_mct_rect::fT(double const&s){
  if( s<= 0.0 ){ return(0.0); } // handle case of no outflow
  
  // solve for height given the cross sectional area
  double A = s/Dx;
  double h_1(0), h_2(0);
  if(A > A_crit){ // then some trapezoid part
    h_1 = A_crit / B0;
    h_2 = -B0 + std::sqrt( std::pow(B0,2.0) + 4*(A-A_crit)*ca ) / (2*ca) ;
  }else{
    h_1 = A/B0;
  }

  // compute flow
  double Wp = B0 + (2*h_1) + 2*h_2/sa;
  double q = beta * std::pow(A, 5/3) / std::pow(Wp, 2/3);
  return( q/s );
};
// internal update
double sfc_mct_rect::fS(double const&q){
  if( q<= 0.0 ){ return(0.0); }
  auto Ay = [&](double y){ return( (B0*y) + std::max(0.0,y-y_crit)*ca*std::max(0.0,y-y_crit) ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
  double y = solve_depth(q);

  return( Ay(y)*Dx );
};


// // //////////////////////////
// // Muskingham Cunge after Todini with two level rectangular channel
// sfc_mct_rect::sfc_mct_rect(std::vector<double> const &param, std::vector<double> const &properties){
//   // store inputs
//   Dx = properties[1];
//   grd = properties[2];
//   // parameters
//   double const& n = param[0];
//   b_lower = param[1];
//   tan_alpha = param[2];
//   q_crit = param[3]; // threshold flow
//   // computed values
//   beta = std::sqrt(grd) / n;
//   y_crit = 1e300; // set very large to start with
//   y_crit = solve_depth(q_crit);
//   sin_alpha = std::sin( std::atan(tan_alpha) );
//   // with R = y approximation
//   //beta_lower = std::sqrt(grd) * std::pow(b_lower, 2.0/3.0) / n_lower; // Q =beta * y^{5/3)
//   //beta_upper = std::sqrt(grd) * std::pow(b_upper-b_lower, 2.0/3.0) / n_upper;
//   //y_crit = std::pow( q_crit / beta_lower , 3.0/5.0 ); // level eqivilent to q_crit
// }

// // q = A*sqrt(s)*(R^2/3)/n;
// // y<yc
// // A=b_lower*y;
// // wp = b_lower+ 2*y;
// // y>yc
// // alpha = atan(b_upper);
// // A = b_lower*y + (y-yc)*(y-c)/b_upper;
// // wp = b_lower + 2*yc + (y-yc)/sin(alpha);
// // celerity
// // beta = sqrt(s)/n;
// // dq/dA = beta*(R^2/3) + A*beta*(2/3)*(R^-1/3)/wp
  
 
// // solve depth for the flow
// double sfc_mct_rect::solve_depth(double const&Q){
//   std::pair<double,double> lbnd, ubnd;
//   double y, qq;
//   if( Q > q_crit ){ // in trapezoid part
//     lbnd.first = y_crit;
//     lbnd.second = q_crit;
//   }else{
//     lbnd.first = 0.0;
//     lbnd.second = 0.0;
//   }
//   int it(0);
//   y = lbnd.first;
//   qq = lbnd.second;
//   while( (it <=100) and qq < Q ){
//     y = 2*y + 0.1;
//     double ytilde = std::max(y-y_crit,0.0);
//     double A = b_lower*y + ytilde*ytilde/tan_alpha;
//     double Wp =  b_lower + 2*std::min(y_crit,y) + 2*( ytilde/sin_alpha );
//     qq = beta * A * std::pow((A/Wp),(2.0/3.0));
//   }
//   ubnd.first = y;
//   ubnd.second = qq;
  
  
//   if( (lbnd.second > Q ) or (ubnd.second < Q) ){
//     Rcpp::Rcout << "error in solving for height at start" << std::endl;
//     Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
//     Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
//     Rcpp::Rcout << y_crit << " " << beta << " " << b_lower << std::endl;
//   }

//   while( (it <= 100) and ( ubnd.second - lbnd.second > 1e-6 ) ){
//     y = (ubnd.first + lbnd.first)/2.0;
//     double ytilde = std::max(y-y_crit,0.0);
//     double A = b_lower*y + ytilde*ytilde/tan_alpha;
//     double Wp =  b_lower + 2*std::min(y_crit,y) + 2*( ytilde/sin_alpha );
//     qq = beta * A * std::pow(A/Wp, 2.0/3.0);
//     if( qq <= Q ){ 
//       lbnd.first = y;
//       lbnd.second = qq;
//     }else{
//       ubnd.first = y;
//       ubnd.second = qq;
//     }
//     it += 1;
//   }
//   y = (ubnd.first + lbnd.first)/2.0;
//   if( (lbnd.second > Q ) or (ubnd.second < Q) ){
//     Rcpp::Rcout << "error in solving for height" << std::endl;
//     Rcpp::Rcout << it << std::endl;
//     Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
//     Rcpp::Rcout << lbnd.second-Q  << " " << Q-ubnd.second << std::endl;
//     Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
//   }
//   if( it > 100 ){
//     Rcpp::Rcout << "max iterations" << std::endl;
//     Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
//     Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
//   }
//   return(y);
// }
// // internal update
// void sfc_mct_rect::update(double const&Q){
  
//   double y = solve_depth(Q);
//   ##
//   double ytilde = std::max(y-y_crit,0.0);
//   double tw = b_lower + ytilde*tan_alpha/sin_alpha;
//   double area =  b_lower*y + ytilde*ytilde/tan_alpha;
//   double Wp =  b_lower + 2*std::max(y-y_crit,0.0) + 2*( ytilde/sin_alpha );
//   if( area <= 1e-6 ){
//     kappa = 0.0; //-999.0;
//     eta = 0.5;
//   }else{
//     double vel = Q/area;
//     double R = area/Wp;
//     // compute celerity...
//     double cel = (beta*std::pow(R,2.0/3.0)) + ( area*(beta*(2/3))*(std::pow(R,-1.0/3.0))/Wp );
//     double D = Q / (2*tw*grd);
//     kappa = Dx / vel;
//     eta = 0.5*( 1.0 -  ((2*D*vel)/(Dx*cel*cel)) );
//     eta = std::max(eta,0.0);
//     //eta = 0; //.5; //PJS test
//     // if( kappa < 0 | eta < 0 | vel > 10 ){ Rcpp::Rcout << "kappa " << kappa << " eta " << eta << "vel " << vel << std::endl; }
//   }
// };


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
  ca = 1.0 / param[1] ; // grad = tan(theta) = sin(theta)/cos(theta) was originally param[1]
  sa = std::sin( std::atan(param[1]) );
  B0 = param[2];
}
// internal update
void sfc_mct::update(double const&Q){
  auto Ay = [&](double y){ return( (B0 + y*ca)*y ); }; // y = x*sin(theta) => x*cos(theta) = y *cos(theta)/sin(theta) = y/grad = y * cot(theta)
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


// //////////////////////////
// Muskingham Cunge after Todini with two level rectangular channel
sfc_mct_rect::sfc_mct_rect(std::vector<double> const &param, std::vector<double> const &properties){
  // store inputs
  Dx = properties[1];
  grd = properties[2];
  // parameters
  double const& n = param[0];
  b_lower = param[1];
  tan_alpha = param[2];
  q_crit = param[3]; // threshold flow
  // computed values
  beta = std::sqrt(grd) / n;
  y_crit = 1e300; // set very large to start with
  y_crit = solve_depth(q_crit);
  sin_alpha = std::sin( std::atan(tan_alpha) );
  // with R = y approximation
  //beta_lower = std::sqrt(grd) * std::pow(b_lower, 2.0/3.0) / n_lower; // Q =beta * y^{5/3)
  //beta_upper = std::sqrt(grd) * std::pow(b_upper-b_lower, 2.0/3.0) / n_upper;
  //y_crit = std::pow( q_crit / beta_lower , 3.0/5.0 ); // level eqivilent to q_crit
}

// q = A*sqrt(s)*(R^2/3)/n;
// y<yc
// A=b_lower*y;
// wp = b_lower+ 2*y;
// y>yc
// alpha = atan(b_upper);
// A = b_lower*y + (y-yc)*(y-c)/b_upper;
// wp = b_lower + 2*yc + (y-yc)/sin(alpha);
// celerity
// beta = sqrt(s)/n;
// dq/dA = beta*(R^2/3) + A*beta*(2/3)*(R^-1/3)/wp
  
 
// solve depth for the flow
double sfc_mct_rect::solve_depth(double const&Q){
  std::pair<double,double> lbnd, ubnd;
  double y, qq;
  if( Q > q_crit ){ // in trapezoid part
    lbnd.first = y_crit;
    lbnd.second = q_crit;
    double step = y_crit + 0.1;
    ubnd = lbnd;
    int it(0);
    y = ubnd.first;
    qq = ubnd.second;
    while( (it <=100) and qq < Q ){
      y =+ step;
      double A = b_lower*y + (y-y_crit)*(y-y_crit)/tan_alpha;
      doubel Wp =  b_lower + 2*y_crit + ( (y-y_crit)/sin_alpha );
      qq = beta * A * ((A/Wp)^(2.0/3.0));
    }
    ubnd.first = y;
    ubnd.second = qq;
  }else{
    lbnd.first = 0.0;
    lbnd.second = 0.0;
    ubnd.first = y_crit;
    ubnd.second = q_crit;
  }
  
  if( (lbnd.second > Q ) or (ubnd.second < Q) ){
    Rcpp::Rcout << "number of iterations is " << it << std::endl;
    Rcpp::Rcout << "error in solving for height at start" << std::endl;
    Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
    Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
    Rcpp::Rcout << y_crit << " " << beta_lower << " " << b_lower << std::endl;
  }
  it = 0;
  while( (it <= 100) and ( ubnd.second - lbnd.second > 1e-6 ) ){
    y = (ubnd.first + lbnd.first)/2.0;
    double ytilde = std::max(y-y_crit,0.0);
    double A = b_lower*y + ytilde*ytilde/tan_alpha;
    doubel Wp =  b_lower + 2*std::min(y_crit,y) + 2*( ytilde/sin_alpha );
    qq = beta * A * ((A/Wp)^(2.0/3.0));
    if( qq <= Q ){ 
      lbnd.first = y;
      lbnd.second = qq;
    }else{
      ubnd.first = y;
      ubnd.second = qq;
    }
    it += 1;
  }
  y = (ubnd.first + lbnd.first)/2.0;
  if( (lbnd.second > Q ) or (ubnd.second < Q) ){
    Rcpp::Rcout << "error in solving for height" << std::endl;
    Rcpp::Rcout << it << std::endl;
    Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
    Rcpp::Rcout << lbnd.second-Q  << " " << Q-ubnd.second << std::endl;
    Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
  }
  if( it > 100 ){
    Rcpp::Rcout << "max iterations" << std::endl;
    Rcpp::Rcout << lbnd.second << " " << Q << " " << ubnd.second << std::endl;
    Rcpp::Rcout << lbnd.first << " " << ubnd.first << std::endl;
  }
  return(y);
}
// internal update
void sfc_mct_rect::update(double const&Q){
  double y = solve_depth(Q);
  double ytilde = std::max(y-y_crit,0.0);
  double tw = b_lower + ytilde*tan_alpha/sin_alpha;
  double area =  b_lower*y + ytilde*ytilde/tan_alpha;
  doubel Wp =  b_lower + 2*std::max(y-y_crit,0.0) + 2*( ytilde/sin_alpha );
  if( area == 0.0 ){
    kappa = -999.0;
    eta = 0.5;
  }else{
    double vel = Q/area;
    double R = area/Wp;
    // compute celerity...
    double cel = (beta*(R^2/3)) + ( A*(beta*(2/3))*(R^(-1/3))/Wp );
    double D = Q / (2*tw*grd);
    kappa = Dx / vel;
    eta = 0.5*( 1.0 -  ((2*D*vel)/(Dx*cel*cel)) );
    eta = std::max(eta,0.0);
    //eta = 0; //.5; //PJS test
    // if( kappa < 0 | eta < 0 | vel > 10 ){ Rcpp::Rcout << "kappa " << kappa << " eta " << eta << "vel " << vel << std::endl; }
  }
};


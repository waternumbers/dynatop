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
  sin_alpha = std::sin( std::atan(tan_alpha) );
  beta_rect = std::sqrt(grd) / n;
  beta_tri = std::sqrt(grd) / n;
  A_crit = 1e300; // set very large to start with
  A_crit = solve_area(q_crit);
  ca = 1/tan_alpha;
  sa = std::sin( std::atan(tan_alpha));
}

// solve for flow and celerity for a given area
// internal update
std::vector<double> sfc_mct_rect::fq(double const&A){
  std::vector<double> out = {0.0,0.0,0.0};
  // comput depth assuming in rectangular part
  double y(A/b_lower), Dy(0.0); // initial assuming in rectangular part
  if( A > A_crit ){
    double y_crit = A_crit/b_lower;
    // include the trapezoid part
    y = y_crit + ( (-b_lower + std::sqrt( std::pow(b_lower,2.0) + 4*(A-A_crit)*ca )) / (2.0*ca) );
    Dy = y-y_crit;
    if( y < y_crit){
      Rcpp::Rcout <<"Negative h_2 " << y << " "<< y_crit << std::endl;
      y = y_crit;
    }
  }
  // compute flux
  double qrect = beta_rect*std::pow(b_lower*y,(5.0/3.0)) /
    std::pow(b_lower+2*(y-Dy),(2.0/3.0));
  double qtri = beta_tri * std::pow(ca/2,5.0/3.0) * std::pow(sa,2.0/3.0) * std::pow(Dy, 8.0/3.0);
  out[0] = qrect + 2*qtri;
  // compute celerity
  out[1] = (5.0/3.)*(qrect/y) - (4.0/3.0)*(qrect/(b_lower+2*y));
  if(Dy > 0.0 ){ out[1] += (16.0/3.0)*(qtri/Dy); }
  out[1] = out[1] / (b_lower + 2*ca*Dy);
  // compute top width
  out[2] = b_lower + 2*Dy*ca;
  
  return(out);
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
 
// solve depth for the flow
double sfc_mct_rect::solve_area(double const&Q){
  if(Q==0){
    return(0.0);
  }
  // this use the ridder algorithm which should converge quickly while being robust
  // search for s_sz
  double lb(0.0), ub(100);
  double qlb = fq(lb)[0];
  double qub = fq(ub)[0];
  
  if( qub > 0.0 ){
    int it(0.0);
    while( (Q > qub) and (it < 1000) ){
      lb = ub;
      qlb = qub;
      ub += ub + 3;
      qub = fq(ub)[0];
      it +=1;
    }
  }
  
  // shrink back to find solution
  int it(0);
  double z((lb+ub)/2.0);
  double qz = fq(z)[0];
  while( (std::abs(Q-qz) > 0.1) and (it < 1000) ){
    // bisection 
    z = (ub+lb)/2.0;
    qz = fq(z)[0];
    if( Q > qz ){
      lb = z;
      qlb = qz;
    }else{
      ub = z;
      qub = qz;
    }
    it += 1; 
  }
  if((it == 1000) and (std::abs(Q-qz) > 0.1)){
    Rcpp::warning("HRU %i SZ: No solution found within %i iterations. Difference between bounds is %d",
		  it, it, ub - lb); //bnd.second - bnd.first);
    //Rcpp::Rcout << "id: " << id << " iter: " << it << " diff: " << ub - lb << " Hzu: " << Hzu << std::endl;
  }
  //Rcpp::Rcout << "id: " << id << " iter: " << it <<std::endl;
  
  return( z );
}

// internal update
void sfc_mct_rect::update(double const&Q){
  
  double A = solve_area(Q);
  std::vector<double> qSum = fq(A);
  
  if( A <= 1e-6 ){
    kappa = 0.0; //-999.0;
    eta = 0.5;
  }else{
    double vel = Q/A;
    double const &cel = qSum[1];
    double const &tw = qSum[2];
    double D = Q / (2*tw*grd);
    //Rcpp::Rcout << Q << " " << A << " " << cel << " " << tw << std::endl;
    kappa = Dx / vel;
    eta = 0.5*( 1.0 -  ((2*D*vel)/(Dx*cel*cel)) );
    eta = std::max(eta,0.0);
  }
};

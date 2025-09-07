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
  double const& n_out_of_bank = param[4];
  // functions
 
  // computed values
  beta_rect = std::sqrt(grd) / n;
  beta_tri = std::sqrt(grd) / n_out_of_bank;
  s_crit = 1e300; // set large so next part stays within the rectangular channel part
  s_crit = solve_storage(q_crit);
  y_crit = s_crit / (Dx*B0);
}
// solve depth for a given outflow
double sfc_mct_rect::solve_storage(double const&Q){
  if(Q==0){
    return(0.0);
  }
  
  
  // this use the ridder algorithm which should converge quickly while being robust
  // search for s_sz
  double lb(0.0), ub(100*Dx);
  double qlb = fq(lb).first;
  double qub = fq(ub).first;
  
  if( qub > 0.0 ){
    int it(0.0);
    while( (Q > qub) and (it < 1000) ){
      lb = ub;
      qlb = qub;
      ub += ub + 3;
      qub = fq(ub).first;
      it +=1;
    }
  }
      
  // shrink back to find solution
  int it(0);
  double z((lb+ub)/2.0);
  double qz = fq(z).first;
  while( (std::abs(Q-qz) > 0.1) and (it < 1000) ){
    // bisection 
    z = (ub+lb)/2.0;
    qz = fq(z).first;
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
  
  
  // newton approach - a bit unstable - maybe gradient is wrong??
  // double s = 1000.0; //initial estimate of storage
  // std::pair<double,double> q = fq(s);
  // double e = Q - q.first;
  // double s_old = 100.2342; // previous guess
  // int it = 0;
  // bool show(false);
  // if( std::abs(Q-45.3) <1e-6 ){
  //   show = true;
  // }
  // if(show){
  //   Rcpp::Rcout << "Q: " << Q << " it: " << it << " s: " << s << " q: " << q.first << " dqds: " << q.second << " e: "<< e << std::endl;
  // }
  // while( (it<100) and (std::abs(s_old - s) > 1e-6) and (std::abs(e)>1e-6) ){
  //   s_old = s;
  //   s = std::max(s + (e/q.second) , 1e-6);
  //   // if(s == 0 ){
  //   //   return(0); //break;
  //   // }
  //   q = fq(s);
  //   e = Q - q.first;
  //   it +=1;
  //   if(show){
  //     Rcpp::Rcout << "it: " << it << " sold: " << s_old <<" s: " << s << " q: " << q.first << " dqds: " << q.second << " e: "<< e << std::endl;
  //   }
  // }
  
  //Rcpp:Rcout << " s " << s << " q " << Qs(s) << " e "<< e << std::endl;
  return( z );
}

// internal update
std::pair<double,double> sfc_mct_rect::fq(double const&s){
  std::pair<double,double> out(0.0,0.0);
  if( s > 0 ){
    // solve for height given the cross sectional area
    double A = s/Dx;
    // comput depth assuming in rectangular part
    double y(A/B0), Dy(0.0), dy_ds(1.0/(B0*Dx)); // initial assuming in rectangular part
    if( s > s_crit ){
      // include the trapezoid part
      double A_crit = s_crit / Dx;  
      y = y_crit + ( (-B0 + std::sqrt( std::pow(B0,2.0) + 4*(A-A_crit)*ca )) / (2.0*ca) );
      Dy = y-y_crit;
      if( y < y_crit){
	Rcpp::Rcout <<"Negative h_2 " << y << " "<< y_crit << std::endl;
	y = y_crit;
      }
      dy_ds = 1.0 / (Dx * std::sqrt( std::pow(B0,2.0) + 4.0*(A-A_crit)*ca )); // TODO check
    }
    // compute flux
    double qrect = beta_rect*std::pow(B0*y,(5.0/3.0)) /
      std::pow(B0+2*(y-Dy),(2.0/3.0));
    double qtri = beta_tri * std::pow(ca/2,5.0/3.0) * std::pow(sa,2.0/3.0) * std::pow(Dy, 8.0/3.0);
    out.first = qrect + 2*qtri;
    // compute gradient
    out.second = (5.0/3.0)*(qrect/y);
    if( s > s_crit ){
      out.second += (16.0/3.0)*qtri/Dy;
    }else{
      out.second -= (4.0/3.0)*(qrect/(B0+2*y));
    }
    
    out.second *= dy_ds;
  }
  return(out);
};
  
    
//     // solve flow in rectangular section
//     double A = B0*y;
//     double P = B0 + 2*std::min(y,y_crit);
//     out.first = beta * std::pow(A,(5.0/3.0)) / std::pow(P,(2.0/3.0));
//     out.second = -999.9;
//     if(y > y_crit){
//       A = 0.5*std::max(0.0,y-y_crit)*ca*std::max(0.0,y-y_crit);
//       P = std::max(0,y-y_crit)/sa;
//       out.first += 2 * beta * std::pow(A,(5.0/3.0)) / std::pow(P,(2.0/3.0));
//       out.second = -9999.9;
//     }

//   }  
//   return( out );
// };

// internal update
double sfc_mct_rect::fs(double const&q){
  //Rcpp:Rcout << "in fS" << std::endl;
  if( q<= 0.0 ){ return(0.0); }
  double y = solve_storage(q);
  return( y ); //Ay(y)*Dx );
};


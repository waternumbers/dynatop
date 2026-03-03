#include "Rcpp.h"
#include <vector>

// [[Rcpp::export]]
std::vector<double> newton(std::vector<double>& states, const double& p, const double& e_p,
			   const double& q_sf_in, const double& q_sz_in,
			   const double& Dt){
  double &s_sf(states[0]), &s_rz(states[1]), &s_uz(states[2]), &s_sz(states[3]);

  // hard wire parameters for testing
  // sf param
  std::vector<double> a = {0.0,0.0,1.0}, b = {0.0,0.0,5.0/3.0}, s = {0.0,0.0};
  // rz param
  double s_rz_max(0.1);
  //uz, sz param
  double K_0(2e-5), m(0.2), D(1000);

  // hard wire properties for testing
  double W_beta(1), A(1);
  
  // vertical fluxes
  double v_sf_rz, v_rz_uz, v_uz_sz;

  // lateral fluxes
  double q_sz, q_sf;

  // hard wire tolerence for testing
  double tol(1e-8);
  
  // ///////////////////////
  // Downward Pass

  v_sf_rz = s_sf + Dt*q_sf_in;
  v_rz_uz = std::min(A*K_0*Dt, std::max(0.0,s_rz + Dt*(p-e_p) + v_sf_rz - s_rz_max));

  // solve saturated zone
  double phi = A*K_0*std::exp(-s_sz/m);
  
  auto Hz = [&](double &z){
    double omega = std::min(1.0, (s_uz + v_rz_uz)/(z + Dt*phi));
    double q = W_beta*m*K_0*std::exp(-z/m); //TODO THID IS WRONG - add D
    // Rcpp::Rcout << "q z " << q << " " << z  << std::endl;
    std::pair<double,double>h = h = {s_sz + Dt*(q - q_sz_in - phi*omega) - z,
				     -(q/m) - 1};
    if(omega < 1.0){
      h.second += omega*((Dt*phi)/(z+Dt*phi));
    }
    return h;
  };
  
  // test 0
  double z = 0;
  std::pair<double,double> h = Hz(z);
  if( h.first <= 0.0 ){
    q_sz = W_beta*m*K_0;
    v_uz_sz = s_sz + Dt*(q_sz - q_sz_in);
    s_sz = z;
  }else{
    // not saturated
    std::pair<double,double> rng = {0,D};
    //Rcpp::Rcout << "range " << rng.first << " " << rng.second << std::endl;
    z = 0.0; //s_sz;
    h = Hz(z);
    if(h.first >= 0.0){rng.first = z;}
    if(h.first <= 0.0){rng.second = z;}
    int it(0);
    //Rcpp::Rcout << "h " << h.first << " " << h.second << std::endl;
    //Rcpp::Rcout << "range " << rng.first << " " << rng.second << std::endl;
    while( ((h.first>0.0) or (-h.first>tol)) and (it<100) ){
      double zz = z;
      z = z - (h.first/h.second);
      if( z < rng.first ){ z = (rng.first + zz)/2; }
      if( z > rng.second ){ z = (rng.second + zz)/2; }
      h = Hz(z);
      if(h.first>=0.0){rng.first = z;}
      if(h.first<=0.0){rng.second = z;}
      it += 1;
      //Rcpp::Rcout << "it " <<it  << std::endl;
      //Rcpp::Rcout << "z " << z  << std::endl;
      //Rcpp::Rcout << "h " << h.first << " " << h.second << std::endl;
      //Rcpp::Rcout << "range " << rng.first << " " << rng.second << std::endl;
    }
    //Rcpp::Rcout << "z " << z  << std::endl;
    v_uz_sz = Dt*phi*std::min(1.0, (s_uz + v_rz_uz)/(z + Dt*phi));
    q_sz = (z - s_sz + v_uz_sz + q_sz_in) / Dt;
    s_sz = z;

    Rcpp::Rcout << "sz iteration " << it << std::endl;
  }

  //Rcpp::Rcout << s_sz << std::endl;
  
  
  
  //upward pass
  z = std::min(s_sz,s_uz + v_rz_uz - v_uz_sz);
  v_rz_uz = z - s_uz + v_uz_sz;
  s_uz = z;

  v_sf_rz = std::min(v_sf_rz,s_rz_max - s_rz - Dt*(p-e_p) + v_rz_uz);
  z = ( s_rz_max / (s_rz_max + e_p*Dt) ) * (s_rz + Dt*p + v_sf_rz - v_rz_uz);
  double v_ep = s_rz + v_sf_rz - v_rz_uz + Dt*p - z;
  s_rz = z;

  Rcpp::Rcout << "v_sf_rz " << v_sf_rz << std::endl;
  
  auto Sw = [&](double w){
    double q(0.0), dq(0.0);
    if(w <= s[0]){
      q = a[0]*std::pow(w,b[0]);
      dq = a[0]*b[0]*std::pow(w,b[0]-1.0);
      //dq = b[0]*q/w;
    }else{
      if(w <= s[1]){
	q = a[0]*std::pow(s[0],b[0]) + a[1]*std::pow((w-s[0]),b[1]);
	dq = a[1]*b[1]*std::pow(w-s[0],b[1]-1.0);
	//dq = b[1]*q/(w-s[0]);
      }else{
	q = a[0]*std::pow(s[0],b[0]) + a[1]*std::pow((s[1]-s[0]),b[1]) + a[2]*std::pow(w-s[1],b[2]);
	dq = a[2]*b[2]*std::pow(w-s[1],b[2]-1.0);
        //dq = b[2]*q /(w-s[1]);
      }
    }
    
    // Rcpp::Rcout << "w is " << w << std::endl;
    // Rcpp::Rcout << "q is " << q << std::endl;
    // Rcpp::Rcout << "s_sf is " << s_sf << std::endl;
    // Rcpp::Rcout << "q_sf_in is " << q_sf_in << std::endl;
    // Rcpp::Rcout << "Dt is " << Dt << std::endl;
    // Rcpp::Rcout << std::endl;
    
    std::pair<double,double> s = {
      s_sf + Dt*(q_sf_in - q) - v_sf_rz - w,
      - Dt*dq - 1
    };
    return s;
  };
    
  z = s_sf;
  h = Sw(z);
  Rcpp::Rcout << "z " << z  << std::endl;
  Rcpp::Rcout << "h " << h.first << " " << h.second << std::endl;
  std::pair<double,double> rng = {0.0,1e30};
  if(h.first>=0.0){rng.first = z;}
  if(h.second<=0.0){rng.second = z;}
  int it(0);
  while( ((h.first>0) and (-h.first>tol)) or (it<10000) ){
    double zz = z;
    z = z - (h.first/h.second);
    if( z < rng.first ){ z = (rng.first + zz)/2; }
    if( z > rng.second ){ z = (rng.second + zz)/2; }
    h = Sw(z);
    if(h.first>=0.0){rng.first = z;}
    if(h.first<=0.0){rng.second = z;}
    it += 1;
  }
  Rcpp::Rcout << "sf iteration " << it << std::endl;
  Rcpp::Rcout << "z " << z  << std::endl;
  Rcpp::Rcout << "h " << h.first << " " << h.second << std::endl;
  Rcpp::Rcout << "range " << rng.first << " " << rng.second << std::endl;
  q_sf = q_sf_in + (s_sf - v_sf_rz - z)/Dt;
  s_sf = z;

  
  return states;
}


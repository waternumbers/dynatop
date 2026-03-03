#include "Rcpp.h"
#include <vector>



class hru {
public:
  hru(double s_sf_, double s_rz_, double s_uz_, double s_sz_,
      double area_, double W_beta_,
      std::vector<double> a_, std::vector<double> b_, std::vector<double> s_,
      double s_rz_max_,
      double K_0_, double m_, double D_):
    s_sf(s_sf_), s_rz(s_rz_), s_uz(s_uz_), s_sz(s_sz_),
    area(area_), W_beta(W_beta_),
    a(a_), b(b_), s(s_),
    s_rz_max(s_rz_max_),
    K_0(K_0_), m(m_), D(D_)
  {}

  void update_inputs(double q_sf_in_, double q_sz_in_, double p_, double e_p_){
    q_sf_in = q_sf_in_;
    q_sz_in = q_sz_in_;
    p = p_;
    e_p = e_p;
  }

  std::vector<double> get_states(void){
    std::vector<double> out = {s_sf,s_rz,s_uz,s_sz};
    return out;
  }
  
  
  void iterate(double const &Dt, double const &tol){
    // ///////////////////////
    // Downward Pass   
    v_sf_rz = s_sf + Dt*q_sf_in;
    v_rz_uz = std::min(area*K_0*Dt, std::max(0.0,s_rz + Dt*(p-e_p) + v_sf_rz - s_rz_max));
    
    // solve saturated zone
    double phi = area*K_0*std::exp(-s_sz/m);
    // test 0
    double z = 0;
    std::pair<double,double> h = Hz(z,phi,Dt);
    if( h.first <= 0.0 ){
      q_sz = W_beta*m*K_0;
      v_uz_sz = s_sz + Dt*(q_sz - q_sz_in);
      s_sz = z;
    }else{
      // not saturated
      std::pair<double,double> rng = {0,D};
      //Rcpp::Rcout << "range " << rng.first << " " << rng.second << std::endl;
      z = s_sz;
      h = Hz(z,phi,Dt);
      if(h.first >= 0.0){rng.first = z;}
      if(h.first <= 0.0){rng.second = z;}
      int it(0);
      while( ((h.first>0.0) or (-h.first>tol)) and (it<100) ){
	double zz = z;
	z = z - (h.first/h.second);
	if( std::isnan(z) ){ z = (rng.first + rng.second)/2.0; }
	if( z < rng.first ){ z = (rng.first + zz)/2.0; }
	if( z > rng.second ){ z = (rng.second + zz)/2.0; }
	h = Hz(z,phi,Dt);
	if(h.first>=0.0){rng.first = z;}
	if(h.first<=0.0){rng.second = z;}
	it += 1;
      }
      //Rcpp::Rcout << "z " << z  << std::endl;
      v_uz_sz = Dt*phi*std::min(1.0, (s_uz + v_rz_uz)/(z + Dt*phi));
      q_sz = (z - s_sz + v_uz_sz + q_sz_in) / Dt;
      s_sz = z;
      
      Rcpp::Rcout << "sz iteration " << it << std::endl;
    }
    
    //upward pass
    z = std::min(s_sz,s_uz + v_rz_uz - v_uz_sz);
    v_rz_uz = z - s_uz + v_uz_sz;
    s_uz = z;
    
    v_sf_rz = std::min(v_sf_rz,s_rz_max - s_rz - Dt*(p-e_p) + v_rz_uz);
    z = ( s_rz_max / (s_rz_max + e_p*Dt) ) * (s_rz + Dt*p + v_sf_rz - v_rz_uz);
    e_a = s_rz + v_sf_rz - v_rz_uz + Dt*p - z;
    s_rz = z;
    
    Rcpp::Rcout << "v_sf_rz " << v_sf_rz << std::endl;
    // solve surface
    z = s_sf;
    h = Sw(z,Dt);
    Rcpp::Rcout << "z " << z  << std::endl;
    Rcpp::Rcout << "h " << h.first << " " << h.second << std::endl;
    std::pair<double,double> rng = {0.0,1000*area};
    if(h.first>=0.0){rng.first = z;}
    if(h.second<=0.0){rng.second = z;}
    int it(0);
    while( ((h.first>0) or (-h.first>tol)) and (it<10000) ){
      double zz = z;
      z = z - (h.first/h.second);
      if( std::isnan(z) ){ z = (rng.first + rng.second)/2.0; }
      if( z < rng.first ){ z = (rng.first + zz)/2.0; }
      if( z > rng.second ){ z = (rng.second + zz)/2.0; }
      h = Sw(z,Dt);
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
  }
  // end of iterate

  // inputs
  double s_sf, s_rz, s_uz, s_sz;
  double area, W_beta;
  std::vector<double> a, b, s;
  double s_rz_max;
  double K_0, m, D;

  // external fluxes - initialise due to stupidity...
  double q_sf_in{0.0}, q_sz_in{0.0}, q_sf{0.0}, q_sz{0.0}, e_a{0.0}, p{0.0}, e_p{0.0};
  // internal fluxes
  double v_sf_rz, v_rz_uz, v_uz_sz;
  
private:

  std::pair<double,double> Hz(double &z, double const &phi, double const Dt){
    double omega = std::min(1.0, (s_uz + v_rz_uz)/(z + Dt*phi));
    double q = W_beta*m*K_0*(std::exp(-z/m) - std::exp(-D/m));
    // Rcpp::Rcout << "q z " << q << " " << z  << std::endl;
    std::pair<double,double>h = {s_sz + Dt*(q - q_sz_in - phi*omega) - z,
				     -(q/m) - 1.0};
    if(omega < 1.0){
      h.second += omega*((Dt*phi)/(z+Dt*phi));
    }
    return h;
  }

  std::pair<double,double> Sw(double &w, double const &Dt){
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
    std::pair<double,double> s = {
      s_sf + Dt*(q_sf_in - q) - v_sf_rz - w,
      - Dt*dq - 1
    };
    return s;
  }
};
  

RCPP_MODULE(hru_module) {
  Rcpp::class_<hru>("hru")
    .constructor<double, double, double, double,
		 double, double,
		 std::vector<double>, std::vector<double>, std::vector<double>,
		 double,
		 double, double, double>()
    
    .field("s_sf", &hru::s_sf)
    .field("s_rz", &hru::s_rz)
    .field("s_uz", &hru::s_uz)
    .field("s_sz", &hru::s_sz)
    
    .field("q_sz", &hru::q_sz)
    .field("q_sf", &hru::q_sf)

    .field("v_sf_rz", &hru::v_sf_rz)
    .field("v_rz_uz", &hru::v_rz_uz)
    .field("v_uz_sz", &hru::v_uz_sz)

    .method("update_inputs", &hru::update_inputs)
    .method("iterate", &hru::iterate)
    .method("get_states", &hru::get_states)
    ;
}

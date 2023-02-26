#include "hru.h"


hru::hru(int const id_,
	 std::vector<double> states_,
	 std::vector<double> const properties_,
	 //double const area_, double const s_bar_, double const width_, 
	 int const sf_type_, std::vector<double> const sf_param_,
	 std::vector<double> const rz_param_,
	 std::vector<double> const uz_param_,
	 int const sz_type_, std::vector<double> const sz_param_,
	 std::vector<int> const precip_lnk_id_, std::vector<double> const precip_lnk_frc_,
	 std::vector<int> const pet_lnk_id_, std::vector<double> const pet_lnk_frc_,
	 std::vector<int> const sf_lnk_id_, std::vector<double> const sf_lnk_frc_,
	 std::vector<int> const sz_lnk_id_, std::vector<double> const sz_lnk_frc_
	 ):
  id(id_),
  //states(states_),
  s_sf(states_[0]), s_rz(states_[1]), s_uz(states_[2]), s_sz(states_[3]),
  properties(properties_),
  //q_sf(states_[4]), q_sz(states_[5]), // these states are only used in the flux limiters
  //area(properties_[0]),width(properties_[1]),Dx(properties_[2]), // Dx passed in this way for future flexibility
  sf_param(sf_param_), s_rzmax(rz_param_[0]),
  t_d(uz_param_[0]), sz_param(sz_param_),
  precip_lnk_id(precip_lnk_id_), precip_lnk_frc(precip_lnk_frc_),
  pet_lnk_id(pet_lnk_id_), pet_lnk_frc(pet_lnk_frc_),
  sf_lnk_id(sf_lnk_id_), sf_lnk_frc(sf_lnk_frc_),
  sz_lnk_id(sz_lnk_id_), sz_lnk_frc(sz_lnk_frc_)
{
  // change depths to volues for storage limits
  //s_rzmax = s_rzmax*area;
  area = properties_[0];
  // initialise the surface flux object
  switch(sf_type_){
  case 1:
    // constant celerity with raf
    sf = std::make_unique<sfc_cnst>( sf_param_, properties );
    break;
  case 2:
    // kinematic with raf
    sf = std::make_unique<sfc_kin>( sf_param_, properties );
  }

  // initialise the saturated flux object
  switch(sz_type_){
  case 1:
    //exp
    sz = std::make_unique<szc_exp>( sz_param_, properties ); // properites_[3] is sbar
    break;
  case 2:
    // bounded exp
    sz = std::make_unique<szc_bexp>( sz_param_, properties ); // properites_[3] is sbar
    break;
  case 3:
    // double exp
    sz = std::make_unique<szc_dexp>( sz_param_, properties ); // properites_[3] is sbar
    break;
  case 4:
    // constant celerity
    sz = std::make_unique<szc_dexp>( sz_param_, properties ); // properites_[3] is sbar
    break;
  }
};

void hru::lateral_redistribution(std::vector<double> &vec_q_sf_in,
				 std::vector<double> &vec_q_sz_in){
  for(uint ii=0; ii<sf_lnk_id.size(); ++ii){
    const int &i = sf_lnk_id[ii];
    const double &f = sf_lnk_frc[ii];
    vec_q_sf_in[i] += f * q_sf;
  }
  for(uint ii=0; ii<sz_lnk_id.size(); ++ii){
    const int &i = sz_lnk_id[ii];
    const double &f = sz_lnk_frc[ii];
    vec_q_sz_in[i] += f * q_sz;
  }
}

void hru::update_met(std::vector<double> &obs){
  precip = 0.0;
  for(uint ii=0; ii<precip_lnk_id.size(); ++ii){
    const int &i = precip_lnk_id[ii];
    const double &f = precip_lnk_frc[ii];
    precip += f * area * obs[i];
  }
  pet = 0.0;
  for(uint ii=0; ii<pet_lnk_id.size(); ++ii){
    const int &i = pet_lnk_id[ii];
    const double &f = pet_lnk_frc[ii];
    pet += f * area * obs[i];
  }
}

void hru::init(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double s_rz_0, double r_uz_sz_0,
	       double const &vtol, double const &etol, int const &max_it){

  if(area == 0.0){ // if the HRU has no area just pass on flow and set states to 0
    q_sf_in = vec_q_sf_in[id];
    q_sf = q_sf_in;
    q_sz_in = vec_q_sz_in[id];
    q_sz = q_sz_in; 
    s_sf = std::numeric_limits<double>::quiet_NaN();
    s_rz = std::numeric_limits<double>::quiet_NaN();
    s_uz = std::numeric_limits<double>::quiet_NaN();
    s_sz = std::numeric_limits<double>::quiet_NaN();
    lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    return;
  }

  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id] - q_sz_in;

  // only water at surface if inflow can't be absorbed so max downward flux is
  double r_sf_rz = q_sf_in;
  
  // if steady state then passed downward flux straight to uz
  double r_rz_uz = r_sf_rz;

  // injected water flux into the unsaturated zone
  double r_inj = area * r_uz_sz_0;

  // evaluate flux from unsaturated zone
  double r_uz_sz = std::min( r_rz_uz + r_inj, area/t_d ); // ensure downward flux is possible
  r_uz_sz = std::min( r_uz_sz , sz->q_szmax - q_sz_in ); // revise down so that outflow from saturated zone is possible
  
  // initialise saturated zone to match q_sz
  q_sz = r_uz_sz + q_sz_in;
  double z = (q_sz + q_sz_in)/2.0;
  s_sz = sz->fs(z);

  s_uz = t_d * r_uz_sz * s_sz / area; // compute unsaturated zone storage

  r_rz_uz = r_uz_sz - r_inj;
  if( (r_sf_rz > 0.0) | (r_rz_uz < 0.0)  ){
    s_rz = s_rzmax * area;
  }else{
    s_rz = s_rzmax * s_rz_0 * area;    
  }

  // balance flux through root zone
  r_sf_rz = std::min( r_sf_rz , r_rz_uz );
  
  // solve surface
  q_sf = q_sf_in - r_sf_rz;
  s_sf = sf->fs(q_sf,q_sf_in);
  
  // redistributed the flows
  //lateral_redistribution(vec_q_sf_in,vec_q_sz_in);  
}


void hru::step(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double const &vtol, double const &etol, int const &max_it, double const &Dt)
{
  if(area == 0.0){ // if the HRU has no area just pass on flow
    q_sf_in = vec_q_sf_in[id];
    q_sf = q_sf_in;
    q_sz_in = vec_q_sz_in[id];
    q_sz = q_sz_in;
    lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    return;
  }

  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id]  - q_sz_in;
  
  // single HRU mass balance for development
  std::vector<double> mass_ballance = {s_sf, s_rz, s_uz, s_sz};
  
  // compute first downward flux estimate from surface zone and root zone
  // limites by
  v_sf_rz = s_sf + Dt*q_sf_in;
  
  // change r_rz_uz to present the maximum downward flux
  v_rz_uz = std::max(0.0 ,
		     s_rz - s_rzmax  + Dt*(precip - pet) + v_sf_rz);

  // search for s_sz
  double z = 0.0;
  v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
  double qq = std::min( sz->q_szmax, std::max(0.0, 2*sz->fq(z) - q_sz_in) );
  double Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - qq);
  std::pair<double,double> bnd(0.0, s_sz + 3*vtol);
  int it(0);
  if( Hz <= 0 ){ // then solve by bisection

    // expand
    //std::pair<double,double> bnd(0.0, s_sz + 3*vtol);
    //int it(0);
    bool flg(true);
    while( (flg == true) and (it < max_it) ){
      z = bnd.second;
      v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
      qq = std::min( sz->q_szmax, std::max(0.0, 2*sz->fq(z) - q_sz_in) );
      Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - qq);
      if( Hz<0 ){
	bnd.first = bnd.second;
	bnd.second = 2* bnd.second;
      }else{
	flg = false;
      }
    }
    
    if( flg==true ){
      Rcpp::warning("SZ: No bound found within %i iterations. Difference between bounds is %d.",
		    it, bnd.second - bnd.first);
    }

    it = 0;  
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) ){
      z = (bnd.first+bnd.second)/2.0;
      v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
      qq = std::min( sz->q_szmax, std::max(0.0, 2*sz->fq(z) - q_sz_in) );
      Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - qq);
      if( Hz <= 0 ){ bnd.first= z; } else { bnd.second=z; }
    }
    if(it > max_it){
      Rcpp::warning("SZ: No solution found within %i iterations. Difference between bounds is %d",
		    it, bnd.second - bnd.first);
    }
    z = bnd.second;
  }

  // upward pass
  q_sz = std::min( sz->q_szmax, std::max(0.0, 2*sz->fq(z) - q_sz_in) );
  v_uz_sz = s_sz + Dt*(q_sz-q_sz_in) -z;
  s_sz = z;

  z = std::min(s_sz, s_uz+v_rz_uz-v_uz_sz);
  v_rz_uz = z + v_uz_sz - s_uz;
  s_uz = z;

  v_sf_rz = std::min( v_sf_rz, s_rzmax - s_rz - Dt*(precip - pet) + v_rz_uz);
  s_rz = (s_rzmax / (s_rzmax + Dt*pet)) * (s_rz + Dt*precip + v_sf_rz - v_rz_uz);
  aet = pet * s_rz / (s_rzmax);
  
  // surface
  double sfmax = s_sf + Dt*q_sf_in - v_sf_rz;
  bnd.first = 0.0;
  bnd.second = sfmax;
  it = 0;
  it = 0;  
  while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) ){
    z = (bnd.first+bnd.second)/2.0;
    qq = sf->fq(z,q_sf_in);
    double Sw = sfmax - Dt*qq - z;
    if( Sw <= 0 ){ bnd.second= z; } else { bnd.first=z; }
  }
  z = bnd.first;
  q_sf = sf->fq(z,q_sf_in);
  s_sf = z;
   
  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    
  //Rcpp::Rcout << s_sf << " " << q_sf << " " << q_sf_in << " " << r_sf_rz << " " << Dt << std::endl;
  // // single HRU mass balance for development
  mass_ballance[0] += Dt*(q_sf_in - q_sf) - v_sf_rz - s_sf; // surface
  mass_ballance[1] += Dt*(precip - aet) + v_sf_rz - v_rz_uz - s_rz; // root zone
  mass_ballance[2] += v_rz_uz - v_uz_sz - s_uz; // unsaturated zone
  mass_ballance[3] += Dt*(q_sz - q_sz_in) - v_uz_sz - s_sz; // saturated zone
  z = 0;
  for(uint ii=0; ii<4; ii++){
    z = std::max( z, std::abs(mass_ballance[ii]));
  }
  if( z > 1e-6){
      Rcpp::Rcout << "At end of " << id << std::endl; //": " << mass_ballance << " : " << std::endl;
      Rcpp::Rcout << "     s_sf:  " << mass_ballance[0] << std::endl;
      Rcpp::Rcout << "     s_rz:  " << mass_ballance[1] << std::endl;
      Rcpp::Rcout << "     s_uz:  " << mass_ballance[2] << std::endl;
      Rcpp::Rcout << "     s_sz:  " << mass_ballance[3] << std::endl;
  }
}


// double hru::fsz(double& s, double &qin, double &q_lower, double &q_upper, double const &Dt){ // compute for saturated zone function
//   double r = std::min(  area / t_d, (s_uz + Dt*r_rz_uz)/( (s*t_d) + area*Dt) ); // compute recharge
//   double q = sz->fq(s,qin);
//   q = std::min( std::max(q ,q_lower),q_upper );
//   //Rcpp::Rcout << q << " " << s << " " << r << std:endl;
//   return( s - s_sz + Dt * (qin + r - q) );
// }

// double hru::fsf(double& s, double &qin, double &q_lower, double &q_upper, double const &Dt){ // compute for saturated zone function
//   double q = sf->fq(s,qin);
//   q = std::min( std::max(q ,q_lower),q_upper );
//   //Rcpp::Rcout << q << " " << s << " " << r << std:endl;
//   return( s - s_sf + Dt * (q + r_sf_rz - qin) );
// }

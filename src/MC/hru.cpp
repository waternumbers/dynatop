#include "hru.h"


hru::hru(int const id_,
	 std::vector<double> states_,
	 std::vector<double> const properties_,
	 int const sf_type_, std::vector<double> const sf_param_,
	 std::vector<double> const rz_param_,
	 std::vector<double> const uz_param_,
	 int const sz_type_, std::vector<double> const sz_param_,
	 std::vector<int> const precip_lnk_id_, std::vector<double> const precip_lnk_frc_,
	 std::vector<int> const pet_lnk_id_, std::vector<double> const pet_lnk_frc_,
	 std::vector<int> const sf_lnk_id_, std::vector<double> const sf_lnk_frc_,
	 std::vector<int> const sz_lnk_id_, std::vector<double> const sz_lnk_frc_
	 ):
  //states(states_),
  //properties(properties_),
  rz_param(rz_param_),
  uz_param(uz_param_),
  precip_lnk_id(precip_lnk_id_), precip_lnk_frc(precip_lnk_frc_),
  pet_lnk_id(pet_lnk_id_), pet_lnk_frc(pet_lnk_frc_),
  sf_lnk_id(sf_lnk_id_), sf_lnk_frc(sf_lnk_frc_),
  sz_lnk_id(sz_lnk_id_), sz_lnk_frc(sz_lnk_frc_),
  id(id_),
  s_sf(states_[0]),s_rz(states_[1]),s_uz(states_[2]),s_sz(states_[3]),
  area(properties_[0])
{
  
  // scale states back up from depths to volumes
  s_sf *= area;
  s_rz *= area;
  s_uz *= area;
  s_sz *= area;
  
  
  // initialise the surface flux object
  switch(sf_type_){
  case 1:
    // two stage constant velocity
    sf = std::make_unique<sfc_cnst>( sf_param_, properties_ );
    break;
  case 2:
    // mannings with raf  <TODO> check properties
    sf =  std::make_unique<sfc_kin>( sf_param_, properties_ );
    break;
  }

  // initialise the saturated flux object
  switch(sz_type_){
  case 1:
    //exp
    sz = std::make_unique<szc_bexp>( sz_param_, properties_ );
    break;
    // case 2:
    //   // bounded exp
    //   sz = std::make_unique<szc_bexp>( sz_param_, properties_ );
    //   break;
  case 3:
      // double exp
    sz = std::make_unique<szc_dexp>( sz_param_, properties_ );
    break;
    // case 4:
    //   // constant celerity
    //   sz = std::make_unique<szc_cnst>( sz_param_, properties_ );
    //   break;
  }
};

void hru::lateral_redistribution(std::vector<double> &vec_q_sf_in,
				 std::vector<double> &vec_q_sz_in){
  for(long unsigned int ii=0; ii<sf_lnk_id.size(); ++ii){
    const int &i = sf_lnk_id[ii];
    const double &f = sf_lnk_frc[ii];
    vec_q_sf_in[i] += f * q_sf;
  }
  for(long unsigned int ii=0; ii<sz_lnk_id.size(); ++ii){
    const int &i = sz_lnk_id[ii];
    const double &f = sz_lnk_frc[ii];
    vec_q_sz_in[i] += f * q_sz;
  }
}

void hru::update_met(std::vector<double> &obs){
  precip = 0.0;
  for(long unsigned int ii=0; ii<precip_lnk_id.size(); ++ii){
    const int &i = precip_lnk_id[ii];
    const double &f = precip_lnk_frc[ii];
    precip += f * obs[i];
  }
  precip *= area; // make into flux
  pet = 0.0;
  for(long unsigned int ii=0; ii<pet_lnk_id.size(); ++ii){
    const int &i = pet_lnk_id[ii];
    const double &f = pet_lnk_frc[ii];
    pet += f * obs[i];
  }
  pet *= area; // make into flux
}

void hru::init(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double s_rz_0, double r_uz_sz_0,
	       int const &max_it){

  double const &s_rzmax = rz_param[0];
  double const &t_d = uz_param[0];
  
  //Rcpp::Rcout << id << std::endl;
  //Rcpp::Rcout << sz->q_szmax << std::endl;
  
  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id] - q_sz_in;
  
  //Rcpp::Rcout << "q_sf_in " << q_sf_in << std::endl;
  //Rcpp::Rcout << "q_sz_in " << q_sz_in << std::endl;
  
  // work out maximum downwrd flux when q_sf=0 so downward flux
  // is the same as inflow
  double r_sf_rz = q_sf_in;

  
  // if steady state then passed downward flux straight to uz
  double r_rz_uz = r_sf_rz;
  
  //Rcpp::Rcout << "before r_inj" << std::endl;

  // injected water flux into the unsaturated zone
  double r_inj = area * r_uz_sz_0;
  //Rcpp::Rcout << "area " << area << std::endl;
  //Rcpp::Rcout << "r_uz_sz_0 " << r_uz_sz_0 << std::endl;
  //Rcpp::Rcout << "r_inj " << r_inj << std::endl;
  //Rcpp::Rcout << "before unsat" << std::endl;
  
  // evaluate flux from unsaturated zone
  double r_uz_sz = std::min( r_rz_uz + r_inj, area/t_d ); // ensure downward flux is possible

  //Rcpp::Rcout << "r_uz_sz " << r_uz_sz << std::endl;
  
  
  //Rcpp::Rcout << "before after unsat" << std::endl;
  
  // work out outflow
  q_sz = std::min(sz->q_szmax, r_uz_sz + q_sz_in);
  // update the saturated zone reference flow
  double q_ref = q_sz; //(q_sz + q_sz_in)/2.0;
  //Rcpp::Rcout << "q_sz " << q_sz << std::endl;
  //Rcpp::Rcout << "q_ref " << q_ref << std::endl;
  sz->update(q_ref);
  double& k = sz->kappa;
  double& eta = sz->eta;
  //Rcpp::Rcout << "kappa " << k << std::endl;
  //Rcpp::Rcout << "eta " << eta << std::endl;
  //Rcpp::Rcout << "sz_max " << sz->s_szmax << std::endl;
  s_sz = sz->s_szmax - (k*( eta*q_sz_in + (1-eta)*q_sz ));
  if(s_sz < 0.0){
    Rcpp::Rcout << "id " << id << std::endl;
    Rcpp::Rcout << "s_sz " << s_sz << std::endl;
    Rcpp::Rcout << "q_sz " << q_sz << std::endl;
    Rcpp::Rcout << "q_sz_in " << q_sz_in << std::endl;
    Rcpp::Rcout << "kappa " << k << std::endl;
    Rcpp::Rcout << "eta " << eta << std::endl;
    Rcpp::Rcout << "storage " << (k*( eta*q_sz_in + (1-eta)*q_sz )) << std::endl;
    Rcpp::Rcout << "sz_max " << sz->s_szmax << std::endl;
  }
  
  //Rcpp::Rcout << "start upward" << std::endl;
  
  r_uz_sz = q_sz - q_sz_in;
    
  s_uz = t_d * r_uz_sz * (s_sz/area) ; //sz->h; // compute unsaturated zone storage
  if( s_uz > s_sz ){
    Rcpp::Rcout << id << " unsaturated" << std::endl;
    Rcpp::Rcout << s_sz << " " << s_uz << " " << r_uz_sz << std::endl;
    Rcpp::Rcout << q_sz << " " << q_sz_in << " " << q_ref << std::endl;
    Rcpp::Rcout << s_uz - s_sz << std::endl;
  }
  
  r_rz_uz = r_uz_sz - r_inj;
  if( (r_sf_rz > 0.0) | (r_rz_uz < 0.0)  ){
    s_rz = s_rzmax * area;
  }else{
    s_rz = s_rzmax * s_rz_0 * area;    
  }

  // balance flux through root zone
  r_sf_rz = std::min( r_sf_rz , r_rz_uz );

  //Rcpp::Rcout << "r_sf_rz upward " << r_sf_rz << std::endl;
  // solve surface
  q_sf = q_sf_in - r_sf_rz;
  //Rcpp::Rcout << "q_sf " << q_sf << std::endl;
  q_ref = (q_sf + q_sf_in) / 2.0;
  //Rcpp::Rcout << "q_ref " << q_ref << std::endl;
  sf->update(q_ref);
  k = sf->kappa;
  eta = sf->eta;
  //Rcpp::Rcout << "kappa " << k << std::endl;
  //Rcpp::Rcout << "eta " << eta << std::endl;
  s_sf = k*( eta*q_sf_in + (1-eta)*q_sf );
  //Rcpp::Rcout << "s_sf " << s_sf << std::endl;
  // redistributed the flows
  //Rcpp::Rcout << "q_sf " << q_sf << std::endl;
  //Rcpp::Rcout << "q_sz " << q_sz << std::endl;
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);

  // debug printing
  double tmp = q_sz_in + q_sf_in + r_inj - q_sz - q_sf;
  if( std::abs(tmp) > 1e-10 ){
    Rcpp::Rcout << "MASS BALANCE" << std::endl;
    Rcpp::Rcout << id << std::endl;
    Rcpp::Rcout << q_sf_in << " " << q_sf << std::endl;
    Rcpp::Rcout << q_sz_in << " " << q_sz << " " << r_inj << std::endl;
    Rcpp::Rcout << tmp << std::endl;
  }
  
}


void hru::step(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       int const &max_it, double const &Dt)
{

  double const &s_rzmax = rz_param[0];
  double const &t_d = uz_param[0];
  
  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id]  - q_sz_in;
 
  if( q_sf_in < 0.0 ){
    Rcpp::Rcout << "q_sf_in error " << q_sf_in << " " << id << std::endl;
  }
  if( q_sz_in < 0.0 ){
    Rcpp::Rcout << "q_sz_in error " << q_sz_in << " " << id << std::endl;
  }

  // single HRU mass balance for development
  std::vector<double> mass_ballance = {s_sf, s_rz, s_uz, s_sz};
  
  // compute first downward flux estimate from surface zone and root zone
  // limites by outflow of 0
  double q_ref = q_sf_in / 2.0;
  sf->update(q_ref);
  v_sf_rz = std::max(0.0, s_sf + Dt*q_sf_in - (sf->kappa*sf->eta*q_sf_in)) ;

  // change r_rz_uz to present the maximum downward flux
  v_rz_uz = std::max(0.0 ,
		     s_rz - (area*s_rzmax)  + Dt*(precip - pet) + v_sf_rz);

  v_uz_sz = area* Dt * std::min( (s_uz+v_rz_uz)/(t_d*s_sz + Dt), 1/t_d ); // could put back into the loop
  // solve for saturated zone
  double q_out = q_sz_in;
  for(long unsigned int ii=0; ii<max_it; ++ii){
    q_ref = q_out; //(q_sz_in + q_out) /2.0;
    sz->update(q_ref);
    //double& h = sz->h; // storage comparamble to q_ref
    double& kappa = sz->kappa;
    double& eta = sz->eta;
    q_out = std::min(sz->q_szmax, std::max(0.0, ( sz->s_szmax - s_sz + v_uz_sz + (Dt-kappa*eta)*q_sz_in ) / ( Dt + kappa*(1.0-eta) ) ));
  }
  double pjs2 = v_uz_sz;
  q_sz = q_out;
  double z = s_sz + Dt*(q_sz - q_sz_in); // max storage
  s_sz = std::max(0.0, z - v_uz_sz);
  v_uz_sz = z - s_sz;
  
  double pjs = mass_ballance[3] + Dt*(q_sz - q_sz_in) - v_uz_sz - s_sz;
  if( std::abs(pjs) > 1e-10 ){
    Rcpp::Rcout << "id = " << id << std::endl;
    Rcpp::Rcout << "error = " << pjs << std::endl;
    Rcpp::Rcout << "Initial s_sz = " << mass_ballance[3] << " " << sz->h << " " << t_d << std::endl;
    Rcpp::Rcout << "Initial v_uz_sz = " << pjs2 << std::endl;
    Rcpp::Rcout << "q_sz = " << q_sz << std::endl;
    Rcpp::Rcout << "q_sz_in = " << q_sz_in << std::endl;
    Rcpp::Rcout << "s_uz = " << s_uz << std::endl;
    Rcpp::Rcout << "s_sz = " << s_sz << std::endl;
    Rcpp::Rcout << "v_rz_uz = " << v_rz_uz << std::endl;
    Rcpp::Rcout << "v_uz_sz = " << v_uz_sz << std::endl;
  }
  if( s_sz > sz->s_szmax ){
    Rcpp::Rcout << "id = " << id << std::endl;
    Rcpp::Rcout << "error = " << s_sz - sz->s_szmax << std::endl;
  }
  if( s_sz <0.0 ){
    Rcpp::Rcout << "id = " << id << std::endl;
    Rcpp::Rcout << "error = " << s_sz  << std::endl;
  }
  // solve unsaturated zone
  z = std::max(0.0, std::min(s_sz, s_uz+v_rz_uz-v_uz_sz));
  v_rz_uz = z + v_uz_sz - s_uz;
  s_uz = z;

  // solve root zone
  v_sf_rz = std::min( v_sf_rz, (area*s_rzmax) - s_rz - Dt*(precip - pet) + v_rz_uz);
  s_rz = ((area*s_rzmax) / ((area*s_rzmax) + Dt*pet)) * (s_rz + Dt*precip + v_sf_rz - v_rz_uz);
  aet = pet * s_rz / (area*s_rzmax);
  
  // surface
  q_out = q_sf_in;
  for(long unsigned int ii=0; ii<max_it; ++ii){
    q_ref = q_out;//(q_sf_in + q_out) /2.0;
    sf->update(q_ref);
    double& kappa = sf->kappa;
    double& eta = sf->eta;
    q_out = std::max(0.0, ( s_sf - v_sf_rz + (Dt - kappa*eta)*q_sf_in ) / ( Dt + kappa*(1.0-eta) ) );
  }
  q_sf = q_out;
  s_sf = s_sf - v_sf_rz + Dt*(q_sf_in-q_sf);
  // Rcpp::Rcout << "q_sf " << q_sf << std::endl;
  // Rcpp::Rcout << "s_sf " << s_sf << std::endl;
  //  Rcpp::Rcout << "v_sf_rz " << s_sf << std::endl;
    
  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    
  //Rcpp::Rcout << s_sf << " " << q_sf << " " << q_sf_in << " " << r_sf_rz << " " << Dt << std::endl;
  // // single HRU mass balance for development
  mass_ballance[0] += Dt*(q_sf_in - q_sf) - v_sf_rz - s_sf; // surface
  mass_ballance[1] += Dt*(precip - aet) + v_sf_rz - v_rz_uz - s_rz; // root zone
  mass_ballance[2] += v_rz_uz - v_uz_sz - s_uz; // unsaturated zone
  mass_ballance[3] += Dt*(q_sz - q_sz_in) - v_uz_sz - s_sz; // saturated zone
  z = 0;
  for(int ii=0; ii<4; ii++){
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

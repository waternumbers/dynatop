#include "hru.h"


hru::hru(int const id_,
	 std::vector<double> states_,
	 std::vector<double> const properties_,
	 int const sf_type_, std::vector<double> const sf_param_,
	 std::vector<double> const rz_param_,
	 std::vector<double> const uz_param_,
	 int const sz_type_, std::vector<double> const sz_param_,
	 std::vector<int> const precip_lnk_id_, std::vector<double> const precip_lnk_area_,
	 std::vector<int> const pet_lnk_id_, std::vector<double> const pet_lnk_area_,
	 std::vector<int> const sf_lnk_id_, std::vector<double> const sf_lnk_frc_,
	 std::vector<int> const sz_lnk_id_, std::vector<double> const sz_lnk_frc_,
	 std::vector<double> const sz_lnk_gradient_//,
	 ):
  s_rzmax(rz_param_[0]),
  t_d(uz_param_[0]),
  precip_lnk_id(precip_lnk_id_), precip_lnk_area(precip_lnk_area_),
  pet_lnk_id(pet_lnk_id_), pet_lnk_area(pet_lnk_area_),
  sf_lnk_id(sf_lnk_id_), sf_lnk_frc(sf_lnk_frc_),
  sz_lnk_id(sz_lnk_id_), sz_lnk_frc(sz_lnk_frc_),
  id(id_), 
  s_sf(states_[0]), s_rz(states_[1]), s_uz(states_[2]), s_sz(states_[3])
{
  // change depths to volues for storage limits - not use hru area no map area
  // area = properties_[0] * properties_[1]; // width * Dx
  //s_rzmax = s_rzmax*area;
  s_sf *= area;
  s_rz *= area;
  s_uz *= area;
  s_sz *= area;
  
  
  // initialise the surface flux object
  for(long unsigned int ii=0; ii<sf_lnk_id.size(); ++ii){
    switch(sf_type_){
    case 1:
      // two stage constant velocity
      sf.push_back( std::make_unique<sfc_cnst>( sf_param_, properties_[2] ) );
      break;
    case 2:
      // mannings with raf  <TODO> check properties
      sf.push_back( std::make_unique<sfc_kin>( sf_param_, properties_[1], properties_[3], properties_[2] ) );
      break;
    }
  }

  // initialise the saturated flux object
  for(long unsigned int ii=0; ii<sz_lnk_id.size(); ++ii){
    switch(sz_type_){
    case 1:
      //exp
      sz.push_back( std::make_unique<szc_exp>( sz_param_, sz_lnk_width_[ii], sz_lnk_gradient_[ii], area ) );
      break;
    case 2:
      // bounded exp
      sz.push_back( std::make_unique<szc_bexp>( sz_param_, sz_lnk_width_[ii], sz_lnk_gradient_[ii], area ) );
      break;
    case 3:
      // double exp
      sz.push_back( std::make_unique<szc_dexp>( sz_param_, sz_lnk_width_[ii], sz_lnk_gradient_[ii], area ) );
      break;
    case 4:
      // constant celerity
      sz.push_back( std::make_unique<szc_cnst>( sz_param_, sz_lnk_width_[ii], sz_lnk_gradient_[ii], area ) );
      break;
    }
  }
};

void hru::lateral_redistribution(std::vector<double> &vec_q_sf_in,
				 std::vector<double> &vec_q_sz_in){
  for(long unsigned int ii=0; ii<sf_lnk_id.size(); ++ii){
    const int &i = sf_lnk_id[ii];
    const double &f = lambda_sf[ii];
    vec_q_sf_in[i] += f * s_sf;
  }
  for(long unsigned int ii=0; ii<sz_lnk_id.size(); ++ii){
    const int &i = sz_lnk_id[ii];
    const double &f = lambda_sz[ii];
    vec_q_sz_in[i] += f * s_sz;
  }
}

void hru::update_met(std::vector<double> &obs){
  precip = 0.0;
  for(long unsigned int ii=0; ii<precip_lnk_id.size(); ++ii){
    const int &i = precip_lnk_id[ii];
    const double &f = precip_lnk_area[ii];
    precip += f * obs[i];
  }
  pet = 0.0;
  for(long unsigned int ii=0; ii<pet_lnk_id.size(); ++ii){
    const int &i = pet_lnk_id[ii];
    const double &f = pet_lnk_area[ii];
    pet += f * obs[i];
  }
}

void hru::init(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double s_rz_0, double r_uz_sz_0,
	       double const &vtol, double const &etol, int const &max_it){

  if(map_area == 0.0){
    // if the HRU has no area then presume it is a channel, return all flow to the surface and set to saturated
    q_sf_in = vec_q_sf_in[id];
    double r_sf_rz = - vec_q_sz_in[id]; // add return flow from the saturated zone;
    q_sf = q_sf_in - r_sf_rz;
    s_sf = sf->fs(q_sf_in,r_sf_rz);
    
    if( std::abs( sf->fq(s_sf,q_sf_in,r_sf_rz) - q_sf ) > 1e-10 ){
      Rcpp::Rcout << id << " surface" << std::endl;
      Rcpp::Rcout << q_sf_in << " " << q_sf << std::endl;
      Rcpp::Rcout << s_sf << " " << sf->fq(s_sf,q_sf_in,r_sf_rz) << std::endl;
    }

    // set the other states
    s_rz = s_rzmax * area;
    s_uz = 0.0;
    s_sz = 0.0;
    q_sz = 0.0;
    
    // redistributed the flows
    lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    return;
  }

  // redivide the inflow so q_sz_in is less then q_szmax
  //q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  //q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id] - q_sz_in;
  q_sz_in = vec_q_sz_in[id];
  q_sf_in = vec_q_sf_in[id];

  // only water at surface if inflow can't be absorbed so max downward flux is
  double r_sf_rz = q_sf_in;
  
  // if steady state then passed downward flux straight to uz
  double r_rz_uz = r_sf_rz;

  // injected water flux into the unsaturated zone
  double r_inj = map_area * r_uz_sz_0;

  
  // evaluate flux from unsaturated zone
  double r_uz_sz = std::min( r_rz_uz + r_inj, area/t_d ); // ensure downward flux is possible

  //  double q_sz_max = sz->fq(0,q_sz_in);
  // r_uz_sz = std::min( r_uz_sz , q_sz_max - q_sz_in );

  //r_uz_sz = std::min( r_uz_sz , sz->q_szmax - q_sz_in ); // revise down so that outflow from saturated zone is possible
  
  // initialise saturated zone to match q_sz
  //q_sz = r_uz_sz + q_sz_in;
  //s_sz = sz->fs(q_sz,q_sz_in);
  // make initial estimate of outflow
  q_sz = r_uz_sz + q_sz_in;
  if(id == 25684){
    Rcpp::Rcout << r_uz_sz << " " << q_sz_in << " " << q_sz << std::endl;
  }

  s_sz = sz->fs(q_sz,q_sz_in);
  q_sz = sz->fq(s_sz,q_sz_in);
  r_uz_sz = q_sz - q_sz_in;
  if(id == 25684){
    Rcpp::Rcout << r_uz_sz << " " << s_sz << " " << q_sz << std::endl;
  }
  
  if( std::abs( sz->fq(s_sz,q_sz_in) - q_sz ) > 1e-10 ){
    Rcpp::Rcout << id << " saturated" << std::endl;
    Rcpp::Rcout << q_sz_in << " " << q_sz << std::endl; 
    Rcpp::Rcout << s_sz << " " << sz->fq(s_sz,q_sf_in) << std::endl;
  }
  
  if(id == 25684){
    Rcpp::Rcout << t_d << " " << r_uz_sz << " " << s_sz << " " << area << std::endl;
  }
  
  s_uz = t_d * r_uz_sz * s_sz / area; // compute unsaturated zone storage
  if( s_uz > s_sz ){
    Rcpp::Rcout << id << " unsaturated" << std::endl;
    Rcpp::Rcout << s_sz << " " << s_uz << " " << r_uz_sz << std::endl;
    Rcpp::Rcout << q_sz << " " << q_sz_in << " " << sz->fq(0.0,q_sz_in) << std::endl;
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
  
  // solve surface
  q_sf = q_sf_in - r_sf_rz;
  //Rcpp::Rcout << "Initialising surface " << id << " " << q_sf << " " << q_sf_in << std::endl;
  //s_sf = sf->fs(q_sf,q_sf_in);
  s_sf = sf->fs(q_sf_in,r_sf_rz);
  if( std::abs( sf->fq(s_sf,q_sf_in,r_sf_rz) - q_sf ) > 1e-10 ){
    Rcpp::Rcout << id << " surface" << std::endl;
    Rcpp::Rcout << q_sf_in << " " << q_sf << std::endl;
    Rcpp::Rcout << s_sf << " " << sf->fq(s_sf,q_sf_in,r_sf_rz) << std::endl;
  }
  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);

  // debug printing
  double tmp = q_sz_in + q_sf_in + r_inj - q_sz - q_sf;
  if( std::abs(tmp) > 1e-10 ){
    Rcpp::Rcout << id << std::endl;
    Rcpp::Rcout << q_sf_in << " " << q_sf << std::endl;
    Rcpp::Rcout << q_sz_in << " " << q_sz << " " << r_inj << std::endl;
    Rcpp::Rcout << tmp << std::endl;
  }
  
}


void hru::step(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double const &vtol, double const &etol, int const &max_it, double const &Dt)
{

  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id]  - q_sz_in;

  // single HRU mass balance for development
  std::vector<double> mass_ballance = {s_sf, s_rz, s_uz, s_sz};
  
  // compute first downward flux estimate from surface zone and root zone
  // limites by outflow of 0
  double q_ref = q_sf_in / 2.0;
  sf->update(q_ref);
  v_sf_rz = s_sf + Dt*q_sf_in - sf.kappa*sf.eta*q_sf_in ;
  
  // change r_rz_uz to present the maximum downward flux
  v_rz_uz = std::max(0.0 ,
		     s_rz - (area*s_rzmax)  + Dt*(precip - pet) + v_sf_rz);

  // solve for saturated zone
  double q_in = sz->q_szmax - q_sz_in;
  double q_out = q_in;
  for(long unsigned int ii=0; ii<max_it; ++ii){
    q_ref = (q_in + q_out) /2.0;
    sz->update(q_ref);
    double& h = sz->s; // storage comparamble to q_ref
    double& kappa = sz->kappa;
    double& eta = sz->eta;
    v_uz_sz = Dt * std::min( (s_uz+v_rz_uz)/(t_d*h + Dt), area/t_d );
    q_out = std::max(0.0, ( s_sz - v_uz_sz + (Dt-kappa*eta)*q_in ) / ( Dt + kappa*(1.0-eta) ) );
  }
  s_max = s_sz + Dt*(q_in - q_out); // max storage
  s_sz = std::max(0.0, s_max - v_uz_sz);
  v_uz_sz = s_max - s_sz;
  q_sf = sz->q_szmax - q_out;

  // solve unsaturated zone
  z = std::min(s_sz, s_uz+v_rz_uz-v_uz_sz);
  v_rz_uz = z + v_uz_sz - s_uz;
  s_uz = z;
  // solve root zone
  v_sf_rz = std::min( v_sf_rz, (area*s_rzmax) - s_rz - Dt*(precip - pet) + v_rz_uz);
  s_rz = ((area*s_rzmax) / ((area*s_rzmax) + Dt*pet)) * (s_rz + Dt*precip + v_sf_rz - v_rz_uz);
  aet = pet * s_rz / (area*s_rzmax);
  
  // surface
  q_out = q_sf_in;
  for(long unsigned int ii=0; ii<max_it; ++ii){
    q_ref = (q_sf_in + q_out) /2.0;
    sf->update(q_ref);
    double& kappa = sz->kappa;
    double& eta = sz->eta;
    q_out = std::max(0.0, ( s_sf - v_sf_rz + (Dt-kappa*eta)*q_in ) / ( Dt + kappa*(1.0-eta) ) );
  }
  q_sf = q_out;
  s_sf = s_sf - v_sf_rz + Dt*(q_in-q_out);
     
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

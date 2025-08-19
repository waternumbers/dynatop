#include "hru.h"


hru::hru(int const id_, int const band_,
	 std::vector<double> states_,
	 std::vector<double> const properties_,
	 int const sf_type_, std::vector<double> const sf_param_,
	 int const rz_type_, std::vector<double> const rz_param_,
	 int const uz_type_, std::vector<double> const uz_param_,
	 int const sz_type_, std::vector<double> const sz_param_,
	 std::vector<int> const precip_lnk_id_, std::vector<double> const precip_lnk_frc_,
	 std::vector<int> const pet_lnk_id_, std::vector<double> const pet_lnk_frc_,
	 std::vector<int> const sf_lnk_id_, std::vector<double> const sf_lnk_frc_,
	 std::vector<int> const sz_lnk_id_, std::vector<double> const sz_lnk_frc_,
	 std::vector<double> const initial_values_,
	 std::vector<double> &vec_q_sf_in_,
	 std::vector<double> &vec_q_sz_in_,
	 double const &vtol_,
	 double const &etol_,
	 int const &max_it_,
	 double const &Dt_
	 ):
  //states(states_),
  //properties(properties_),
  //q_sf(states_[4]), q_sz(states_[5]), // these states are only used in the flux limiters
  //area(properties_[0]),width(properties_[1]),Dx(properties_[2]), // Dx passed in this way for future flexibility
  sf_param(sf_param_), rz_param(rz_param_),
  uz_param(uz_param_), sz_param(sz_param_),
  precip_lnk_id(precip_lnk_id_), precip_lnk_frc(precip_lnk_frc_),
  pet_lnk_id(pet_lnk_id_), pet_lnk_frc(pet_lnk_frc_),
  sf_lnk_id(sf_lnk_id_), sf_lnk_frc(sf_lnk_frc_),
  sz_lnk_id(sz_lnk_id_), sz_lnk_frc(sz_lnk_frc_),
  initial_values(initial_values_),
  vec_q_sf_in(vec_q_sf_in_), vec_q_sz_in(vec_q_sz_in_),
  vtol(vtol_), etol(etol_), max_it(max_it_),
  Dt(Dt_),
  id(id_),
  band(band_),
  s_sf(states_[0]), s_rz(states_[1]), s_uz(states_[2]), s_sz(states_[3])
  {
    // change depths to volues for storage limits - use hru area not map area
    area = properties_[0];
    s_sf *= area;
    s_rz *= area;
    s_uz *= area;
    s_sz *= area;
    
    Dx = properties_[1];
    
    // initialise the surface flux object
  switch(sf_type_){
  case 1:
    // constant celerity with raf
    sf = std::make_unique<sfc_cnst>( sf_param_, properties_ );
    break;
  case 2:
    // kinematic with raf
    sf = std::make_unique<sfc_kin>( sf_param_, properties_ );
    break;
  case 3:
    // compound channel
    sf = std::make_unique<sfc_comp>( sf_param_, properties_ );
    break;
  case 4:
    // manning with raf solved as tank
    sf = std::make_unique<sfc_power_law>( sf_param_, properties_ );
    break;
  case 5:
    // arbitary area flow relationship
    sf = std::make_unique<sfc_arb_kin>( sf_param_, properties_ );
    break;
  case 6:
    // Muskingham-Cunge-Todini
    sf = std::make_unique<sfc_mct>( sf_param_, properties_ );
    break;
  case 7:
    sf = std::make_unique<sfc_mct_rect>( sf_param_, properties_ );
    break;
  }

  // initialise the saturated flux object
  switch(sz_type_){
  case 1:
    //exp
    sz = std::make_unique<szc_exp>( sz_param_, properties_ ); // properites_[3] is sbar
    break;
  case 2:
    // bounded exp
    sz = std::make_unique<szc_bexp>( sz_param_, properties_ ); // properites_[3] is sbar
    break;
  case 3:
    // double exp
    sz = std::make_unique<szc_dexp>( sz_param_, properties_ ); // properites_[3] is sbar
    break;
  case 4:
    // constant celerity
    sz = std::make_unique<szc_cnst>( sz_param_, properties_ ); // properites_[3] is sbar
    break;
  }
};

void hru::lateral_redistribution(){
  // std::vector<double> &vec_q_sf_in,
  // 				 std::vector<double> &vec_q_sz_in){
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
    precip += f * area * obs[i];
  }
  pet = 0.0;
  for(long unsigned int ii=0; ii<pet_lnk_id.size(); ++ii){
    const int &i = pet_lnk_id[ii];
    const double &f = pet_lnk_frc[ii];
    pet += f * area * obs[i];
  }
}

void hru::init(){
  // std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
  // 	       double s_rz_0, double r_uz_sz_0,
  // 	       double const &vtol, double const &etol, int const &max_it){

  double const &s_rzmax = rz_param[0];
  double const &t_d = uz_param[0];

  double const &s_rz_0 = initial_values[0];
  double const &r_uz_sz_0 = initial_values[1];
  
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id];
  q_sz_in = std::min( sz->q_szmax, vec_q_sz_in[id]);
  q_sf_in -= q_sz_in;
  
  // only water at surface if inflow can't be absorbed so max downward flux is
  double r_sf_rz = q_sf_in;
  
  // if steady state then passed downward flux straight to uz
  double r_rz_uz = r_sf_rz;

  // injected water flux into the unsaturated zone
  double r_inj = area * r_uz_sz_0; // use map area since input is standardised by map area
  
  // evaluate flux from unsaturated zone
  double r_uz_sz = std::min( r_rz_uz + r_inj, area/t_d ); // ensure downward flux is possible

  // make initial estimate of outflow
  q_sz = std::min( sz->q_szmax, r_uz_sz + q_sz_in );
  // muskingham cunge s_sz = sz->fs( (q_sz+q_sz_in)/2.0 );
  s_sz = sz->fs(q_sz);

  r_uz_sz = q_sz - q_sz_in;
  //if( std::abs( sz->fq(s_sz,q_sz_in) - q_sz ) > 1e-10 ){
  // if( std::abs( sz->fq(s_sz) - q_sz ) > 1e-10 ){
  //   Rcpp::Rcout << id << " saturated" << std::endl;
  //   Rcpp::Rcout << q_sz_in << " " << q_sz << std::endl; 
  //   Rcpp::Rcout << s_sz << " " << sz->fq(s_sz) << std::endl; //,q_sf_in) << std::endl;
  // }
  
  s_uz = t_d * r_uz_sz * s_sz / area; // compute unsaturated zone storage
  // if( s_uz > s_sz ){
  //   Rcpp::Rcout << id << " unsaturated" << std::endl;
  //   Rcpp::Rcout << s_sz << " " << s_uz << " " << r_uz_sz << std::endl;
  //   //Rcpp::Rcout << q_sz << " " << q_sz_in << " " << sz->fq(0.0) << std::endl; //,q_sz_in) << std::endl;
  //   Rcpp::Rcout << s_uz - s_sz << std::endl;
  // }
  

  
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
  s_sf = sf->fS(q_sf);
  // s_sf = sf->fs(q_sf_in,r_sf_ );
  // if( std::abs( sf->fq(s_sf) - q_sf ) > 1e-10 ){ //,q_sf_in,r_sf_rz) - q_sf ) > 1e-10 ){
  //   Rcpp::Rcout << id << " surface" << std::endl;
  //   Rcpp::Rcout << q_sf_in << " " << q_sf << std::endl;
  //   Rcpp::Rcout << s_sf << " " << sf->fq(s_sf) << std::endl; //,q_sf_in,r_sf_rz) << std::endl;
  // }
  // redistributed the flows
  //lateral_redistribution(); //vec_q_sf_in,vec_q_sz_in);

  // debug printing
  double tmp = q_sz_in + q_sf_in + r_inj - q_sz - q_sf;
  // if( std::abs(tmp) > 1e-10 ){
  //   Rcpp::Rcout << "Intiailisation Error: " << id << std::endl;
  //   Rcpp::Rcout << q_sf_in << " " << q_sf << std::endl;
  //   Rcpp::Rcout << q_sz_in << " " << q_sz << " " << r_inj << std::endl;
  //   Rcpp::Rcout << tmp << std::endl;
  // }
  
}


void hru::step(){
  // std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
  // 	       double const &vtol, double const &etol, int const &max_it, double const &Dt)
  // {

  double const &s_rzmax = rz_param[0];
  double const &t_d = uz_param[0];

  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id];
  q_sz_in = std::min( sz->q_szmax, vec_q_sz_in[id]);
  q_sf_in -= q_sz_in;

  // single HRU mass balance for development
  std::vector<double> mass_ballance = {s_sf, s_rz, s_uz, s_sz};
  
  // compute first downward flux estimate from surface zone and root zone
  // limites by
  v_sf_rz = s_sf + Dt*q_sf_in;
  
  // change r_rz_uz to present the maximum downward flux
  v_rz_uz = std::max(0.0 ,
		     s_rz - (area*s_rzmax)  + Dt*(precip - pet) + v_sf_rz);

  // update the saturated zone
  // semi-implicit
  v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*s_sz + area*Dt), 1/t_d );
  double eta(0.0);
  if( s_sz > 0 ){
    eta = (sz->q_szmax - q_sz) / s_sz;
  }
  q_sz = (sz->q_szmax - eta*std::max(0.0,s_sz - v_uz_sz - Dt*q_sz_in)) / (1.0 + eta*Dt);
  double z = std::max(0.0,s_sz + Dt*(q_sz - q_sz_in) - v_uz_sz);
  
  
  // Musk-Cunge solution
  // v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*s_sz + area*Dt), 1/t_d );
  // double Qref = sz->fq( s_sz );
  // double vel = (sz->q_szmax - Qref)  / (s_sz/Dx);
  // double eta = Dx/(2*vel);
  // q_sz = ( (2*eta*sz->q_szmax) + (Dt-eta)*q_sz_in - s_sz + v_uz_sz ) / (Dt + eta);
  // q_sz = std::max( 0.0, std::min(sz->q_szmax,q_sz) );
  // double z = std::max(0.0, s_sz + Dt*(q_sz - q_sz_in) - v_uz_sz);

  // iterative MCT
  // std::pair<double,double> ubnd(0.0, 9999.9); // wettest saturated zone
  // double z = ubnd.first;
  // double Qref = sz->fq( z ); // reference flow
  // q_sz = std::min(sz->q_szmax, std::max(0.0,2.0*Qref - q_sz_in)); ///2.0;
  // v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
  // ubnd.second = z - s_sz + Dt*(q_sz_in - q_sz) + v_uz_sz; // should be -ve
  
  				  
  // upward pass
  v_uz_sz = s_sz + Dt*(q_sz-q_sz_in) - z;
  s_sz = z;
  
  z = std::max(0.0,std::min(s_sz, s_uz+v_rz_uz-v_uz_sz)); // to stop negative values appearing
  v_rz_uz = z + v_uz_sz - s_uz;
  s_uz = z;

  v_sf_rz = std::min( v_sf_rz, (area*s_rzmax) - s_rz - Dt*(precip - pet) + v_rz_uz);
  s_rz = ((area*s_rzmax) / ((area*s_rzmax) + Dt*pet)) * (s_rz + Dt*precip + v_sf_rz - v_rz_uz);
  aet = pet * s_rz / (area*s_rzmax);
  
  // surface
  z = s_sf + (Dt*q_sf_in) - v_sf_rz; // max surface storage
  if( z == 0.0 ){ // no stroage so no outflow
    s_sf = 0;
    q_sf = 0.0;
  }else{
    // semi-explicit solution
    double T_sf = sf->fT(z); // q_sf = T_sf * s_sf
    s_sf = z / (1 + Dt*T_sf) ;
    q_sf = (z - s_sf)/Dt;
    if( std::isnan(s_sf) | std::isnan(q_sf) | (s_sf<0) | (q_sf<0) ){
      Rcpp::Rcout << "id: " << id << " T_sf: " << T_sf << " s_sf: " << s_sf << " q_sf_in " << q_sf_in << " v_sf_rz: " << v_sf_rz << " z: " << z << std::endl;
      Rcpp::Rcout << "s_sf: " << s_sf << " q_sf: " << q_sf << std::endl;
    }
   
    // Full numerical search solution
    // std::pair<double,double> ubnd(z/Dt, 9999.9); // wettest surface
    // q_sf= ubnd.first;
    // Qref = (q_sf+q_sf_in)/2.0;
    // sf->update( Qref );
    // ubnd.second = q_sf - std::max(0.0, (z - (sf->kappa*sf->eta)*q_sf_in) / (Dt + sf->kappa*(1.0 - sf->eta)) ); // negative
    
    // std::pair<double,double> lbnd(0, 9999.9); // driest surface
    // q_sf = lbnd.first;
    // Qref = (q_sf+q_sf_in)/2.0;
    // sf->update( Qref );
    // lbnd.second = q_sf - std::max(0.0, (z - (sf->kappa*sf->eta)*q_sf_in) / (Dt + sf->kappa*(1.0 - sf->eta)) );

    // // Rcpp::Rcout << "id is " << id << std::endl;
    // // Rcpp::Rcout << "upper bound " << ubnd.first << " " << ubnd.second << std::endl;
    // // Rcpp::Rcout << "lower bound " << lbnd.first << " " << lbnd.second << std::endl;
    // int it = 0;
    // while( (lbnd.second < -1e-6) & (it < 100) ){
    //   //for(int it=0; it<max_it; ++it){
    //   q_sf = (ubnd.first + lbnd.first)/ 2.0;
    //   Qref = (q_sf+q_sf_in)/2.0;
    //   sf->update( Qref );
    //   double e = q_sf - std::max(0.0, (z - (sf->kappa*sf->eta)*q_sf_in) / (Dt + sf->kappa*(1.0 - sf->eta)) );
    //   if( e <= 0.0 ){
    // 	lbnd.first = q_sf;
    // 	lbnd.second = e;
    //   }else{
    // 	ubnd.first = q_sf;
    // 	ubnd.second = e;
    //   }
    //   it +=1;
    // }
    // q_sf = lbnd.first;
    // Qref = (q_sf+q_sf_in)/2.0;
    // sf->update( Qref );
    
    // Iterative solution a la Todini paper
    // for(int it=0; it<max_it; ++it){
    //   double Qref = (q_sf + q_sf_in)/2.0;
    //   sf->update( Qref );
    //   if( id == 0 ){
    // 	Rcpp::Rcout << it << " " << Qref << " " << sf->kappa << " " << sf->eta << std::endl;
    //   }
    //   q_sf = std::max(0.0, (z - (sf->kappa*sf->eta)*q_sf_in) / (Dt + sf->kappa*(1.0 - sf->eta)) );
    //   //q_sf = std::max(0.0, (2*sf->Cs*s0 - (1-sf->Ds)*qin) / (1+sf->Ds+2*sf->Cs*Dt) );
    // }
    //s_sf = z - Dt*q_sf;
    //z = std::max(0.0, sf->kappa*(sf->eta*q_sf_in + (1.0-sf->eta)*q_sf) );
    // if( std::abs(s_sf - z) > 1e-6 ){
    //   Rcpp::Rcout << "At end of s_sf update" << id << std::endl;
    //   Rcpp::Rcout << "     s_sf:  " << s_sf << std::endl;
    //   Rcpp::Rcout << "     s_sf alt:  " << z << std::endl;
    //   Rcpp::Rcout << "     kappa:  " << sf->kappa << std::endl;
    //   Rcpp::Rcout << "     eta:  " << sf->eta << std::endl;
    //   Rcpp::Rcout << "     q_sf_in:  " << q_sf_in << std::endl;
    //   Rcpp::Rcout << "     q_sf:  " << q_sf << std::endl;
    //   Rcpp::Rcout << "     v_sf_rz:  " << v_sf_rz << std::endl;
    // }
  }
     
  // redistributed the flows
  //lateral_redistribution(); //vec_q_sf_in,vec_q_sz_in);
    
  // single HRU mass balance for development
  mass_ballance[0] += Dt*(q_sf_in - q_sf) - v_sf_rz - s_sf; // surface
  mass_ballance[1] += Dt*(precip - aet) + v_sf_rz - v_rz_uz - s_rz; // root zone
  mass_ballance[2] += v_rz_uz - v_uz_sz - s_uz; // unsaturated zone
  mass_ballance[3] += Dt*(q_sz - q_sz_in) - v_uz_sz - s_sz; // saturated zone
  z = 0;
  for(int ii=0; ii<4; ii++){
    z = std::max( z, std::abs(mass_ballance[ii]));
  }
  // if( z > 1e-10){
  //     Rcpp::Rcout << "At end of " << id << std::endl; //": " << mass_ballance << " : " << std::endl;
  //     Rcpp::Rcout << "     s_sf:  " << mass_ballance[0] << std::endl;
  //     Rcpp::Rcout << "     s_rz:  " << mass_ballance[1] << std::endl;
  //     Rcpp::Rcout << "     s_uz:  " << mass_ballance[2] << std::endl;
  //     Rcpp::Rcout << "     s_sz:  " << mass_ballance[3] << std::endl;
  // }
}

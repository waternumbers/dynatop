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
  s_sf(states_[0]), s_rz(states_[1]), s_uz(states_[2]), s_sz(states_[3]),
  q_sf(states_[4]), q_sz(states_[5]), // these states are only used in the flux limiters
  area(properties_[0]),width(properties_[1]),Dx(properties_[2]), // Dx passed in this way for future flexibility
  sf_param(sf_param_), s_rzmax(rz_param_[0]),
  t_d(uz_param_[0]), sz_param(sz_param_),
  precip_lnk_id(precip_lnk_id_), precip_lnk_frc(precip_lnk_frc_),
  pet_lnk_id(pet_lnk_id_), pet_lnk_frc(pet_lnk_frc_),
  sf_lnk_id(sf_lnk_id_), sf_lnk_frc(sf_lnk_frc_),
  sz_lnk_id(sz_lnk_id_), sz_lnk_frc(sz_lnk_frc_)
{

  // initialise the surface flux object
  switch(sf_type_){
  case 1:
    // constant celerity kinematic wave
    sf = std::make_unique<sfc_cnstCD>( sf_param_, Dx );
    break;
  }

  // initialise the saturated flux object
  switch(sz_type_){
  case 1:
    //exp
    sz = std::make_unique<szc_exp>( sz_param_, properties_[3], area, width, Dx); // properites_[3] is sbar
    break;
  }

}

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

  //Rcpp::Rcout << "id: " << id << std::endl;
  
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

  //Rcpp::Rcout << vec_q_sf_in[id] << " " << vec_q_sz_in[id] << " " << sz->q_szmax << std::endl;
  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id] - q_sz_in;

  //Rcpp::Rcout << "Process inflows" << std::endl;
  //Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << a_sf_in << " " << a_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;

  // only water at surface if inflow can't be absorbed so max downward flux is
  r_sf_rz = q_sf_in;
  
  // if steady state then passed downward flux straight to uz
  r_rz_uz = r_sf_rz;

  double r_inj = area * std::min( r_uz_sz_0, 1/t_d ); // amount of water injected, can't be more then r_uzmax else can't be caused by storage in uz
  Rcpp::Rcout << r_inj << " " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  r_uz_sz = std::min( r_inj + r_rz_uz, area/t_d ); // ensure downward flux is possible
  Rcpp::Rcout << r_inj << " " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  r_uz_sz = std::min( r_uz_sz , sz->q_szmax - q_sz_in ); // revise down so that outflow from saturated zone is possible
  Rcpp::Rcout << r_inj << " " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  r_rz_uz = r_uz_sz - r_inj; // revise down if enough to reach maximum inflow to saturated zone
  Rcpp::Rcout << r_inj << " " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;

  //Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << a_sf_in << " " << a_sz_in << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;

  
  // initialise saturated zone by bisection to match q_sz
  q_sz = r_uz_sz + q_sz_in;

  double z,e;
  std::pair<double,double> bnd(0.0, sz->sz_max), fbnd(0.0,0.0);
  fbnd.first = sz->fq(bnd.first,q_sz) - q_sz;
  fbnd.second = sz->fq(bnd.second,q_sz) - q_sz;
  Rcpp::Rcout << "HRU " << id << " Initialise saturated zone " << q_sz << " " << r_inj << std::endl;
  Rcpp::Rcout << "Lower limit " << bnd.first << " " << fbnd.first << std::endl;
  Rcpp::Rcout << "Upper limit " << bnd.second << " " << fbnd.second << std::endl;
  // bisection
  z = (bnd.second + bnd.first)/ 2.0 ;
  e = sz->fq(z,q_sz) - q_sz;
  //Rcpp::Rcout << "first guess " << z << " " << e << std::endl;
  double it(0.0);
  while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
    if(e<=0.0){ bnd.first = z; fbnd.first = e;}
    if(e>=0.0){ bnd.second = z; fbnd.second = e;}
    z = (bnd.second + bnd.first)/ 2.0 ;
    e = sz->fq(z,q_sz) - q_sz;
    //Rcpp::Rcout << z << " " << e << std::endl;
    it +=1;
  }
  if(it > max_it){
    Rcpp::warning("SZ: No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		  it, bnd.second - bnd.first,e);
  }
  s_sz = z;

  s_uz = t_d*r_uz_sz*s_sz/area; // compute unsaturated zone storage
  if( r_uz_sz < 0.0 ){ r_rz_uz =  r_uz_sz; } // should only occur when saturated limited
  
  // solve rz
  if( r_rz_uz == 0.0 ){
    s_rz = s_rzmax * s_rz_0 * area;
  }else{
    s_rz = s_rzmax * area;
  }
  // balance flux through root zone
  r_sf_rz = r_rz_uz;

  // solve surface
  q_sf = q_sf_in - r_sf_rz;
  bnd.first = 0.0;
  bnd.second = 0.1*area;
  fbnd.first = sf->fq(bnd.first,q_sf_in) - q_sf;
  fbnd.second = sf->fq(bnd.second,q_sf_in) - q_sf;

  if( fbnd.first==0.0 ){
    s_sf = 0.0;
  }else{
    // expand
    int it(0);
    while( (fbnd.second) < 0.0 and (it <= max_it) ) {
      bnd.first=bnd.second;
      fbnd.first = fbnd.second;
      bnd.second += bnd.second;
      fbnd.second = sf->fq(bnd.second,q_sf_in) - q_sf;
      //Rcpp::Rcout << "estimated bounds: " << bnd.first << " " << bnd.second << std::endl;
      //Rcpp::Rcout << "estimated bound values: " << fbnd.first << " " << fbnd.second << std::endl;
      it+=1;
    }
    //Rcpp::Rcout << "Initial bounds after wide: " << bnd.first << " " << bnd.second << std::endl;
    //Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;
 
    // bisection
    z = (bnd.first + bnd.second) / 2.0; // (bnd.first + bnd.second) / 2.0;
    e = sf->fq(bnd.second,q_sf_in) - q_sf;
    //Rcpp::Rcout << "first guess " << z << " " << e << std::endl;
    it = 0;
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      z = (bnd.second + bnd.first)/ 2.0 ;
      e = sf->fq(bnd.second,q_sf_in) - q_sf;
      //Rcpp::Rcout << z << " " << e << std::endl;
      it +=1;
    }
    s_sf = z;
  }

  // Rcpp::Rcout << sf->fa(q_sf_in) << " " << sz->fa(q_sz_in) << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << precip << " " << pet << " " << aet << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  // Rcpp::Rcout << sf->fa(q_sf) << " " << sz->fa(q_sz) << std::endl;

  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);  
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
  
  //Rcpp::Rcout << "Process inflows" << std::endl;
  //Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << a_sf_in << " " << a_sz_in << std::endl;
  
  // single HRU mass balance for development
  double mass_ballance = s_sf + s_rz + s_uz - s_sz + Dt*( precip + q_sz_in + q_sf_in );
  
  // compute first downward flux estimate from surface zone and root zone
  // limites by 
       r_sf_rz = (s_sf/Dt) + q_sf_in;
  
  //Rcpp::Rcout << "Max from sf" << std::endl;
  //Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  //Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // change r_rz_uz to present the maximum downward flux
  r_rz_uz = std::max(0.0 ,
		     (s_rz - (area*s_rzmax))/Dt  + precip - pet + r_sf_rz);

  //Rcpp::Rcout << "Max from rz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;

  // bounds for search of s_sz
  double q_lower = std::min( q_sz, q_sz_in); //
  double q_upper = std::min( sz->q_szmax, std::max( q_sz, q_sz_in) + area/t_d ); // a bit ropey using area/t_d
  std::pair<double,double> bnd(0.0, sz->sz_max), fbnd(0.0,0.0);
  fbnd.first = fsz(bnd.first,q_sz_in,q_lower,q_upper,Dt);
  fbnd.second = fsz(bnd.first,q_sz_in,q_lower,q_upper,Dt);

  //Rcpp::Rcout << "Initial bounds: " << bnd.first << " " << bnd.second << std::endl;
  //Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;
  //Rcpp::Rcout << "Limits: " << vtol<< " " << etol << " " << max_it << std::endl;
  

  
  double z,e;
  if( fbnd.first >= 0.0 ){
    // then reached point where there is no defict
    z = bnd.first;
  }else{
    
    // // expand
    // int it(0);
    // while( (fbnd.second) < 0.0 and (it <= max_it) ) {
    //   bnd.first=bnd.second;
    //   fbnd.first = fbnd.second;
    //   bnd.second += bnd.second;
    //   fbnd.second = fsz(bnd.second,Dt,q_lower,q_upper);
    //   //Rcpp::Rcout << "estimated bounds: " << bnd.first << " " << bnd.second << std::endl;
    //   //Rcpp::Rcout << "estimated bound values: " << fbnd.first << " " << fbnd.second << std::endl;
    //   it+=1;
    // }
    // //Rcpp::Rcout << "Initial bounds after wide: " << bnd.first << " " << bnd.second << std::endl;
    // //Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;

    // bisection
    z = s_sz; // (bnd.first + bnd.second) / 2.0;
    e = fsz(z,q_sz_in,q_lower,q_upper,Dt);
    //Rcpp::Rcout << "first guess " << z << " " << e << std::endl;
    int it(0);
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      z = (bnd.second + bnd.first)/ 2.0 ;
      e = fsz(z,q_sz_in,q_lower,q_upper,Dt);
      //Rcpp::Rcout << z << " " << e << std::endl;
      it +=1;
    }
    if(it > max_it){
      Rcpp::warning("SZ: No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		    it, bnd.second - bnd.first,e);
    }
    z = bnd.second; // this should ensure that q_sz is positive in the following mass balance
  }
  // update values cased on uz formulationa and mass balance ;
  //Rcpp::Rcout << "after bisection..." << std::endl;
  //Rcpp::Rcout << bnd.first << " " << bnd.second << " " << z << std::endl;
  //Rcpp::Rcout << fbnd.first << " " << fbnd.second << " " << fsz(z,Dt) << std::endl;
  
  
  r_uz_sz = std::min(  area / t_d, (s_uz + Dt*r_rz_uz)/( (z*t_d) + area*Dt) );
  q_sz = ((z - s_sz)/Dt) + q_sz_in + r_uz_sz;
  s_sz = z;
  
  if( (s_sz < 0) ){
    Rcpp::Rcout << "HRU " << id << ": s_sz value " << s_sz << " out of bounds" << std::endl;
  }
  if( (q_sz < 0) | (q_sz > sz->q_szmax) ){
    Rcpp::Rcout << "HRU " << id << ": q_sz value " << q_sz << " is out of bounds" << std::endl;
  }
  
  //Rcpp::Rcout << "Solved sz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // solve unsaturated zone
  z = std::min( s_sz, s_uz + Dt*(r_rz_uz - r_uz_sz) );
  r_rz_uz = ( (z - s_uz)/Dt ) + r_uz_sz;
  s_uz = z;
  
  //Rcpp::Rcout << "Solved uz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // solve root zone
  r_sf_rz = std::min( r_sf_rz , ( ((area*s_rzmax)-s_rz) / Dt ) - precip + r_rz_uz + pet );
  //z = ( s_rzmax / (s_rzmax + pet*Dt) ) * ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) );
  s_rz = ( (area*s_rzmax) / ( (area*s_rzmax) + pet*Dt) ) * ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) );
  aet = pet * s_rz / (s_rzmax*area);
  
  //Rcpp::Rcout << "Solved rz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  //Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  //Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  
  // solve for surface by bisection
  q_lower = std::max( 0.0, std::min( q_sf, q_sf_in) - std::max(r_sf_rz,0.0) );
  q_upper = std::max( q_sf, q_sf_in) - std::min(r_sf_rz,0.0);
    
  bnd.first = 0.0;
  bnd.second = s_sf+1000;
  fbnd.first = fsf(bnd.first,q_sf_in,q_lower,q_upper,Dt);
  fbnd.second = fsf(bnd.second,q_sf_in,q_lower,q_upper,Dt);

  if( fbnd.first==0.0 ){
    s_sf = 0.0;
  }else{
    // expand
    int it(0);
    while( (fbnd.second) < 0.0 and (it <= max_it) ) {
      bnd.first=bnd.second;
      fbnd.first = fbnd.second;
      bnd.second += bnd.second;
      fbnd.second = fsf(bnd.second,q_sf_in,q_lower,q_upper,Dt);
      //Rcpp::Rcout << "estimated bounds: " << bnd.first << " " << bnd.second << std::endl;
      //Rcpp::Rcout << "estimated bound values: " << fbnd.first << " " << fbnd.second << std::endl;
      it+=1;
    }
    //Rcpp::Rcout << "Initial bounds after wide: " << bnd.first << " " << bnd.second << std::endl;
    //Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;
 
    // bisection
    z = s_sf; // (bnd.first + bnd.second) / 2.0;
    e = fsf(z,q_sf_in,q_lower,q_upper,Dt);
    //Rcpp::Rcout << "first guess " << z << " " << e << std::endl;
    it = 0;
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      z = (bnd.second + bnd.first)/ 2.0 ;
      e = fsf(z,q_sf_in,q_lower,q_upper,Dt);
      //Rcpp::Rcout << z << " " << e << std::endl;
      it +=1;
    }
    if(it > max_it){
      Rcpp::warning("SF: No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		    it, bnd.second - bnd.first,e);
    }
    z = bnd.first; // this should ensure that q_sz is positive in the following mass balance
  }
  q_sf = ((z-s_sf)/Dt) + q_sf_in - r_sf_rz;
  
  //Rcpp::Rcout << "Solved all" << std::endl;
  // Rcpp::Rcout << "inflow " << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << "external input: " << precip << " " << pet << " " << aet << std::endl;
  //Rcpp::Rcout << "vertical flux " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  //Rcpp::Rcout << "states " << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  //Rcpp::Rcout << "outflow " << q_sf << " " << q_sz << std::endl;
 
  if( (s_sf < -1e-10) ){
    Rcpp::Rcout << "HRU " << id << ": s_sf value " << s_sf << " is less then 0" << std::endl;
  }
  if( (q_sf < 0) ){
    Rcpp::Rcout << "HRU " << id << ": q_sf value " << q_sf << " is less then 0" << std::endl;
  }
 
  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
  
  //Rcpp::Rcout << "Outflows after partitioning" << std::endl;
  //Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  
  // Rcpp::Rcout << s_sf << " " << q_sf << " " << q_sf_in << " " << r_sf_rz << " " << Dt << std::endl;
  // // single HRU mass balance for development
  mass_ballance = mass_ballance - Dt*(q_sz + q_sf + aet);
  mass_ballance = mass_ballance - (s_sf + s_rz + s_uz - s_sz);
  if( std::abs(mass_ballance) > 1e-9){
    Rcpp::Rcout << "At end: " << mass_ballance << " : " << std::endl;
    Rcpp::Rcout << s_sf  << " : " << s_rz << " : " << s_uz << " : " << s_sz << " : " << q_sz <<  " : " << q_sf << std::endl;
  }
}


double hru::fsz(double& s, double &qin, double &q_lower, double &q_upper, double const &Dt){ // compute for saturated zone function
  double r = std::min(  area / t_d, (s_uz + Dt*r_rz_uz)/( (s*t_d) + area*Dt) ); // compute recharge
  double q = sz->fq(s,qin);
  q = std::min( std::max(q ,q_lower),q_upper );
  //Rcpp::Rcout << q << " " << s << " " << r << std:endl;
  return( s - s_sz + Dt * (qin + r - q) );
}

double hru::fsf(double& s, double &qin, double &q_lower, double &q_upper, double const &Dt){ // compute for saturated zone function
  double q = sf->fq(s,qin);
  q = std::min( std::max(q ,q_lower),q_upper );
  //Rcpp::Rcout << q << " " << s << " " << r << std:endl;
  return( s - s_sz + Dt * (q + r_sf_rz - qin) );
}

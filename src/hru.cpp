#include "hru.h"


hru::hru(int const id_,
	 double s_sf_,double s_rz_,double s_uz_,double s_sz_,
	 double const area_, double const s_bar_, double const width_, 
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
  s_sf(s_sf_), s_rz(s_rz_), s_uz(s_uz_), s_sz(s_sz_),
  area(area_),
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
    sf = std::make_unique<sfc_cnstC>( sf_param_, area, width_);
    break;
  case 2:
    // constant celerity and diffusivity diffuse wave
    sf = std::make_unique<sfc_cnstCD>( sf_param_, area, width_);
    break;
  // case 2:
  //   // generic storage discharge
  //   sf = std::make_unique<sfc_sd>( sf_param_, area, width_);
  //   break;
  //   //  case 3:
  //   // diffusive routing with rectangular channel
  //   // sf_func = std::make_unique<qsf_dfr>( sf_param_, channel_length_);
  //   //break;
  // case 4:
  //   sf = std::make_unique<sfc_msk>( sf_param_, area, width_ );
  }

  // initialise root zone
  //t_d = uz_param_[0];
  //s_rzmax = rz_param[0];

  // initialise the saturated flux object
  switch(sz_type_){
  case 1:
    //bexp
    sz = std::make_unique<szc_bexp>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*sz_param_[2];
    break;
  case 2:
    //dexp
    //sz = std::make_unique<szc_dexp>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*1e34;
    break;
  case 3:
    //cnst
    sz = std::make_unique<szc_cnst>( sz_param_, area, width_);
    //s_szmax = area*sz_param_[1];
    break;
  }

};

void hru::lateral_redistribution(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in){
  for(uint ii=0; ii<sf_lnk_id.size(); ++ii){
    const int &i = sf_lnk_id[ii];
    const double &f = sf_lnk_frc[ii];
    vec_q_sf_in[i] += area * f * q_sf;
  }
  for(uint ii=0; ii<sz_lnk_id.size(); ++ii){
    const int &i = sz_lnk_id[ii];
    const double &f = sz_lnk_frc[ii];
    vec_q_sz_in[i] += area * f * q_sz;
  }
};

void hru::update_met(std::vector<double> &obs){
  precip = 0.0;
  for(uint ii=0; ii<precip_lnk_id.size(); ++ii){
    const int &i = precip_lnk_id[ii];
    const double &f = precip_lnk_frc[ii];
    precip += f * obs[i];
  }
  pet = 0.0;
  for(uint ii=0; ii<pet_lnk_id.size(); ++ii){
    const int &i = pet_lnk_id[ii];
    const double &f = pet_lnk_frc[ii];
    pet += f * obs[i];
  }
};

void hru::init(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double s_rz_0, double r_uz_sz_0,
	       double const &vtol, double const &etol, int const &max_it){

  //Rcpp::Rcout << "id: " << id << std::endl;
  
  if(area == 0.0){ // if the HRU has no area just pass on flow and set states to 0
    q_sf_in = vec_q_sf_in[id] / area;
    q_sf = q_sf_in;
    q_sz_in = vec_q_sz_in[id] / area;
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
  q_sz_in = std::min( vec_q_sz_in[id] / area, sz->q_szmax) ;
  q_sf_in = (vec_q_sf_in[id] / area) + (vec_q_sz_in[id] / area) - q_sz_in;

  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << kappa_sf << " " << eta_sf << std::endl;
  // Rcpp::Rcout << kappa_sz << " " << eta_sz << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;

  // only water at surface if inflow can't be absorbed so max downward flux is
  r_sf_rz = q_sf_in/(1-eta_sf);
  
  // if steady state then passed downward flux straight to uz
  r_rz_uz = r_sf_rz;

  // work out max flux from us to sz
  r_uz_sz = std::min( 1/t_d, r_rz_uz + r_uz_sz_0);

  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << kappa_sf << " " << eta_sf << std::endl;
  // Rcpp::Rcout << kappa_sz << " " << eta_sz << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  
  // initialise saturated zone
  // at steady state the grdient of d(D-s_sz)/dt is zero
  std::pair<double,double> bnd(0.0,sz->D), fbnd(0.0,0.0);
  fzi(bnd.first, fbnd.first); // changes kappa_sz and eta_sz
  fzi(bnd.second, fbnd.second); // changes kappa_sz and eta_sz
  
  // Rcpp::Rcout << "Initial bounds: " << bnd.first << " " << bnd.second << std::endl;
  // Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;			
  if( fbnd.first >= 0.0 ){
    // then saturated
    s_sz = 0.0;
    sz->fke( kappa_sz, eta_sz, s_sz);
  }else{
    if(fbnd.second <= 0.0){
      // then fully drained
      s_sz = bnd.second;
      sz->fke( kappa_sz, eta_sz, s_sz);
    }else{
      
      // bisection
      s_sz = (bnd.first + bnd.second) /2.0;
      double e;
      fzi(s_sz, e);
      int it(0.0);
      while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
	if(e<=0.0){ bnd.first = s_sz; fbnd.first = e;}
	if(e>=0.0){ bnd.second = s_sz; fbnd.second = e;}
	//e = fbnd.first / (fbnd.second - fbnd.first);
	//z = (bnd.first * (1 + e)) - (e * bnd.second);
	s_sz = (bnd.second + bnd.first)/ 2.0 ;
	fzi(s_sz, e);
	it +=1;
      }
      if(it > max_it){
	Rcpp::warning("No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		      max_it, bnd.second - bnd.first,e);
      }
    }
    // update values cased on uz formulationa and mass balance ;
    // Rcpp::Rcout << "after bisection..." << std::endl;
    // Rcpp::Rcout << bnd.first << " " << bnd.second << " " << z << std::endl;
    // Rcpp::Rcout << fbnd.first << " " << fbnd.second << " " << std::endl;
  }

  // this point s_sz, kappa_sz and eta_sz are all updated
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << kappa_sf << " " << eta_sf << std::endl;
  // Rcpp::Rcout << kappa_sz << " " << eta_sz << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  
  // compute r_uz_sz
  r_uz_sz =  (1.0/(kappa_sz*(1.0-eta_sz)))*(sz->D - s_sz) - q_sz_in/(1.0-eta_sz);
  // compute s_uz
  s_uz = r_uz_sz * t_d * s_sz;
  // by mass balance
  r_rz_uz = std::min(r_rz_uz, r_uz_sz);
  r_sf_rz = r_rz_uz;
  // for root zone
  if( r_rz_uz == 0.0 ){
    s_rz = s_rzmax*s_rz_0;
  }else{
    s_rz = s_rzmax;
  }
  // solve surface - assumes constant kappa and eta
  sf->fke(kappa_sf,eta_sf,s_sf);
  s_sf = kappa_sf*( q_sf_in - (1-eta_sf)*r_sf_rz );

  // solve for lateral fluxes using mass balance
  q_sf = q_sf_in - r_sf_rz;
  q_sz = q_sz_in + r_uz_sz;
  
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << kappa_sf << " " << eta_sf << std::endl;
  // Rcpp::Rcout << kappa_sz << " " << eta_sz << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  
  // Repartition the outflows
  q_sf = q_sf + q_sz;
  q_sz = std::min( q_sz, sz->q_szmax) ;
  q_sf -= q_sz;

  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);  
}


void hru::step(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double const &vtol, double const &etol, int const &max_it, double const &Dt)
{
  if(area == 0.0){ // if the HRU has no area just pass on flow
    q_sf_in = vec_q_sf_in[id] / area;
    q_sf = q_sf_in;
    q_sz_in = vec_q_sz_in[id] / area;
    q_sz = q_sz_in;
    lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    return;
  }

  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] / area, sz->q_szmax) ;
  q_sf_in = (vec_q_sf_in[id] / area) + (vec_q_sz_in[id] / area) - q_sz_in;

  // update a_sf_in and a_sz_in
  a_sf_in = sf->fa(q_sf_in);
  a_sz_in = sz->fa(q_sf_in);
  
  // single HRU mass balance for development
  double mass_ballance = s_sf + s_rz + s_uz - s_sz + Dt*( precip + q_sz_in + q_sf_in );
  
  // compute first downward flux estimate from surface zone and root zone
  r_sf_rz =  s_sf/Dt + q_sf_in;
  
  // Rcpp::Rcout << "Max from sf" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // change r_rz_uz to present the maximum downward flux
  r_rz_uz = std::max(0.0 ,
		     (s_rz - s_rzmax)/Dt  + precip - pet + r_sf_rz);

  // Rcpp::Rcout << "Max from rz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;

  // bounds for search of s_sz
  std::pair<double,double> bnd(0.0, sz->q_szmax), fbnd(0.0,0.0);
  fbnd.first = fz(bnd.first,Dt);
  fbnd.second = fz(bnd.second,Dt);

  // Rcpp::Rcout << "Initial bounds: " << bnd.first << " " << bnd.second << std::endl;
  // Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;

  doubel z = q_sz;
  if( fbnd.second >= 0.0 ){
    // then reached point of maximum outflow
    z = fbnd.second;
  }else{
    // bisection
    double e = fz(z,Dt);
    int it(0.0);
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      z = (bnd.second + bnd.first)/ 2.0 ;
      e = fz(z,Dt);
      it +=1;
    }
    if(it > max_it){
      Rcpp::warning("No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		    max_it, bnd.second - bnd.first,e);
    }
  }
  // update values cased on uz formulationa and mass balance ;
  // Rcpp::Rcout << "after bisection..." << std::endl;
  // Rcpp::Rcout << bnd.first << " " << bnd.second << " " << z << std::endl;
  // Rcpp::Rcout << fbnd.first << " " << fbnd.second << " " << std::endl;
  
  
  //r_uz_sz = std::min( 1/t_d , (s_uz + Dt*r_rz_uz) / ( (z*t_d) + Dt ) );
  q_sz = z;
  z = sz->fsz(q_sz_in,q_sz);
  r_uz_sz = ((s_sz - z)/Dt) - q_sz_in + q_sz;
  
  
  if( (s_sz < 0) ){
    Rcpp::Rcout << "HRU " << id << ": s_sz value " << s_sz << " out of bounds" << std::endl;
  }
  if( (q_sz < -1e-10) | (q_sz > sz->q_szmax) ){
    Rcpp::Rcout << "HRU " << id << ": q_sz value " << q_sz << " is out of bounds" << std::endl;
  }
  
  // Rcpp::Rcout << "Solved sz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // solve unsaturated zone
  z = std::min( s_sz, s_uz + Dt*(r_rz_uz - r_uz_sz) );
  r_rz_uz = ( (z - s_uz)/Dt ) + r_uz_sz;
  s_uz = z;
  
  // Rcpp::Rcout << "Solved uz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // solve root zone
  r_sf_rz = std::min( r_sf_rz , ( (s_rzmax-s_rz) / Dt ) - precip + r_rz_uz + pet );
  //z = ( s_rzmax / (s_rzmax + pet*Dt) ) * ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) );
  s_rz = ( s_rzmax / (s_rzmax + pet*Dt) ) * ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) );
  aet = pet * s_rz / s_rzmax;
  
  // Rcpp::Rcout << "Solved rz" << std::endl;
  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // solve for surface
  bnd.first=0.0;
  bnd.second=std::min(s_sf, 0.01);
  fbnd.first = fsf(bnd.first,Dt);
  fbnd.second = fsf(bnd.second,Dt);

  double z = bnd.second;
  if( fbnd.first == 0.0 ){
    z = bnd.first;
  }else{
    // expand
    while( fbnd.second < 0.0) {
      bnd.first=bnd.second;
      fbnd.first = fbnd.second;
      bnd.second += bnd.second;
      fbnd.second = fz(bnd.second,Dt);
    }
    // bisect
    z = (bnd.second + bnd.first)/ 2.0 ;
    double e = fsf(z,Dt);
    int it(0.0);
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      z = (bnd.second + bnd.first)/ 2.0 ;
      e = fsf(z,Dt);
      it +=1;
    }
    if(it > max_it){
      Rcpp::warning("No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		    max_it, bnd.second - bnd.first,e);
    }
  }

  q_sf = (s_sf - z)/ Dt + q_sf_in - r_sf_rz;
  s_sf = z;
  	
  if( (s_sf < 0) ){
    Rcpp::Rcout << "HRU " << id << ": s_sf value " << s_sf << " is less then 0" << std::endl;
  }
  if( (q_sf < -1e-10) ){
    Rcpp::Rcout << "HRU " << id << ": q_sf value " << q_sf << " is less then 0" << std::endl;
  }
 
  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
  
  // Rcpp::Rcout << "Outflows after partitioning" << std::endl;
  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  
  // Rcpp::Rcout << s_sf << " " << q_sf << " " << q_sf_in << " " << r_sf_rz << " " << Dt << std::endl;
  // // single HRU mass balance for development
  mass_ballance = mass_ballance - Dt*(q_sz + q_sf + aet);
  mass_ballance = mass_ballance - (s_sf + s_rz + s_uz - s_sz);
  if( std::abs(mass_ballance) > 1e-10){
    Rcpp::Rcout << "At end: " << mass_ballance << " : " << std::endl;
    Rcpp::Rcout << s_sf  << " : " << s_rz << " : " << s_uz << " : " << s_sz << " : " << q_sz <<  " : " << q_sf << std::endl;
  }
}



void hru::fzi(double& x, double &e){ // compute for saturated zone function
  // update using fke
  sz->fke(kappa_sz,eta_sz,x);
  e = q_sf_in/(1-eta_sz) + r_uz_sz - (sz->D-x)/(kappa_sz*(1.0-eta_sz));
}
 


void hru::fsz(double& q, double const &Dt){ // compute for saturated zone function
  double s = sz -> fs(q,q_sz_in); // compute storage
  double r = std::min(  1 / t_d, (s_uz + Dt*r_rz_uz)/( (s*t_d) + Dt) ); // compute recharge
  return( s - s_sz + Dt * (q_sz_in + r - q) );
}

void hru::fsf(double& s, double const &Dt){ // compute for saturated zone function
  double q = sf -> fq(s,q_sz_in); // compute storage
  return( s - s_sf - Dt * (q_sf_in - r_sf_rz - q) );
}


  q_sz = ( (sz->D - s_sz) + (Dt+kappa_sz*eta_sz)*q_sz_in + Dt*r_uz_sz ) /
    (kappa_sz*(1-eta_sz)+Dt);
  if(q_sz > sz->q_szmax){
    // take off r_uz_sz component to q_sz
    q_sz -= (Dt*r_uz_sz) /  (kappa_sz*(1-eta_sz)+Dt);
    // work out new r_uz_sz
    r_uz_sz = (kappa_sz*(1-eta_sz)+Dt) * ( sz->q_szmax - q_sz ) / Dt;
    // set q_sz to maximum value
    r_uz_sz = 
  r_uz_sz = std::min(  r_uz_sz, rmax);
  
  e = x - s_sz - Dt*( (1.0/(kappa_sz*(1.0-eta_sz)))*(sz->D - x) - q_sz_in/(1.0-eta_sz) - r_uz_sz);
}

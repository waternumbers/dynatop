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
    // constant velocity with raf tank:
    sf = std::make_unique<sfc_cnst>( sf_param_, area, width_);
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
    sz = std::make_unique<szc_cnst>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*sz_param_[1];
    break;
  }

  // compute properties
  Dx = area/ width_;
  q_szmax = sz -> fq_szmax();
  D_sz = sz -> fD();
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

  // redivide the inflow so q_sz_in is less then q_szmax
  q_sz_in = std::min( vec_q_sz_in[id] / area, q_szmax) ;
  q_sf_in = (vec_q_sz_in[id] / area) + (vec_q_sz_in[id] / area) - q_sz_in;

  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  
  // only water at surface so inflow can't be absorbed so max downward flux is
  r_sf_rz = q_sf_in;
  
  // if steady state then passed downward flux straight to uz
  r_rz_uz = r_sf_rz;

  // work out max flux from us to sz
  r_uz_sz = std::min( 1/t_d, r_rz_uz + r_uz_sz_0);

  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  
  // initialise saturated zone
  //Rcpp::Rcout << "Prior to call: " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  sz -> init(s_sz, q_sz, q_sz_in, r_uz_sz); // revises s_sz, q_sz and possibly r_uz_sz (if saturated)
  //Rcpp::Rcout << "After call: " << r_uz_sz << std::endl;

  //Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  
  // initialise unsaturated zone
  s_uz = r_uz_sz * t_d * s_sz;

  //Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // correct flux from rz to uz
  r_rz_uz = r_uz_sz - ( r_rz_uz + r_uz_sz_0);
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;

  // solve root zone
  if( (r_rz_uz == 0) & ( r_sf_rz==0) ){
    s_rz = s_rz_0 * s_rzmax;
  }else{
    s_rz = s_rzmax;
    r_sf_rz = r_rz_uz;
  }
  //Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  //Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;

  // initialise the surface zone
  sf -> init(s_sf, q_sf, q_sf_in, r_sf_rz); // revises s_sf and q_sf
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // repartition the outflows
  q_sf = q_sf + q_sz;
  q_sz = std::min( q_sz, q_szmax) ;
  q_sf -= q_sz;
  //Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;

  // redistributed the flows
  lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
};


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
  q_sz_in = std::min( vec_q_sz_in[id] / area, q_szmax) ;
  q_sf_in = (vec_q_sz_in[id] / area) + (vec_q_sz_in[id] / area) - q_sz_in;

  // Rcpp::Rcout << q_sf_in << " " << q_sz_in << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;
  
  // single HRU mass balance for development
  double mass_ballance = s_sf + s_rz + s_uz - s_sz + Dt*( precip + q_sz_in + q_sf_in );
  

  // compute first downward flux estimate from surface zone and root zone
  r_sf_rz = sf -> max_down(s_sf, q_sf_in, Dt);

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
  
  
  // bounds for search of outflow
  double z = sz -> fD(); // max storage volume
  std::pair<double,double> bnd(0.0,z), fbnd(0.0,0.0);
  fbnd.first = fz(bnd.first,Dt);
  fbnd.second = fz(bnd.second,Dt);
  
  // Rcpp::Rcout << "Initial bounds: " << bnd.first << " " << bnd.second << std::endl;
  // Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;				
  


  if( fbnd.first >= 0.0 ){
    // then saturated
    z = 0.0;
  }else{
    // bisection
    z = s_sz;
    double e = fz(z,Dt);
    int it(0.0);
    while( (it <= max_it) and ( (bnd.second - bnd.first)>vtol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      //e = fbnd.first / (fbnd.second - fbnd.first);
      //z = (bnd.first * (1 + e)) - (e * bnd.second);
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
  // // Rcpp::Rcout << "after bisection..." << std::endl;
  // Rcpp::Rcout << bnd.first << " " << bnd.second << " " << z << std::endl;
  // Rcpp::Rcout << fbnd.first << " " << fbnd.second << " " << std::endl;
  
  r_uz_sz = std::min( 1/t_d , (s_uz + Dt*r_rz_uz) / ( (z*t_d) + Dt ) );
  q_sz = ( (z - s_sz) / Dt ) + r_uz_sz + q_sz_in;
  s_sz = z;

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
  sf -> solve(s_sf, q_sf, q_sf_in, r_sf_rz, Dt); // changes q_sf and s_sf

  // Rcpp::Rcout << "Solved sf" << std::endl;
  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;
  // Rcpp::Rcout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  // Rcpp::Rcout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << std::endl;

  // repartition the outflows
  q_sf = q_sf + q_sz;
  q_sz = std::min( q_sz, q_szmax) ;
  q_sf -= q_sz;

  // Rcpp::Rcout << q_sf << " " << q_sz << std::endl;

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

double hru::fz(double& x, double const &Dt){ // compute for saturated zone function
  double c = sz -> fc(x);
  double r = std::min(  1 / t_d, (s_uz + Dt*r_rz_uz)/( (x*t_d) + Dt) );
  return x - s_sz - Dt*( (2*c*(D_sz-x)/Dx) - 2*q_sz_in - r);
}

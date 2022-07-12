#include "hru.h"


hru::hru(double &s_sf_,double &s_rz_,double &s_uz_,double &s_sz_,
	 double &q_sf_,double &q_sz_, double &q_sf_in_, double &q_sz_in_,
	 double &precip_, double &pet_, double &e_a_,
	 double &r_sf_sp_, double &r_sf_rz_, double &r_sp_rz_, double &r_rz_uz_, double &r_uz_sz_,
	 double &area_, double &s_bar_, double &width_,
	 int const sf_type_, std::vector<double> const sf_param_,
	 int const sp_type_, std::vector<double> const sp_param_,
	 int const rz_type_, std::vector<double> const rz_param_,
	 int const uz_type_, std::vector<double> const uz_param_,
	 int const sz_type_, std::vector<double> const sz_param_, // saturated flux type and parameters
	 double const &Dt_, double const &ztol_, double const &etol_, int const &max_it_  
	 ):
  s_sf(s_sf_), s_rz(s_rz_), s_uz(s_uz_), s_sz(s_sz_),
  q_sf(q_sf_), q_sz(q_sz_),q_sf_in(q_sf_in_), q_sz_in(q_sz_in_),
  precip(precip_), pet(pet_), e_a(e_a_),
  r_sf_rz(r_sf_rz_), r_rz_uz(r_rz_uz_), r_uz_sz(r_uz_sz_),
  area(area_),
  s_rzmax(s_rzmax_), t_d(t_d_),
  area(area_),
  Dt(Dt_), ztol(ztol_), etol(etol_), max_it(max_it_)
{
  // initialise the surface flux object
  switch(sf_type_){
  case 1:
    // constant velocity with raf tank:
    sf_func = std::make_unique<qsf_cnst>( sf_param_, area, width_);
    break;
  case 2:
    // generic storage discharge
    sf_func = std::make_unique<qsf_sd>( sf_param_, area, width_);
    break;
    //  case 3:
    // diffusive routing with rectangular channel
    // sf_func = std::make_unique<qsf_dfr>( sf_param_, channel_length_);
    //break;
  case 4:
    sf_func = std::make_unique<qsf_msk>( sf_param_, area, width_ );
  }

  // initialise the spill flux object
  switch(sp_type_){
  case 1:
    // constant velocity with raf tank:
    sp_func = std::make_unique<qsf_cnst>( sp_param_, area, width_);
    break;
  case 2:
    // generic storage discharge
    sp_func = std::make_unique<qsf_sd>( sp_param_, area, width_);
    break;
    //  case 3:
    // diffusive routing with rectangular channel
    // sf_func = std::make_unique<qsf_dfr>( sf_param_, channel_length_);
    //break;
  case 4:
    sp_func = std::make_unique<qsf_msk>( sp_param_, area, width_ );
  }

  // since multiple options for rz and uz are not implimented
  s_rzmax = rz_param_[0]; 
  r_uzmax = area / uz_param_[0]; // uz_param_[0] is t_d
					
  // initialise the saturated flux object
  switch(sz_type_){
  case 1:
    //exp
    sz_func = std::make_unique<qsz_exp>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*1e34;
    break;
  case 2:
    //bexp
    sz_func = std::make_unique<qsz_bexp>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*sz_param_[2];
    break;
  case 3:
    //dexp
    sz_func = std::make_unique<qsz_dexp>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*1e34;
    break;
  case 4:
    //cnst
    sz_func = std::make_unique<qsz_cnst>( sz_param_, s_bar_, area, width_);
    //s_szmax = area*sz_param_[1];
    break;
  }


  
  // std::vector<double> tmp(3);
  // tmp[0]=0.0; tmp[1] = 0.1; tmp[2] = 1.0;
  // Rcpp::Rcout << "fq(0): " << sz_func -> fq(tmp[0]) << " fq(0.1): " << sz_func -> fq(tmp[1]) << " fq(1.0): " << sz_func -> fq(tmp[2]) <<std::endl;
  // Rcpp::Rcout << "ADt "<< ADt << " As_rz_max: " << As_rz_max << " Ar_sf_max: " << Ar_sf_max << " s_szmax: " << s_szmax << std::endl;
};

void hru::init(double &s_rz_0, double &r_uz_sz_0){

  if(area == 0.0){ // if the HRU has no area just pass on flow and set states to 0
    q_sf = q_sf_in;
    q_sz = q_sz_in;
    s_sf = 0.0;
    s_rz = 0.0;
    s_uz = 0.0;
    s_sz = 0.0;
    return;
  }

  // work out maximum downward flux
  double dummy(-999.9);
  sf_func -> max_down(q_sf_in, r_sf_sp, r_sf_rz); // alter r_sf_rz and r_sf_sp to represent max downward flux
  sp_func -> max_down(r_sf_sp, dummy, r_sp_rz); // alter r_sp_rz to represent max downward flux - summy should eb unaltered

  // if steady state then passed downward flux straight to uz
  // flux from uz to sz is at most min of this + initial downflow and max value
  r_uz_sz = std::min( r_uzmax, r_sp_rz + r_sf_rz + r_uz_sz_0*area );

  // initialise saturated zone
  //Rcpp::Rcout << "Prior to call: " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
  sz_func -> init(s_sz,q_sz,q_sz_in,r_uz_sz); // revises s_sz, q_sz and r_uz_sz (if saturated)
  //Rcpp::Rcout << "After call: " << r_uz_sz << std::endl;

  s_uz = r_uz_sz*s_sz/r_uzmax; // compute unsaturated zone storage
  r_rz_uz = r_uz_sz;
  // solve rz
  if( r_rz_uz == 0.0 ){ // not the best test if r_rz_uz == 0 but downward flux from surface
    s_rz = s_rzmax * s_rz_0;
  }else{
    s_rz = s_rzmax;
  }
  // balance flux through root zone
  r_sf_rz = std::min(r_sf_rz,r_rz_uz);
  r_sp_rz = std::min(r_sp_rz,r_rz_uz - r_sf_rz);

  q_sf = 0.0;
  sf_func -> solve(q_sf, q_sf_in, r_sf_sp, r_sf_rz); // alters q_sf and r_sf_sp
  sp_func -> solve(q_sf, r_sf_sp, dummy, r_sp_rz); // alters q_sf and dummy (should not be altered)

};


void hru::step()
{
  if(area == 0.0){ // if the HRU has no area just pass on flow
    q_sf = q_sf_in;
    q_sz = q_sz_in;
    return;
  }
  
  // single HRU mass balance for development
  double mass_ballance = s_sf + s_sp + s_rz + s_uz - s_sz + Dt*( precip + q_sz_in + q_sf_in );
  

  // compute first downward flux estimate from surface and root zone
  // these are given as \check{r} in documentation
  double dummy(-999.9);
  sf_func -> max_down(q_sf_in, r_sf_sp, r_sf_rz); // alter r_sf_rz and r_sf_sp to represent max downward flux
  sp_func -> max_down(r_sf_sp, dummy, r_sp_rz); // alter r_sp_rz to represent max downward flux - summy should eb unaltered
  
  r_rz_uz = std::max(0.0 ,
		     (s_rz - s_rzmax)/Dt  + precip - pet + r_sf_rz + r_sp_rz);

  //Rcpp::Rcout << "Computed maximum downward fluxes" << std::endl;
  
  double z = sz_func -> D(); // max storage volume
  
  // bounds for search of outflow
  std::pair<double,double> bnd(0.0,z), fbnd(0.0,0.0);
  fbnd.first = fz(bnd.first);
  fbnd.second = fz(bnd.second);
  
  //Rcpp::Rcout << "Initial bounds: " << bnd.first << " " << bnd.second << std::endl;
  //Rcpp::Rcout << "Initial bounds values: " << fbnd.first << " " << fbnd.second << std::endl;							 


  if( fbnd.first >= 0.0 ){
    // then saturated
    z = 0.0;
  }else{
    // bisection
    z = s_sz;
    double e = fz(z);
    while( (it <= max_it) and ( (bnd.second - bnd.first)>ztol ) and (std::abs(e)>etol) ){
      if(e<=0.0){ bnd.first = z; fbnd.first = e;}
      if(e>=0.0){ bnd.second = z; fbnd.second = e;}
      //e = fbnd.first / (fbnd.second - fbnd.first);
      //z = (bnd.first * (1 + e)) - (e * bnd.second);
      z = (bnd.second + bnd.first)/ 2.0 ;
      e = fz(z);
      it +=1;
    }
    if(it > max_it){
      Rcpp::warning("No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		    max_it, bnd.second - bnd.first,e);
    }
  }
  // update values cased on uz formulationa and mass balance ;
  // Rcpp::Rcout << "after bisection..." << std::endl;
  //Rcpp::Rcout << bnd.first << " " << bnd.second << " " << z << std::endl;
  //Rcpp::Rcout << fbnd.first << " " << fbnd.second << " " << e << std::endl;
  r_uz_sz = std::min( r_uzmax, (s_uz + Dt*r_rz_uz) / ( (z/r_uzmax) + Dt ) );
  q_sz = ( (z - s_sz) / Dt ) + r_uz_sz + q_sz_in;
  s_sz = z;
      
  
  // Rcpp::Rcout << "s_szmax = " << s_szmax << std::endl;
  //Rcpp::Rcout << z  << " " << q_sz << " " << s_sz << " " << r_uz_sz << std::endl;

  // solve unsaturated zone
  z = std::min( s_sz, s_uz + Dt*(r_rz_uz - r_uz_sz) );
  r_rz_uz = ( (z - s_uz)/Dt ) + r_uz_sz;
  s_uz = z;
  
  // solve root zone
  //Rcpp::Rcout << precip << " " << pet << " " << s_rzmax << " " << r_sf_rz << " " << s_rz << std::endl;
  
  z = std::min( r_sf_rz + r_sp_rz , ( (s_rzmax-s_rz) / Dt ) - precip + r_rz_uz + pet );
  //z = ( s_rzmax / (s_rzmax + pet*Dt) ) * ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) );
  s_rz = ( s_rzmax / (s_rzmax + pet*Dt) ) * ( s_rz + Dt*(precip + z - r_rz_uz) );
  e_a = pet * s_rz / s_rzmax;
  
  // solve for surface
  if( z < r_sf_rz + r_sp_rz ){
    if( z < 0.0 ){
      // then all inflow to the surface non-spill but no downward flux from s_sp
      r_sf_rz = z;
      r_sp_rz = 0.0;
    }else{
      z = z / (r_sf_rz + r_sp_rz);
      r_sf_rz = z * r_sf_rz;
      r_sp_rz = z * r_sp_rz;
    }
  }
  
  // Rcpp::Rcout << s_sf << " " << q_sf << " " << q_sf_in << " " << r_sf_rz << " " << Dt << std::endl;
  q_sf = 0.0;
  sf_func -> solve(s_sf, q_sf, q_sf_in, r_sf_rz, r_sf_sp, Dt); // computes s_sf and adds to q_sf
  sp_func -> solve(s_sf, q_sf, r_sf_sp, r_sp_rz, z, Dt); // computes s_sp and adds to q_sf z is a dummy and should be unchanged
  
  // Rcpp::Rcout << s_sf << " " << q_sf << " " << q_sf_in << " " << r_sf_rz << " " << Dt << std::endl;
  // // single HRU mass balance for development
  mass_ballance = mass_ballance - Dt*(q_sz + q_sf + e_a);
  mass_ballance = mass_ballance - (s_sf + s_rz + s_uz - s_sz);
  if( std::abs(mass_ballance) > 1e-10){
      Rcpp::Rcout << "At end: " << mass_ballance << " : " << std::endl;
     Rcpp::Rcout << s_sf  << " : " << s_rz << " : " << s_uz << " : " << s_sz << " : " << q_sz <<  " : " << q_sf << std::endl;
  }
}

double hru::fz(double& x){ // compute for saturated zone function
  double c = sz_func -> fc(x);
  double r = std::min( r_uzmax, (s_uz + Dt*r_rz_uz)/( (x/r_uzmax) + Dt) );
  // GOT HERE
  return x - s_sz - Dt*(q - q_sz_in - r);
}

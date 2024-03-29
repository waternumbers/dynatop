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
   //states(states_),
  properties(properties_),
  //q_sf(states_[4]), q_sz(states_[5]), // these states are only used in the flux limiters
  //area(properties_[0]),width(properties_[1]),Dx(properties_[2]), // Dx passed in this way for future flexibility
  sf_param(sf_param_), s_rzmax(rz_param_[0]),
  t_d(uz_param_[0]), sz_param(sz_param_),
  precip_lnk_id(precip_lnk_id_), precip_lnk_frc(precip_lnk_frc_),
  pet_lnk_id(pet_lnk_id_), pet_lnk_frc(pet_lnk_frc_),
  sf_lnk_id(sf_lnk_id_), sf_lnk_frc(sf_lnk_frc_),
  sz_lnk_id(sz_lnk_id_), sz_lnk_frc(sz_lnk_frc_),
  id(id_), 
  s_sf(states_[0]), s_rz(states_[1]), s_uz(states_[2]), s_sz(states_[3])
{
  // change depths to volues for storage limits
  area = properties_[0];
  //s_rzmax = s_rzmax*area;
  s_sf *= area;
  s_rz *= area;
  s_uz *= area;
  s_sz *= area;
  
  
  // initialise the surface flux object
  switch(sf_type_){
  case 1:
    // constant celerity with raf
    sf = std::make_unique<sfc_cnst>( sf_param_, properties );
    break;
  case 2:
    // kinematic with raf
    sf = std::make_unique<sfc_kin>( sf_param_, properties );
    break;
  case 3:
    // compound channel
    sf = std::make_unique<sfc_comp>( sf_param_, properties );
    break;
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
    sz = std::make_unique<szc_cnst>( sz_param_, properties ); // properites_[3] is sbar
    break;
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
    precip += f * area * obs[i];
  }
  pet = 0.0;
  for(long unsigned int ii=0; ii<pet_lnk_id.size(); ++ii){
    const int &i = pet_lnk_id[ii];
    const double &f = pet_lnk_frc[ii];
    pet += f * area * obs[i];
  }
}

void hru::init(std::vector<double> &vec_q_sf_in, std::vector<double> &vec_q_sz_in,
	       double s_rz_0, double r_uz_sz_0,
	       double const &vtol, double const &etol, int const &max_it){

  if(area == 0.0){ // if the HRU has no area just pass on flow and set states to NaN
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
  //q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  //q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id] - q_sz_in;
  q_sz_in = vec_q_sz_in[id];
  q_sf_in = vec_q_sf_in[id];

  // only water at surface if inflow can't be absorbed so max downward flux is
  double r_sf_rz = q_sf_in;
  
  // if steady state then passed downward flux straight to uz
  double r_rz_uz = r_sf_rz;

  // injected water flux into the unsaturated zone
  double r_inj = area * r_uz_sz_0;

  
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
  s_sz = sz->fs(q_sz,q_sz_in);
  q_sz = sz->fq(s_sz,q_sz_in);
  r_uz_sz = q_sz - q_sz_in;
  
  if( std::abs( sz->fq(s_sz,q_sz_in) - q_sz ) > 1e-10 ){
    Rcpp::Rcout << id << " saturated" << std::endl;
    Rcpp::Rcout << q_sz_in << " " << q_sz << std::endl; 
    Rcpp::Rcout << s_sz << " " << sz->fq(s_sz,q_sf_in) << std::endl;
  }
  
  //double z = (q_sz + q_sz_in)/2.0;
  //Rcpp::Rcout << "Initialising saturated zone " << z << std::endl;
  //s_sz = sz->fs(z);

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
  if(area == 0.0){ // if the HRU has no area just pass on flow
    q_sf_in = vec_q_sf_in[id];
    q_sf = q_sf_in;
    q_sz_in = vec_q_sz_in[id];
    q_sz = q_sz_in;
    lateral_redistribution(vec_q_sf_in,vec_q_sz_in);
    return;
  }
    
  // redivide the inflow so q_sz_in is less then q_szmax
  //q_sz_in = std::min( vec_q_sz_in[id] , sz->q_szmax) ;
  //q_sf_in = vec_q_sf_in[id] + vec_q_sz_in[id]  - q_sz_in;
  q_sz_in = vec_q_sz_in[id];
  q_sf_in = vec_q_sf_in[id];

  // single HRU mass balance for development
  std::vector<double> mass_ballance = {s_sf, s_rz, s_uz, s_sz};
  
  // compute first downward flux estimate from surface zone and root zone
  // limites by
  v_sf_rz = s_sf + Dt*q_sf_in;
  
  // change r_rz_uz to present the maximum downward flux
  v_rz_uz = std::max(0.0 ,
		     s_rz - (area*s_rzmax)  + Dt*(precip - pet) + v_sf_rz);

  // search for s_sz
  double lb(0.0), ub(0.0), Hzu(2*vtol), Hzl(2*vtol);
  double z(-999.9);
  
  // test ub=0.0 to see if saturated
  v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*ub + area*Dt), 1/t_d );
  Hzu = ub - s_sz + v_uz_sz + Dt*(q_sz_in - sz->fq(ub,q_sz_in));
  Hzl = Hzu; // since both are at 0 from initialisation
  
  if( Hzu < 0.0 ){
    // need  a numerical solution

    // scale out upper limit until positive
    ub = s_sz + 3*vtol;
    v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*ub + area*Dt), 1/t_d );
    Hzu = ub - s_sz + v_uz_sz + Dt*(q_sz_in - sz->fq(ub,q_sz_in));
    int it(0.0);
    while( (Hzu < 0.0) and (it < max_it) ){
      lb = ub;
      Hzl = Hzu;
      ub += ub;
      v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*ub + area*Dt), 1/t_d );
      Hzu = ub - s_sz + v_uz_sz + Dt*(q_sz_in - sz->fq(ub,q_sz_in));
      it +=1;
    }
    if( Hzu < 0 ){
      Rcpp::warning("SZ: No upper bound found within %i iterations. Difference between bounds is %d.",
		    it, ub - lb); //bnd.second - bnd.first);
    }

    // shrink back to find solution
    it = 0;
    double Hz(2*vtol), iW(0.001);
    while( (Hzu > vtol) and (it < max_it) ){ //((ub -lb) > vtol) and (it < max_it) ){
      iW = Hzu / (Hzu - Hzl);
      iW = std::max(0.001,std::min(iW,0.999));
      z = (iW*lb) + (1.0-iW)*ub;
      v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
      Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - sz->fq(z,q_sz_in));
      
      if( Hz < 0 ){
	lb = z;
	Hzl = Hz;
      }else{
	ub = z;
	Hzu = Hz;
      }
      it += 1; 
    }
    if(it > max_it){
      Rcpp::warning("HRU %i SZ: No solution found within %i iterations. Difference between bounds is %d",
		    id, it, ub - lb); //bnd.second - bnd.first);
    }
  }
      
  // std::pair<double,double> ubnd = lbnd;// make sure have same solution is saturated
  
  // if( lbnd.second < 0 ){ // need a numeric solution
  //   // find upper limits of search
  //   ubnd.first = s_sz+3*vtol;
  //   z = ubnd.first;
  //   v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
  //   qq = sz->fq(z,q_sz_in);
  //   Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - qq);
  //   ubnd.second = Hz;
  //   int it(0);
  //   while( (ubnd.second < 0) and (it < max_it) ){
  //     ubnd.first *= 2.0;
  //     double z = ubnd.first;
  //     v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
  //     qq = sz->fq(z,q_sz_in);
  //     Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - qq);
  //     ubnd.second = Hz;
  //     it +=1;
  //   }
  //   if( ubnd.second < 0 ){
  //     Rcpp::warning("SZ: No bound found within %i iterations. Difference between bounds is %d.",
  // 		    it, ubnd.first - lbnd.first); //bnd.second - bnd.first);
  //   }
  //   // shrink
  //   it = 0;  
  //   while( (it <= max_it) and ( ubnd.second > vtol ) ){
  //     //z = (lbnd.first+ubnd.first)/2.0;
  //     double iW = ubnd.second / (ubnd.second-lbnd.second);
  //     iW = std::max(0.001,std::min(iW,0.999));
  //     z = (iW*lbnd.first) + (1.0-iW)*ubnd.first;
  //     v_uz_sz = area * Dt * std::min( (s_uz+v_rz_uz)/(t_d*z + area*Dt), 1/t_d );
  //     qq = sz->fq(z,q_sz_in);
  //     Hz = z - s_sz + v_uz_sz + Dt*(q_sz_in - qq);
      
  //     if( Hz < 0 ){
  // 	lbnd.first = z;
  // 	lbnd.second = Hz;
  //     }else{
  // 	ubnd.first = z;
  // 	ubnd.second = Hz;
  //     }
  //     it += 1;  
  //   }
  //   if(it > max_it){
  //     Rcpp::warning("HRU %i SZ: No solution found within %i iterations. Difference between bounds is %d",
  // 		    id, it, ubnd.first - lbnd.first); //bnd.second - bnd.first);
  //   }
  // }

  // upward pass
  q_sz = sz->fq(ub,q_sz_in);
  v_uz_sz = s_sz + Dt*(q_sz-q_sz_in) - ub;
  s_sz = ub;
  
  // q_sz = sz->fq(ubnd.first,q_sz_in);
  // v_uz_sz = s_sz + Dt*(q_sz-q_sz_in) - ubnd.first;
  // s_sz = ubnd.first;

  z = std::min(s_sz, s_uz+v_rz_uz-v_uz_sz);
  v_rz_uz = z + v_uz_sz - s_uz;
  s_uz = z;

  v_sf_rz = std::min( v_sf_rz, (area*s_rzmax) - s_rz - Dt*(precip - pet) + v_rz_uz);
  s_rz = ((area*s_rzmax) / ((area*s_rzmax) + Dt*pet)) * (s_rz + Dt*precip + v_sf_rz - v_rz_uz);
  aet = pet * s_rz / (area*s_rzmax);
  
  // surface
  
  sf->update(s_sf, q_sf, q_sf_in, v_sf_rz, Dt, vtol, max_it);
     
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

#include "Rcpp.h"
#include <boost/integer/common_factor.hpp>  
#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>
// Functions for the exponential profile

// START HERE:
// (i) Need to check output
// (ii) fix reorder warning in complitation
// (iii) Can we generalise?

// recall: Rcpp deep copies inputs unless they are the correct type!

// template class for solving the saturated zone - must have unique name
template <class T>
struct fsz_exp {
  fsz_exp(double const& uz_, double const& sz_,
	  double const& lmax_, double const& cbm_,
	  double const& td_,
	  double const& dt_, double const& dx_,
	  double const& lin_, double const& rc_):
    lmax(lmax_), cbm(cbm_),
    td(td_), dt(dt_), dx(dx_),
    rc(rc_), lin(lin_), uz(uz_), sz(sz_)
  { /* Constructor just stores values to use in solution */ }
  double operator()(double const& x){
    double l = lmax*std::exp(-x*cbm);
    double r = std::min( 1.0/td , (uz + dt*rc)/(td*x + dt) );
    return x - sz - dt*(l/dx - lin/dx - r);
  }
private:
  double lmax, cbm, dt, td, rc, lin, uz, sz, dx; // values to use in calc
};

struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 1e-32;
  }
};

//* ======================================================================================= *//
//* ======================================================================================= *//
//* ======================================================================================= *//

// Function for initialising
// [[Rcpp::export]]
void dt_exp_init(std::vector<int> id, // hillslope id
		 Rcpp::NumericMatrix states, // hillslope states
		 Rcpp::NumericMatrix attr, // hillslope attributes
		 Rcpp::NumericMatrix param, // hillslope parameters
		 std::vector<int> channel_id, // need to allow full passing
		 std::vector<int> flow_from,
		 std::vector<int> flow_to,
		 std::vector<double> flow_frc,
		 std::vector<double> r_uz_sz_0 // initial recharge as flux per unit area
		 ){
  // work out some dimensions
  int nhillslope = states.nrow();
  int nlink = flow_from.size();
  
  // create vectors for storing lateral fluxes
  int maxid = std::max( *std::max_element(std::begin(id), std::end(id)),
			*std::max_element(std::begin(channel_id),std::end(channel_id)) );
  std::vector<double> q_sf_in(maxid+1,0.0), q_sz_in(maxid+1,0.0);
  // std::vector<double> q_sf_out(maxid,0.0), q_sz_out(maxid,0.0);
  
  // seperate out states to the hillslope
  Rcpp::NumericMatrix::Column s_sf = states.column(0);
  Rcpp::NumericMatrix::Column s_rz = states.column(1);
  Rcpp::NumericMatrix::Column s_uz = states.column(2);
  Rcpp::NumericMatrix::Column s_sz = states.column(3);

  // seperate out attributes of the hillslope - transmissivity dependent?
  // Must match model description in R code
  Rcpp::NumericMatrix::Column atb_bar = attr.column(0); // average topographic index
  Rcpp::NumericMatrix::Column s_bar = attr.column(1); // average gradient
  Rcpp::NumericMatrix::Column area = attr.column(2); // surface area (plan)
  Rcpp::NumericMatrix::Column width = attr.column(3); // contour length of outflow

  // seperate out parameters of the hillslope - transmissivity dependent
  // Must match the order in model description in R code
  Rcpp::NumericMatrix::Column r_sf_max = param.column(0); // max downflow rate from surface
  Rcpp::NumericMatrix::Column s_rz_max = param.column(1); // max soil moisture depth
  Rcpp::NumericMatrix::Column s_rz_0 = param.column(2); // initial root zone depth as fraction
  Rcpp::NumericMatrix::Column c_sf = param.column(3); // surface flow celerity
  Rcpp::NumericMatrix::Column t_d = param.column(4); // unsaturated zone time constant
  Rcpp::NumericMatrix::Column ln_t0 = param.column(5); // log of saturated transmissivity
  Rcpp::NumericMatrix::Column m = param.column(6); // transmissivity decay parameter
  

  // compute the property summaries for the hillslope required
  std::vector<double> l_sz_max(nhillslope,-999.0); // max saturated zone flux
  std::vector<double> log_l_sz_max(nhillslope,-999.0); // log of max saturated zone flux
  std::vector<double> beta (nhillslope,-999.0); // slope angle
  std::vector<double> cosbeta_m (nhillslope,-999.0); // cos(beta)/m

  for(int i=0;i<nhillslope;++i){
    beta[i] = std::atan(s_bar[i]);
    l_sz_max[i] = std::exp(ln_t0[i])*std::sin(beta[i]);
    log_l_sz_max[i] = ln_t0[i] + std::log( std::sin(beta[i]) );
    cosbeta_m[i] = std::cos(beta[i]) /m[i];
  }


  
  // variable use with loop
  int cid; // current id
  double l_sz; // saturated zone outflow per unit width
  double l_sz_in; // saturated zone inflow per unit width
  double Dx; // effective length
  double r_uz_sz; // actual vertical recharge
  double q_sf_out, q_sz_out; // surface and saturated zone outflow (total)
  
  // variables for handling flow linkages
  int link_cntr = 0; // counter for flow links
  int link_from_id = flow_from[link_cntr];
  
  for(int ii=0; ii<nhillslope; ++ii){
    cid = id[ii]; // current id used to reference longer vectors
    Dx = area[ii]/width[ii]; // compute effective length
    //Rcpp::Rcout << cid << " " <<q_sz_in.size() << std::endl;
    
    l_sz_in = q_sz_in[cid] / width[ii]; // standardise inflows by width
    s_sf[ii] = 0.0; // initialise surface store
    s_rz[ii] = s_rz_max[ii]*s_rz_0[ii]; // initialise root zone store
    l_sz = std::min(l_sz_max[ii], l_sz_in + Dx*r_uz_sz_0[ii]); // outflow flux under steady state
    r_uz_sz = (l_sz - l_sz_in)/Dx; // solve to find actual r_uz_sz
    // compute saturated zone storeage deficit
    s_sz[ii] = std::max(0.0, (log_l_sz_max[ii] - std::log(l_sz))/cosbeta_m[ii]);
    s_uz[ii] = std::min( s_sz[ii], r_uz_sz*t_d[ii]*s_sz[ii] ); // compute unsaturated zone storage

    // tranfer on the outflow
    q_sf_out = width[ii]*c_sf[ii]*s_sf[ii];
    q_sz_out = width[ii]*l_sz;
    while( (link_from_id == cid) & (link_cntr<nlink) ){
      int& j = flow_to[link_cntr];
      double& f = flow_frc[link_cntr];
      q_sf_in[j] += f*q_sf_out;
      q_sz_in[j] += f*q_sz_out;
      link_cntr +=1;
      if(link_cntr < nlink) {
  	link_from_id = flow_from[link_cntr];
      }
    }
  }
}


//* ======================================================================================= *//
//* ======================================================================================= *//
//* ======================================================================================= *//

// Function for initialising
// [[Rcpp::export]]
void dt_exp_implicit(std::vector<int> id, // hillslope id
		     Rcpp::NumericMatrix states, // hillslope states
		     Rcpp::NumericMatrix attr, // hillslope attributes
		     Rcpp::NumericMatrix param, // hillslope parameters
		     std::vector<int> channel_id, // need to allow full passing
		     Rcpp::NumericMatrix channel_attr, // attributes of the channel
		     std::vector<int> flow_from, // from part of linkages
		     std::vector<int> flow_to, // to part of linkages
		     std::vector<double> flow_frc, // fraction part of linkages
		     std::vector<int> precip_col, // column number of precipitation series to use
		     std::vector<int> precip_id, // HRU id of precipitation input
		     std::vector<double> precip_frc, // fraction of precip to add
		     std::vector<int> pet_col, // column number of pet series to use
		     std::vector<int> pet_id, // HRU id of pet input
		     std::vector<double> pet_frc, // fraction of pet to add to input
		     Rcpp::NumericMatrix obs, // hillslope statesexternal series
		     Rcpp::NumericMatrix channel_inflow, // channel_inflow
		     Rcpp::NumericMatrix mass_balance, // mass balance for each timestep
		     std::vector<bool> keep_states,
		     Rcpp::List state_rec,
		     double timestep,
		     int n_sub_step
		     ){

  //Rcpp::Rcout << "Arrived in function...." << std::endl;
  // work out some dimensions
  int nhillslope = id.size();
  int nchannel = channel_id.size();
  int nlink = flow_from.size();
  
  //Rcpp::Rcout << "Set dimensions" << std::endl;
 
  
  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;
  
  //Rcpp::Rcout << "Computed timestep" << std::endl;

  // create vectors for storing lateral fluxes
  int maxid = std::max( *std::max_element(std::begin(id), std::end(id)),
			*std::max_element(std::begin(channel_id),std::end(channel_id)) );
  std::vector<double> q_sf_in(maxid+1,0.0), q_sz_in(maxid+1,0.0);
  std::vector<double> precip(maxid+1,0.0), pet(maxid+1,0.0);
  // std::vector<double> q_sf_out(maxid,0.0), q_sz_out(maxid,0.0);

  //Rcpp::Rcout << "Created lateral flux vector" << std::endl;
  
  // seperate out states to the hillslope
  Rcpp::NumericMatrix::Column s_sf = states.column(0);
  Rcpp::NumericMatrix::Column s_rz = states.column(1);
  Rcpp::NumericMatrix::Column s_uz = states.column(2);
  Rcpp::NumericMatrix::Column s_sz = states.column(3);

  //Rcpp::Rcout << "Seperated states" << std::endl;
  
  // seperate out attributes of the hillslope - transmissivity dependent?
  // Must match model description in R code
  Rcpp::NumericMatrix::Column atb_bar = attr.column(0); // average topographic index
  Rcpp::NumericMatrix::Column s_bar = attr.column(1); // average gradient
  Rcpp::NumericMatrix::Column area = attr.column(2); // surface area (plan)
  Rcpp::NumericMatrix::Column width = attr.column(3); // contour length of outflow

  //Rcpp::Rcout << "Seperated hillslope attributes" << std::endl;
  
  // seperate out parameters of the hillslope - transmissivity dependent
  // Must match the order in model description in R code
  Rcpp::NumericMatrix::Column r_sf_max = param.column(0); // max downflow rate from surface
  Rcpp::NumericMatrix::Column s_rz_max = param.column(1); // max soil moisture depth
  Rcpp::NumericMatrix::Column s_rz_0 = param.column(2); // initial root zone depth as fraction
  Rcpp::NumericMatrix::Column c_sf = param.column(3); // surface flow celerity
  Rcpp::NumericMatrix::Column t_d = param.column(4); // unsaturated zone time constant
  Rcpp::NumericMatrix::Column ln_t0 = param.column(5); // log of saturated transmissivity
  Rcpp::NumericMatrix::Column m = param.column(6); // transmissivity decay parameter
 

  //Rcpp::Rcout << r_sf_max[1] << std::endl;
  //Rcpp::Rcout << s_rz_max[1] << std::endl;
  //Rcpp::Rcout << s_rz_0[1] << std::endl;
  //Rcpp::Rcout << c_sf[1] << std::endl;
  //Rcpp::Rcout << ln_t0[1] << std::endl;
  //Rcpp::Rcout << m[1] << std::endl;
  //Rcpp::Rcout << t_d[1] << std::endl;
  
  //Rcpp::Rcout << "Seperated hillslope param" << std::endl;
  
  // seperate out the attributes of the channel
  Rcpp::NumericMatrix::Column channel_area = channel_attr.column(0); // channel surface area

  //Rcpp::Rcout << "Seperated channel attributes" << std::endl;
  
  //Rcpp::Rcout << "End of unpacking the columns" << std::endl;
  
  // compute the property summaries for the hillslope required
  std::vector<double> l_sz_max(nhillslope,-999.0); // max saturated zone flux
  std::vector<double> beta (nhillslope,-999.0); // slope angle
  std::vector<double> cosbeta_m (nhillslope,-999.0); // cos(beta)/m
  std::vector<double> Dx (nhillslope,-999.0); // Effective length

  for(int i=0;i<nhillslope;++i){
    beta[i] = std::atan(s_bar[i]);
    l_sz_max[i] = std::exp(ln_t0[i])*std::sin(beta[i]);
    cosbeta_m[i] = std::cos(beta[i]) /m[i];
    Dx[i] = area[i]/width[i];
  }

  //Rcpp::Rcout << "End of calculation of summary variables" << std::endl;
  
  // boost optimisation parameters
  const boost::uintmax_t opt_maxit = 1000;
  int digits = std::numeric_limits<double>::digits;
  int get_digits = digits-1;
  // get_digits - highest tolerence (most accurate) is digits-1,
  // but this can be problematic with numeric rounding in
  // the functions being optimised
  // To low a tolerence, for example (digits-3)/4. produces significant mass balance errors
  boost::math::tools::eps_tolerance<double> tol(get_digits);
  boost::uintmax_t opt_it = opt_maxit;
  std::pair<double, double> opt_res; // solution for output

  //Rcpp::Rcout << "Set optimisation parameters" << std::endl;
  
  // variable use with loop
  int cid; // current id
  double r_sf_rz,r_rz_uz; // fluxes between stores
  double l_sz_in; // influxes per unit width
  double l_sz; // outflow per unit width
  double q_sf_out,q_sz_out; // out flow fluxes
  double chn_in; // flow volume to channel
  std::vector<double> mbv(5,0.0); // mass balance vector
  Rcpp::Rcout << "nchannel " << nchannel << std::endl;
  std::vector<double> ch_in(nchannel,0.0); // channel inflow vector

  // variables for handling links
  int link_cntr = 0; // counter for flow links
  int link_from_id = flow_from[link_cntr];

  // start loop of time steps
  for(int it = 0; it < obs.nrow(); ++it) {

    // initialise mass balance for the timestep
    std::fill(mbv.begin(), mbv.end(), 0.0);
    for(int ii=0; ii<nhillslope; ++ii){
      mbv[0] += area[ii]*(s_sf[ii] + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }

    // initialise channel inflow
    std::fill(ch_in.begin(), ch_in.end(), 0.0);

    // compute the precipitation input
    std::fill(precip.begin(), precip.end(),0.0);
    for(uint ii=0; ii<precip_id.size(); ++ii){
      int& i = precip_id[ii];
      int& c = precip_col[ii];
      double& f = precip_frc[ii];
      Rcpp::Rcout << i << std::endl;
      Rcpp::Rcout << c << std::endl;
      Rcpp::Rcout << f << std::endl;
      Rcpp::Rcout << obs(it,c) << std::endl;
      Rcpp::Rcout << timestep << std::endl;
      precip[i] += f*obs(it,c)/timestep;
    }

    // compute the pet input
    std::fill(pet.begin(), pet.end(),0.0);
    for(uint ii=0; ii<pet_id.size(); ++ii){
      int& i = pet_id[ii];
      int& c = pet_col[ii];
      double& f = pet_frc[ii];
      pet[i] += f*obs(it,c)/timestep;
    }
   
    //Rcpp::Rcout << "Computed precip + pet" << std::endl;
    
    // start loop of substeps
    for(int nn = 0; nn < n_sub_step; ++nn){

      // set flow passing records to 0
      std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
      std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);

      // set all counters to initial value
      link_cntr = 0; // counter for flow links
      link_from_id = flow_from[link_cntr];
      
      // loop HSUs
      for(int ii=0; ii<nhillslope; ++ii){

    	//Rcpp::Rcout << it << " " << nn << " " << ii << std::endl;
    	// current id used to reference longer vectors
    	cid = id[ii];

    	//Rcpp::Rcout << cid << " " << q_sf_in.size() << " " << q_sz_in.size() << " " << precip.size() << " " << pet.size() << std::endl;
    	// standardise inflows by width
    	// l_sf_in = q_sf_in[cid] / w[ii];
    	l_sz_in = q_sz_in[cid] / width[ii];

    	// compute first downward flux estimate from surface ans root zone
    	// these are given as \check{r} in documentation	
    	r_sf_rz = std::min( r_sf_max[ii] , (s_sf[ii] + Dt*q_sf_in[cid]/area[ii])/Dt );
    	r_rz_uz = std::max( 0.0 ,
    			    (s_rz[ii] + Dt*(precip[cid] + r_sf_rz - pet[cid]) - s_rz_max[ii])/Dt);

    	//Rcpp::Rcout << "Downward flux complete" << std::endl;
	// setup template function
	//Rcpp::Rcout << "s_uz " << s_uz[ii] << std::endl;
	//Rcpp::Rcout << "s_sz " << s_sz[ii] << std::endl;
	//Rcpp::Rcout << "l_sz_max " << l_sz_max[ii] << std::endl;
	//Rcpp::Rcout << "cosbeta_m " << cosbeta_m[ii] << std::endl;
	//Rcpp::Rcout << "t_d " << t_d[ii] << std::endl;
	//Rcpp::Rcout << "Dt " << Dt << std::endl;
	//Rcpp::Rcout << "Dx " << Dx[ii] << std::endl;
	//Rcpp::Rcout << "l_sz_in " << l_sz_in << std::endl;
	//Rcpp::Rcout << "r_rz_uz " << r_rz_uz << std::endl;
    	fsz_exp<double> fnc(s_uz[ii], s_sz[ii],
    			    l_sz_max[ii], cosbeta_m[ii],
    			    t_d[ii], Dt,Dx[ii],  l_sz_in, r_rz_uz);
	
    	//Rcpp::Rcout << "Set template funcion" << std::endl;
	
    	//double tmp = 0.0;
    	//tmp = fnc(0.0);
    	//Rcpp::Rcout << tmp << std::endl;
	
    	// test for saturation
    	if( fnc(0.0) > 0.0 ){
    	  // then saturated
	  //Rcpp::Rcout << "saturated" << std::endl;
	  l_sz = l_sz_max[ii];
	  r_rz_uz = ( s_sz[ii] + (Dt/Dx[ii])*(l_sz - l_sz_in) - s_uz[ii] )/Dt;
	  s_sz[ii] = 0.0;
    	}else{
	  //Rcpp::Rcout << "unsaturated" << std::endl;
	  
	  
    	  //tmp = fnc(s_sz[ii]);
    	  //Rcpp::Rcout << "s_sz[ii] " << s_sz[ii] << " fnc(s_sz[ii]) " << tmp << std::endl;
    	  // not saturated need to solve
    	  if( fnc(s_sz[ii]) >= 0.0 ){
    	    Rcpp::Rcout << "wetting" << std::endl;
    	    // then wetting - solution between current s_sz and 0
	    opt_res = boost::math::tools::bisect(fnc, 0.0, s_sz[ii], TerminationCondition(),opt_it);
    	    //opt_res = boost::math::tools::toms748_solve(fnc, 0.0, s_sz[ii], TerminationCondition(),opt_it);
	    
    	  }else{
    	    Rcpp::Rcout << "drying" << std::endl;
    	    // drying s_sz getting bigger
    	    double upr = 10.0; //2.0*s_sz[ii];
	    //tmp = fnc(upr);
	    //Rcpp::Rcout << "upper: " << upr << " fnc(upr) "<< tmp << std::endl;
    	    //while( fnc(upr) <= 0.0 ){
    	    //  upr += upr;
    	    //}
	    opt_res = boost::math::tools::bisect(fnc, s_sz[ii],upr, TerminationCondition(),opt_it);
    	    //opt_res = boost::math::tools::toms748_solve(fnc, s_sz[ii], upr, TerminationCondition(),opt_it);
    	  }
    	  if(opt_it >= opt_maxit){
    	    Rcpp::Rcout << "Unable to locate solution in chosen iterations:" <<
    	      " Current best guess is between " << opt_res.first << " and " <<
    	      opt_res.second << std::endl;
    	  }
	  
	  // l_sz = l_sz_max[ii];
	  // r_rz_uz = ( s_sz[ii] + (Dt/Dx[ii])*(l_sz - l_sz_in) - s_uz[ii] )/Dt;
	  // s_sz[ii] = 0.0;
    	  s_sz[ii] = (opt_res.second + opt_res.first)/2.0;
    	  r_rz_uz = std::min( r_rz_uz, (s_sz[ii] + Dt/t_d[ii] - s_uz[ii])/Dt );
    	  l_sz = l_sz_max[ii]*exp(-s_sz[ii]*cosbeta_m[ii]);
    	}
	
    	//Rcpp::Rcout << "computed optimisation" << std::endl;
	
    	// solve unsaturated zone
    	s_uz[ii] = ( (t_d[ii]*s_sz[ii])/(t_d[ii]*s_sz[ii] + Dt) ) * (s_uz[ii]+Dt*r_rz_uz);
    	// compute revised r_sf_rz
    	r_sf_rz = std::min( r_sf_rz,
    			    (s_rz_max[ii] - s_rz[ii] - Dt*(precip[cid]-pet[cid]-r_rz_uz))/Dt
    			    );

	Rcpp::Rcout << "step HRU " << it << " " << ii << std::endl;
	Rcpp::Rcout << "r_rz_uz " << r_rz_uz << std::endl;
	Rcpp::Rcout << "r_sf_rz " << r_sf_rz << std::endl;
	Rcpp::Rcout << "l_sz_in " << l_sz_in << std::endl;
	Rcpp::Rcout << "l_sz " << l_sz << std::endl;
	double tmp = fnc(s_sz[ii]);
	Rcpp::Rcout << "fnc " << tmp << std::endl;
	Rcpp::Rcout << "iterations " << opt_it << std::endl;
	// solve for root zone
	s_rz[ii] = ( s_rz[ii] + Dt*(precip[cid] + r_sf_rz - r_rz_uz) ) * ( s_rz_max[ii] / (s_rz_max[ii] + Dt*pet[ii]) );
    	// solve for surface
    	s_sf[ii] = ( Dx[ii] / (Dx[ii] + Dt*c_sf[ii]) ) * ( s_sf[ii] + Dt*( (q_sf_in[cid]/area[ii]) - r_sf_rz ));
	
    	//Rcpp::Rcout << "Solved for stores" << std::endl;
	
    	// apply pet loss and precip input to mass balance
    	mbv[1] -= pet[ii]*(s_rz[ii]/s_rz_max[ii])*area[ii]*Dt;
    	mbv[2] += precip[ii]*area[ii]*Dt;

    	//Rcpp::Rcout << "Adjusted mass balance" << std::endl;
	
    	// tranfer on the outflow
    	q_sf_out = width[ii]*c_sf[ii]*s_sf[ii];
    	q_sz_out = width[ii]*l_sz;
	//Rcpp::Rcout << link_cntr << " " << link_from_id << std::endl;
    	while( (link_from_id == cid) & (link_cntr<nlink) ){
    	  int j = flow_to[link_cntr];
    	  double f = flow_frc[link_cntr];
    	  q_sf_in[j] += f*q_sf_out; //w[ii]*l_sf[ii];
    	  q_sz_in[j] += f*q_sz_out; //w[ii]*l_sz[ii];
    	  link_cntr +=1;
    	  if(link_cntr < nlink) {
    	    link_from_id = flow_from[link_cntr];
    	  }
    	}
	//Rcpp::Rcout << link_cntr << " " << link_from_id << std::endl;
    	// end of hillslope loop
    	
      }
      
      // loop channels for the sub step
      for(int ii=0; ii < nchannel; ++ii){
      	cid = channel_id[ii];
      	chn_in = (channel_area[ii]*precip[cid]+ q_sf_in[cid] + q_sz_in[cid])*Dt;
      	mbv[2] += precip[cid]*channel_area[ii]*Dt; 
      	ch_in[ii] += chn_in; // volume of flow to channel
      	mbv[3] -= chn_in; // volume lost from hillslope to channel
      }
      // end of substep loop
    }
    //Rcpp::Rcout << "after substep" << std::endl;
    // final mass balance states
    for(int ii=0; ii<nhillslope; ++ii){
      mbv[4] -= area[ii]*(s_sf[ii] + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }

    
    // copy mass balance to record
    for(int ii=0; ii < 5; ++ii){
      mass_balance(it,ii) = mbv[ii];
    }
    //Rcpp::Rcout << "Copied MB record" << std::endl;
    
    // convert channel inflow to rate
    for(int ii=0; ii < nchannel; ++ii){
      channel_inflow(it,ii) = ch_in[ii]/timestep; //channel_inflow(it,ii)/timestep;
    }
    //Rcpp::Rcout << "Copied channel_inflow" << std::endl;
    
    // keep states if required
    if( keep_states[it] ){
      state_rec(it) = Rcpp::clone(states);
    }
    //Rcpp::Rcout << "Done keep states" << std::endl;
    // check user interupt
    Rcpp::checkUserInterrupt(); 
  }
}


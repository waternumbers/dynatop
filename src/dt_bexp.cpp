#include "Rcpp.h"
#include <boost/math/tools/roots.hpp>

// Functions for the bounded (maximumd epth) exponential profile

// recall: Rcpp deep copies inputs unless they are the correct type!

// template class for solving the saturated zone - must have unique name
template <class T>
struct fsz_bexp {
  fsz_bexp(double const& uz_, double const& sz_,
	   double const& lsc_, double const& cbm_,
	   double const lcnst_, double const& td_,
	   double const& dt_, double const& dx_,
	   double const& lin_, double const& rc_):
    lsc(lsc_), cbm(cbm_), lcnst(lcnst_),
    td(td_), dt(dt_), dx(dx_),
    rc(rc_), lin(lin_), uz(uz_), sz(sz_)
  { /* Constructor just stores values to use in solution */ }
  double operator()(double const& x){
    double l = lsc * ( std::exp(-x*cbm) - lcnst );
    double r = std::min( 1.0/td , (uz + dt*rc)/(td*x + dt) );
    return x - sz - dt*(l/dx - lin/dx - r);
  }
private:
  double lsc, cbm, lcnst, td, dt, dx, rc, lin, uz, sz; // values to use in calc
};

struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 1e-16;
  }
};

//* ======================================================================================= *//
//* ======================================================================================= *//
//* ======================================================================================= *//

// Function ofr computing the courant numbers for the surface and unsaturated zones
// [[Rcpp::export]]
void dt_bexp_courant(Rcpp::DataFrame hillslope, // hillslope data frame
		    Rcpp::NumericMatrix courant, // courant numbers
		    double timestep, // time step
		    int n_sub_step // number of sub steps
		    ){
  // seperate out attributes of the hillslope
  // recall: NumericVector are references to the DataFrame
  // Rcpp::NumericMatrix::Column atb_bar = attr.column(0); // average topographic index
  Rcpp::IntegerVector id = hillslope["id"]; 
  Rcpp::NumericVector s_bar = hillslope["s_bar"]; // average gradient
  Rcpp::NumericVector area = hillslope["area"]; // surface area (plan)
  Rcpp::NumericVector width = hillslope["width"]; // contour length of outflow
  Rcpp::NumericVector c_sf = hillslope["c_sf"]; // surface flow celerity
  Rcpp::NumericVector ln_t0 = hillslope["ln_t0"]; // log of saturated transmissivity
  Rcpp::NumericVector m = hillslope["m"]; // transmissivity decay parameter

  // work out some dimensions
  int nhillslope = id.size();
  
  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;

  double c_sz_max(0.0), beta(0.0);
  
  for(int i=0;i<nhillslope;++i){
    // surface number
    courant(i,0) = c_sf[i]*Dt*width[i]/area[i];
    beta = std::atan(s_bar[i]);
    c_sz_max = std::exp(ln_t0[i])*std::sin(2.0*beta)/(2.0*m[i]);
    courant(i,1) = c_sz_max*Dt*width[i]/area[i];
  }
}
  
//* ======================================================================================= *//
//* ======================================================================================= *//
//* ======================================================================================= *//

// Function for initialising
// [[Rcpp::export]]
void dt_bexp_init(Rcpp::DataFrame hillslope, // hillslope data frame
		 Rcpp::DataFrame channel, // channel data frame
		 Rcpp::DataFrame flow_direction, // flow directions data frame
		 std::vector<double> r_uz_sz_0 // initial recharge as flux per unit area
		 ){

  // seperate out hillslope to vector
  // recall NumericVector are references in this case
  Rcpp::IntegerVector id = hillslope["id"];
  Rcpp::NumericVector s_sf = hillslope["s_sf"];
  Rcpp::NumericVector s_rz = hillslope["s_rz"];
  Rcpp::NumericVector s_uz = hillslope["s_uz"];
  Rcpp::NumericVector s_sz = hillslope["s_sz"];
  Rcpp::NumericVector s_bar = hillslope["s_bar"]; // average gradient
  Rcpp::NumericVector area = hillslope["area"]; // surface area (plan)
  Rcpp::NumericVector width = hillslope["width"]; // contour length of outflow
  Rcpp::NumericVector s_rz_max = hillslope["s_rzmax"]; // max soil moisture depth
  Rcpp::NumericVector s_rz_0 = hillslope["s_rz0"]; // initial root zone depth as fraction
  Rcpp::NumericVector c_sf = hillslope["c_sf"]; // surface flow celerity
  Rcpp::NumericVector t_d = hillslope["t_d"]; // unsaturated zone time constant
  Rcpp::NumericVector ln_t0 = hillslope["ln_t0"]; // log of saturated transmissivity
  Rcpp::NumericVector m = hillslope["m"]; // transmissivity decay parameter
  Rcpp::NumericVector D = hillslope["D"]; // maximum storage depth

  // seperate out channel to vector
  // recall NumericVector are references in this case
  Rcpp::IntegerVector channel_id = channel["id"];

  // seperate out the flow directions
  Rcpp::IntegerVector flow_from = flow_direction["from"];
  Rcpp::IntegerVector flow_to = flow_direction["to"];
  Rcpp::NumericVector flow_frc = flow_direction["frc"];
  
  // work out some dimensions
  int nhillslope = id.size();
  int nlink = flow_from.size();
  int maxid = std::max( *std::max_element(std::begin(id), std::end(id)),
			*std::max_element(std::begin(channel_id),std::end(channel_id)) );

  // create vectors for passing flux
  // +1 since R indexs from 1 and C++ 0
  std::vector<double> q_sf_in(maxid+1,0.0), q_sz_in(maxid+1,0.0);
    
  // compute the property summaries for the hillslope required
  std::vector<double> l_sc(nhillslope,-999.0); // saturated zone flux scale valie
  std::vector<double> l_cnst(nhillslope,-999.0); // saturated zone constant offset
  std::vector<double> l_sz_max(nhillslope,-999.0); // max saturated zone flux
  std::vector<double> log_l_sz_max(nhillslope,-999.0); // log of max saturated zone flux
  std::vector<double> beta (nhillslope,-999.0); // slope angle
  std::vector<double> cosbeta_m (nhillslope,-999.0); // cos(beta)/m
  
  for(int i=0;i<nhillslope;++i){
    beta[i] = std::atan(s_bar[i]);
    cosbeta_m[i] = std::cos(beta[i]) /m[i];
    l_sc[i] = std::exp(ln_t0[i])*std::sin(beta[i]);
    l_cnst[i] = std::exp( -D[i] *cosbeta_m[i] );
    l_sz_max[i] = l_sc[i] * (1 - l_cnst[i]);
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
    l_sz_in = q_sz_in[cid] / width[ii]; // standardise inflows by width
    s_sf[ii] = 0.0; // initialise surface store
    s_rz[ii] = s_rz_max[ii]*s_rz_0[ii]; // initialise root zone store
    l_sz = std::min(l_sz_max[ii], l_sz_in + Dx*r_uz_sz_0[ii]); // outflow flux under steady state
    r_uz_sz = (l_sz - l_sz_in)/Dx; // solve to find actual r_uz_sz
    // compute saturated zone storeage deficit
    s_sz[ii] = std::max(0.0, -std::log( (l_sz/l_sc[ii]) + l_cnst[ii] ) / cosbeta_m[ii] );
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

// Function for solving
// [[Rcpp::export]]
void dt_bexp_implicit(Rcpp::DataFrame hillslope, // hillslope data frame
		     Rcpp::DataFrame channel, // channel data frame
		     Rcpp::DataFrame flow_direction, // flow directions data frame
		     Rcpp::DataFrame precip_input, // precipitation input data frame
		     Rcpp::DataFrame pet_input, // PET input data frame
		     Rcpp::NumericMatrix obs, // external series
		     Rcpp::NumericMatrix channel_inflow, // channel_inflow to compute
		     Rcpp::NumericMatrix mass_balance, // mass balance for each timestep
		     std::vector<bool> keep_states,
		     Rcpp::List state_rec,
		     double timestep,
		     int n_sub_step
		     ){

  // seperate out hillslope to vector
  // recall NumericVector are references in this case
  Rcpp::IntegerVector id = hillslope["id"];
  Rcpp::NumericVector s_sf = hillslope["s_sf"];
  Rcpp::NumericVector s_rz = hillslope["s_rz"];
  Rcpp::NumericVector s_uz = hillslope["s_uz"];
  Rcpp::NumericVector s_sz = hillslope["s_sz"];
  Rcpp::NumericVector s_bar = hillslope["s_bar"]; // average gradient
  Rcpp::NumericVector area = hillslope["area"]; // surface area (plan)
  Rcpp::NumericVector width = hillslope["width"]; // contour length of outflow
  Rcpp::NumericVector r_sf_max = hillslope["r_sfmax"]; // max flux down from surface
  Rcpp::NumericVector s_rz_max = hillslope["s_rzmax"]; // max soil moisture depth
  Rcpp::NumericVector c_sf = hillslope["c_sf"]; // surface flow celerity
  Rcpp::NumericVector t_d = hillslope["t_d"]; // unsaturated zone time constant
  Rcpp::NumericVector ln_t0 = hillslope["ln_t0"]; // log of saturated transmissivity
  Rcpp::NumericVector m = hillslope["m"]; // transmissivity decay parameter
  Rcpp::NumericVector D = hillslope["D"]; // maximum storage deficit
 
  // rebuild the states as a data frame for saving
  Rcpp::DataFrame states =
    Rcpp::DataFrame::create(Rcpp::Named("id") = id ,
			    Rcpp::Named("s_sf") = s_sf,
			    Rcpp::Named("s_rz") = s_rz,
			    Rcpp::Named("s_uz") = s_uz,
			    Rcpp::Named("s_sz") = s_sz);
  
  // seperate out channel to vector
  // recall NumericVector are references in this case
  Rcpp::IntegerVector channel_id = channel["id"];
  Rcpp::NumericVector channel_area = channel["area"]; // channel surface area

  // seperate out the flow directions
  Rcpp::IntegerVector flow_from = flow_direction["from"];
  Rcpp::IntegerVector flow_to = flow_direction["to"];
  Rcpp::NumericVector flow_frc = flow_direction["frc"];

  // seperate out the precip_input
  Rcpp::IntegerVector precip_id = precip_input["id"];
  Rcpp::IntegerVector precip_col = precip_input["col_idx"];
  Rcpp::NumericVector precip_frc = precip_input["frc"];

  // seperate out the pet_input
  Rcpp::IntegerVector pet_id = pet_input["id"];
  Rcpp::IntegerVector pet_col = pet_input["col_idx"];
  Rcpp::NumericVector pet_frc = pet_input["frc"];


  // work out some dimensions
  int nhillslope = id.size();
  int nchannel = channel_id.size();
  int nlink = flow_from.size();
  int maxid = std::max( *std::max_element(std::begin(id), std::end(id)),
			*std::max_element(std::begin(channel_id),std::end(channel_id)) );
   
  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;

  // create vectors for storing lateral fluxes, precip and pet
  std::vector<double> q_sf_in(maxid+1,0.0), q_sz_in(maxid+1,0.0);
  std::vector<double> precip(maxid+1,0.0), pet(maxid+1,0.0);

  // compute the property summaries for the hillslope required
  std::vector<double> l_sc(nhillslope,-999.0); // saturated zone flux scale valie
  std::vector<double> l_cnst(nhillslope,-999.0); // saturated zone constant offset
  std::vector<double> beta (nhillslope,-999.0); // slope angle
  std::vector<double> cosbeta_m (nhillslope,-999.0); // cos(beta)/m
  std::vector<double> Dx (nhillslope,-999.0); // Effective length
  
  for(int i=0;i<nhillslope;++i){
    beta[i] = std::atan(s_bar[i]);
    cosbeta_m[i] = std::cos(beta[i]) /m[i];
    l_sc[i] = std::exp(ln_t0[i])*std::sin(beta[i]);
    l_cnst[i] = std::exp( -D[i] *cosbeta_m[i] );
    Dx[i] = area[i]/width[i];
  }
  
  // boost optimisation parameters
  const boost::uintmax_t opt_maxit = 1000;
  boost::uintmax_t opt_it = opt_maxit;
  std::pair<double, double> opt_res; // solution for output
  
  // variable use with loop
  int cid; // current id
  double r_sf_rz,r_rz_uz; // fluxes between stores
  double l_sz_in; // influxes per unit width
  double l_sz; // outflow per unit width
  double q_sf_out,q_sz_out; // out flow fluxes
  double chn_in; // flow volume to channel
  std::vector<double> mbv(5,0.0); // mass balance vector
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
    for(unsigned int ii=0; ii<precip_id.size(); ++ii){
      int& i = precip_id[ii];
      int& c = precip_col[ii];
      double& f = precip_frc[ii];
      precip[i] += f*obs(it,c)/timestep;
    }

    // compute the pet input
    std::fill(pet.begin(), pet.end(),0.0);
    for(unsigned int ii=0; ii<pet_id.size(); ++ii){
      int& i = pet_id[ii];
      int& c = pet_col[ii];
      double& f = pet_frc[ii];
      pet[i] += f*obs(it,c)/timestep;
    }
      
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

	// setup template function
    	fsz_bexp<double> fnc(s_uz[ii], s_sz[ii],
			     l_sc[ii], cosbeta_m[ii],l_cnst[ii],
			     t_d[ii], Dt,Dx[ii],  l_sz_in, r_rz_uz);
	
    	// test for saturation
    	if( fnc(0.0) >= 0.0 ){
    	  // then saturated
	  l_sz = l_sc[ii] * (1 - l_cnst[ii]);
	  r_rz_uz = ( s_sz[ii] + (Dt/Dx[ii])*(l_sz - l_sz_in) - s_uz[ii] )/Dt;
	  s_sz[ii] = 0.0;
    	}else{
	  // not saturated test to see if wetting or drying
    	  if( fnc(s_sz[ii]) >= 0.0 ){
    	    // then wetting - solution between current s_sz and 0
	    opt_res = boost::math::tools::bisect(fnc, 0.0, s_sz[ii], TerminationCondition(),opt_it);
    	  }else{
	    // drying - need to work lower depth to bracket
	    double upr = 2.0*s_sz[ii];
	    double fupr = fnc(upr);
	    while( (fupr < 0.0) & (upr < D[ii])){
	      upr += upr;
	      fupr = fnc(upr);
	    }
	    opt_res = boost::math::tools::bisect(fnc, s_sz[ii],upr, TerminationCondition(),opt_it);
    	  }
    	  if(opt_it >= opt_maxit){
    	    Rcpp::Rcout << "Unable to locate solution in chosen iterations:" <<
    	      " Current best guess is between " << opt_res.first << " and " <<
    	      opt_res.second << std::endl;
    	  }
	  
    	  s_sz[ii] = (opt_res.second + opt_res.first)/2.0;
    	  r_rz_uz = std::min( r_rz_uz, (s_sz[ii] + Dt/t_d[ii] - s_uz[ii])/Dt );
    	  l_sz = l_sc[ii]*( exp(-s_sz[ii]*cosbeta_m[ii]) - l_cnst[ii]);
    	}
	
    	// solve unsaturated zone
    	s_uz[ii] = ( (t_d[ii]*s_sz[ii])/(t_d[ii]*s_sz[ii] + Dt) ) * (s_uz[ii]+Dt*r_rz_uz);
    	// compute revised r_sf_rz
    	r_sf_rz = std::min( r_sf_rz,
    			    (s_rz_max[ii] - s_rz[ii] - Dt*(precip[cid]-pet[cid]-r_rz_uz))/Dt
    			    );

	// solve for root zone
	s_rz[ii] = ( s_rz[ii] + Dt*(precip[cid] + r_sf_rz - r_rz_uz) ) * ( s_rz_max[ii] / (s_rz_max[ii] + Dt*pet[cid]) );
    	// solve for surface
    	s_sf[ii] = ( Dx[ii] / (Dx[ii] + Dt*c_sf[ii]) ) * ( s_sf[ii] + Dt*( (q_sf_in[cid]/area[ii]) - r_sf_rz ));
	
    	// apply pet loss and precip input to mass balance
    	mbv[1] -= pet[cid]*(s_rz[ii]/s_rz_max[ii])*area[ii]*Dt;
    	mbv[2] += precip[cid]*area[ii]*Dt;

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
    // final mass balance states
    for(int ii=0; ii<nhillslope; ++ii){
      mbv[4] -= area[ii]*(s_sf[ii] + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }

    // copy mass balance to record
    for(int ii=0; ii < 5; ++ii){
      mass_balance(it,ii) = mbv[ii];
    }
    
    // convert channel inflow to rate
    for(int ii=0; ii < nchannel; ++ii){
      channel_inflow(it,ii) = ch_in[ii]/timestep; //channel_inflow(it,ii)/timestep;
    }
    
    // keep states if required
    if( keep_states[it] ){
      state_rec(it) = Rcpp::clone(states);
    }

    // check user interupt
    Rcpp::checkUserInterrupt(); 
  }
}


#include "Rcpp.h"
#include "hillslope_hru.h"


// [[Rcpp::export]]
void dt_implicit_sim(Rcpp::DataFrame hillslope, // hillslope data frame
		     Rcpp::DataFrame channel, // channel data frame
		     Rcpp::DataFrame flow_direction, // flow directions data frame
		     Rcpp::DataFrame precip_input, // precipitation input data frame
		     Rcpp::DataFrame pet_input, // PET input data frame
		     Rcpp::NumericMatrix obs, // external series
		     Rcpp::NumericMatrix channel_inflow_sf, // channel_inflow from surface - to compute
		     Rcpp::NumericMatrix channel_inflow_sz, // channel_inflow from saturated - to compute
		     Rcpp::NumericMatrix mass_balance, // mass balance for each timestep
		     std::vector<bool> keep_states,
		     Rcpp::List state_rec,
		     double timestep,
		     int n_sub_step,
		     double tol,
		     int max_it,
		     double ftol
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
  Rcpp::NumericVector c_sz = hillslope["c_sz"]; // constant celerity of saturated zone
  Rcpp::NumericVector m = hillslope["m"]; // decay parameter of transmissivity profile
  Rcpp::NumericVector D = hillslope["D"]; // depth parameter of transmissivity profile
  Rcpp::NumericVector m_2 = hillslope["m_2"]; // second transmissivity decay parameter
  Rcpp::NumericVector omega = hillslope["omega"]; // second transmissivity decay parameter
  Rcpp::IntegerVector opt = hillslope["opt"]; // type of saturated zone
  
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
  std::vector<double> q_sf_out(maxid+1,0.0), q_sz_out(maxid+1,0.0);
  std::vector<double> e_a(maxid+1,0.0);
  std::vector<double> precip(maxid+1,0.0), pet(maxid+1,0.0);

  // variable use with loops
  int cid; // current id
  std::vector<double> mbv(6,0.0); // mass balance vector
  std::vector<double> ch_in_sf(nchannel,0.0); // channel inflow vector for surface
  std::vector<double> ch_in_sz(nchannel,0.0); // channel inflow vector for saturated
 
  // loop to create hillslope class objects
  std::vector<hillslope_hru> hs_hru;
  for(int ii=0; ii<nhillslope; ++ii){
    cid = id[ii];
    hs_hru.push_back(hillslope_hru(id[ii],
			       s_sf[ii], s_rz[ii], s_uz[ii], s_sz[ii],
			       s_bar[ii],   area[ii],   width[ii],
			       q_sf_in[cid], q_sf_out[cid], // surface zone lateral fluxes
			       q_sz_in[cid], q_sz_out[cid], // saturated zone lateral fluxes
			       e_a[cid], // actual evapotranspiration as a rate [m/s]
			       r_sf_max[ii],   c_sf[ii], // surface store parameters
			       s_rz_max[ii], // root zone store parameters
			       t_d[ii], // unsaturated zone parameters
			       ln_t0[ii], c_sz[ii], m[ii], D[ii], m_2[ii], omega[ii],// saturated zone parameters
			       opt[ii]   ) //type of saturated zone
		     );
  }

  // loop time steps

  // variables for handling links
  int link_cntr = 0; // counter for flow links
  int link_from_id = flow_from[link_cntr];

  // start loop of time steps
  for(int it = 0; it < obs.nrow(); ++it) {
    //Rcpp::Rcout << "Interation " << it << std::endl;

    // initialise mass balance for the timestep
    std::fill(mbv.begin(), mbv.end(), 0.0);
    for(int ii=0; ii<nhillslope; ++ii){
      mbv[0] += area[ii]*(s_sf[ii] + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }

    // initialise channel inflow
    std::fill(ch_in_sf.begin(), ch_in_sf.end(), 0.0);
    std::fill(ch_in_sz.begin(), ch_in_sz.end(), 0.0);
    
    // compute the precipitation input
    std::fill(precip.begin(), precip.end(),0.0);
    for(unsigned int ii=0; ii<precip_id.size(); ++ii){
      int& i = precip_id[ii];
      int& c = precip_col[ii];
      double& f = precip_frc[ii];
      precip[i] += f*obs(it,c)/timestep; // precip as rate
    }

    // compute the pet input
    std::fill(pet.begin(), pet.end(),0.0);
    for(unsigned int ii=0; ii<pet_id.size(); ++ii){
      int& i = pet_id[ii];
      int& c = pet_col[ii];
      double& f = pet_frc[ii];
      pet[i] += f*obs(it,c)/timestep; // pet as rate
    }
      
    // start loop of substeps
    for(int nn = 0; nn < n_sub_step; ++nn){

      // set flow passing records to 0
      std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
      std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);
      std::fill(q_sf_out.begin(), q_sf_out.end(), 0.0);
      std::fill(q_sz_out.begin(), q_sz_out.end(), 0.0);
      std::fill(e_a.begin(), e_a.end(), 0.0);

      // set all counters to initial value
      link_cntr = 0; // counter for flow links
      link_from_id = flow_from[link_cntr];
      
      // loop HSUs
      for(int ii=0; ii<nhillslope; ++ii){

	// current id used to reference longer vectors
    	cid = id[ii];

	// evolve hillslope
	//int max_it = 10000;
	hs_hru[ii].implicit_step(pet[cid], precip[cid], Dt, tol, max_it, ftol);

	
    	while( (link_from_id == cid) & (link_cntr<nlink) ){
    	  int j = flow_to[link_cntr];
    	  double f = flow_frc[link_cntr];
    	  q_sf_in[j] += f*q_sf_out[cid]; //w[ii]*l_sf[ii];
    	  q_sz_in[j] += f*q_sz_out[cid]; //w[ii]*l_sz[ii];
    	  link_cntr +=1;
    	  if(link_cntr < nlink) {
    	    link_from_id = flow_from[link_cntr];
    	  }
    	}
	
	// apply pet loss and precip input to mass balance
    	mbv[1] += e_a[cid]*area[ii]*Dt;
    	mbv[2] += precip[cid]*area[ii]*Dt;
	
	// end of hillslope loop
      }
      
      
      // loop channels for the sub step
      for(int ii=0; ii < nchannel; ++ii){
      	cid = channel_id[ii];
	// add volumes to correct channel contributions
      	ch_in_sf[ii] += (channel_area[ii]*precip[cid]+ q_sf_in[cid])*Dt; // volume of rainfall and flow from surface
	ch_in_sz[ii] += q_sz_in[cid]*Dt; // volume of flow to channel from saturated zone

	// mass balance contribution
      	mbv[2] += precip[cid]*channel_area[ii]*Dt; // rainfall into system 
      	mbv[3] += ( channel_area[ii]*precip[cid] + q_sf_in[cid] + q_sz_in[cid] ) * Dt; // volume lost to channel
      }
      // end of substep loop
    }
    
    // final state volumes for mass balance
    for(int ii=0; ii<nhillslope; ++ii){
      mbv[4] += area[ii]*(s_sf[ii] + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }
    mbv[5] = mbv[0] - mbv[1] + mbv[2] - mbv[3] - mbv[4];// mass balance error
    // copy mass balance to record
    for(int ii=0; ii < 6; ++ii){
      mass_balance(it,ii) = mbv[ii];
    }
    
    // convert channel inflow to rate
    for(int ii=0; ii < nchannel; ++ii){
      channel_inflow_sf(it,ii) = ch_in_sf[ii]/timestep;
      channel_inflow_sz(it,ii) = ch_in_sz[ii]/timestep;
    }
    
    // keep states if required
    if( keep_states[it] ){
      state_rec(it) = Rcpp::clone(states);
    }

    // check user interupt
    Rcpp::checkUserInterrupt(); 
  }
}


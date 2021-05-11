#include "Rcpp.h"
#include <boost/math/tools/roots.hpp>

//* ======================================================================================= *//
//* ======================================================================================= *//
//* ======================================================================================= *//

// Function for solving
// [[Rcpp::export]]
void dt_exp_explicit(Rcpp::DataFrame hillslope, // hillslope data frame
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

  
  // safety factor for computing computational timestep
  double alpha = 0.7;
  double min_step = 1;
  
  double Dt = timestep;
  
  
  
  // create vectors for storing lateral fluxes, precip and pet
  std::vector<double> q_sf_in(maxid+1,0.0), q_sz_in(maxid+1,0.0);
  std::vector<double> q_sf_out(maxid+1,0.0), q_sz_out(maxid+1,0.0);
  std::vector<double> precip(maxid+1,0.0), pet(maxid+1,0.0);
  
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
  
  
  
  // variable use with loop
  int cid; // current id
  double r_sf_rz,r_rz_uz,r_uz_sz; // fluxes between stores
  double tDt(0.0); // total time evaluated in step
  
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
    
    // set tDt to 0
    tDt = 0.0;
    
    // start loop of substeps
    while( tDt < timestep ){
  
      //Rcpp::Rcout << "Time step "<< it << " time passed " << tDt << " of " << timestep << std::endl;
											  
      // set flow passing records to 0
      std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
      std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);
      std::fill(q_sf_out.begin(), q_sf_out.end(), 0.0);
      std::fill(q_sz_out.begin(), q_sz_out.end(), 0.0);

      // work out outflow fluxes and timestep
      Dt = timestep - tDt;

      for(int ii=0; ii<nhillslope; ++ii){
	// current id used to reference longer vectors
    	cid = id[ii];
	// apply Dt estimate relating to surface
	double cr = c_sf[ii];
	Dt = std::min( Dt, alpha*Dx[ii]/cr );
	// Apply Dt estimate relating to subsurface
	cr = cosbeta_m[ii]*l_sz_max[ii]*std::exp(-cosbeta_m[ii]*s_sz[ii]);
	Dt = std::min( Dt, alpha*Dx[ii]/cr );
      }
      // Rcpp::Rcout << Dt << std::endl;
			     
      Dt = std::max(Dt,min_step);
      
      
      // set all counters to initial value
      link_cntr = 0; // counter for flow links
      link_from_id = flow_from[link_cntr];
      
      // loop HSUs
      for(int ii=0; ii<nhillslope; ++ii){
	
    	//Rcpp::Rcout << it << " " << nn << " " << ii << std::endl;
    	// current id used to reference longer vectors
    	cid = id[ii];
	
	// apply pet loss and precip input to mass balance
    	mbv[1] -= pet[cid]*(s_rz[ii]/s_rz_max[ii])*area[ii]*Dt;
    	mbv[2] += precip[cid]*area[ii]*Dt;
	
    	// compute first downward flux estimate from surface and root zone
    	// these are given as \check{r} in documentation	
    	r_sf_rz = std::min( r_sf_max[ii] , (s_sf[ii] + Dt*q_sf_in[cid]/area[ii])/Dt );
    	r_rz_uz = std::max( 0.0 ,
    			    (s_rz[ii] + Dt*(precip[cid] + r_sf_rz - pet[cid]*(s_rz[ii]/s_rz_max[ii])) - s_rz_max[ii])/Dt);

	if( s_sz[ii] > 0.0 ){
	  r_uz_sz = std::min( s_uz[ii]/(t_d[ii]*s_sz[ii]), (s_uz[ii] + Dt*r_rz_uz)/Dt );	  
	}else{
	  r_uz_sz = std::min( 1/t_d[ii], (s_uz[ii] + Dt*r_rz_uz)/Dt );
	}

	// if(cid==734){
	//   Rcpp::Rcout << "Inflow " << q_sf_in[cid] << " " << q_sz_in[cid] << std::endl;
	//   Rcpp::Rcout << "uzcalc " << s_uz[ii] << " " << s_sz[ii] << " " << t_d[ii]*s_sz[ii] << std::endl;
	//   Rcpp::Rcout << "Before " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
	// }
	
	// solve saturated zone
	q_sz_out[cid] = l_sz_max[ii]*std::exp(-cosbeta_m[ii]*s_sz[ii]);
	r_uz_sz = std::min( r_uz_sz ,
			    (s_sz[ii] - Dt*(q_sz_in[cid]-q_sz_out[cid])/area[ii])/Dt ); 
	s_sz[ii] = s_sz[ii] - Dt*r_uz_sz -Dt*(q_sz_in[cid]-q_sz_out[cid])/area[ii];
	
	// solve unsat
	r_rz_uz = std::min( r_rz_uz, (s_sz[ii] - s_uz[ii] +Dt*r_uz_sz)/Dt );
	s_uz[ii] += Dt*(r_rz_uz-r_uz_sz);
	//s_uz[ii] = std::max(0.0,s_uz[ii]);
	
	// solve root zone
	r_sf_rz = std::min( r_sf_rz,
			    (s_rz_max[ii] - (s_rz[ii] + Dt*(precip[cid] - r_rz_uz - pet[cid]*(s_rz[ii]/s_rz_max[ii]))))/Dt );
	s_rz[ii] += Dt*(precip[cid] + r_sf_rz - r_rz_uz - pet[cid]*(s_rz[ii]/s_rz_max[ii]));
	//s_rz[ii] = std::max(0.0,s_rz[ii]);
	// solve surface
	q_sf_out[cid] = std::min( width[ii]*s_sf[ii]*c_sf[ii],
				  (area[ii]*s_sf[ii] - Dt*area[ii]*r_sf_rz + Dt*q_sf_in[cid])/Dt );
	s_sf[ii] += Dt*(q_sf_in[cid]-q_sf_out[cid])/area[ii] - Dt*r_sf_rz;
	//s_sf[ii] = std::max(0.0,s_sf[ii]);
	
	// if(cid==734){
	//   Rcpp::Rcout << "After " << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
	//   Rcpp::Rcout << "Outflow " << q_sf_out[cid] << " " << q_sz_out[cid] << std::endl;
	// }
	
	//Rcpp::Rcout << link_cntr << " " << link_from_id << std::endl;
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
    	// end of hillslope loop
    	
      }
      
      // loop channels for the sub step
      for(int ii=0; ii < nchannel; ++ii){
      	cid = channel_id[ii];
      	double chn_in = (channel_area[ii]*precip[cid]+ q_sf_in[cid] + q_sz_in[cid])*Dt;
      	mbv[2] += precip[cid]*channel_area[ii]*Dt; 
      	ch_in[ii] += chn_in; // volume of flow to channel
      	mbv[3] -= chn_in; // volume lost from hillslope to channel
      }

      // add to time
      tDt += Dt;
      
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


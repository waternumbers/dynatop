#include "Rcpp.h"
#include <vector>
#include "hru.h"
#include "helpers.h"

// ///////////////////////////////////////
// Initialisation
// ///////////////////////////////////////
// [[Rcpp::export]]
void dt_init(Rcpp::List mdl, // hru data frame
	     double const vtol,
	     double const etol,
	     int const max_it
	    ){

  // dimensions and constants
  int nhru = mdl.size(); // number of HRUs

  // storage for inflow fluxes
  std::vector<double> q_sf_in(nhru,0.0);// vector to surface inflow fluxes
  std::vector<double> q_sz_in(nhru,0.0);// vector to saturated zone inflow fluxes
  
  // make HRUs
  std::vector<hru> hrus = makeHRUs(mdl);

  // start loop of hrus
  for(int ii=0; ii<nhru; ++ii){
    //Rcpp::Rcout << "hru " << ii << std::endl;
    Rcpp::List L = mdl[ii];
    Rcpp::NumericVector ivec = L["initialisation"];
    
    hrus[ii].init(q_sf_in,q_sz_in,ivec["s_rz_0"], ivec["r_uz_sz_0"], vtol, etol, max_it);
  }
  
  // copy back states
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List tmp = mdl[ii];
    Rcpp::NumericVector svec = tmp["states"];
    svec["s_sf"] = hrus[ii].s_sf;
    svec["s_rz"] = hrus[ii].s_rz;
    svec["s_uz"] = hrus[ii].s_uz;
    svec["s_sz"] = hrus[ii].s_sz;
  };

  //end of dt_init
}

// ////////////////////////////////////////
// Simulation
// ///////////////////////////////////////

// [[Rcpp::export]]
void dt_sim(Rcpp::List mdl, // list of HRUs
	    Rcpp::DataFrame out_dfn,
	    std::vector<bool> keep_states,
	    Rcpp::NumericMatrix obs_matrix, // external series
	    Rcpp::NumericMatrix mass_balance, // mass balance for each timestep
	    Rcpp::NumericMatrix out_matrix, // output series to populate
	    Rcpp::List state_rec,
	    double const timestep,
	    int const n_sub_step,
	    double const vtol,
	    double const etol,
	    int const max_it
	    ){
  
  // dimensions
  int nhru = mdl.size(); // number of HRUs

  // constant used in simulation
  double dbl_n_sub_step = (double)n_sub_step;
  double const Dt = timestep / dbl_n_sub_step;


  // create storage for input and output series
  std::vector<double> obs(obs_matrix.ncol(),0.0);// vector to store observed values
  std::vector<double> out(out_matrix.ncol(),0.0);// vector to store observed values
  std::vector<double> mbv(6,0.0);// vector to store mass balance calulations
  
  // sotrage for inflow fluxes
  std::vector<double> q_sf_in(nhru,0.0);// vector to surface inflow fluxes
  std::vector<double> q_sz_in(nhru,0.0);// vector to saturated zone inflow fluxes
  
  // make HRUs
  std::vector<hru> hrus = makeHRUs(mdl);
  
  // create output flux object
  outFlux out_flux(out_dfn["name_idx"], out_dfn["id_idx"], out_dfn["flux_int"]);
  
  // start loop of time steps
  for(int tt = 0; tt < obs_matrix.nrow(); ++tt) {
    // Rcpp::Rcout << "Time step " << tt << std::endl;

    // copy obs values
    for(unsigned int ii=0; ii < obs.size(); ++ii){
      obs[ii] = obs_matrix(tt,ii) / timestep ; // convert to rate
    }
    
    // clear vectors valid for all time step
    std::fill( mbv.begin(), mbv.end(), 0.0) ;
    std::fill( out.begin(), out.end(), 0.0) ;

    // update precip and pet in hrus
    for(int ii=0; ii<nhru; ++ii){
      hrus[ii].update_met(obs);
    }
    
    // compute the mass balance initial storage
    for(int ii=0; ii<nhru; ++ii){
      mbv[0] += (hrus[ii].s_sf + hrus[ii].s_rz + hrus[ii].s_uz - hrus[ii].s_rz) * hrus[ii].area; // initial state volume
      mbv[1] += hrus[ii].precip * hrus[ii].area; // precip volume
    }
    mbv[1] = mbv[1] * timestep;
    
    // Rcpp::Rcout << "mbv after initialised of time step inputs " << mbv[0] << std::endl;
    // Rcpp::Rcout << "summed Precip: " << std::accumulate(precip.begin(), precip.end(), 0.0) << std::endl;
    // Rcpp::Rcout << "summed pet: " << std::accumulate(pet.begin(), pet.end(), 0.0) << std::endl;
    
    // start loop of substeps
    for(int nn = 0; nn < n_sub_step; ++nn){

      std::fill( q_sf_in.begin(), q_sf_in.end(), 0.0) ;
      std::fill( q_sz_in.begin(), q_sz_in.end(), 0.0) ;
      //Rcpp::Rcout << "cleared flux" << std::endl;
  
      // start loop of hrus
      for(int ii=0; ii<nhru; ++ii){
  	//Rcpp::Rcout << "hru " << ii << " at timestep " << tt << std::endl;
  	hrus[ii].step(q_sf_in,q_sz_in,vtol,etol,max_it,Dt);

	// mass balance components
	mbv[2] += hrus[ii].aet * hrus[ii].area * Dt ; // actual evapotranspiration
	mbv[3] += Dt * (hrus[ii].q_sf + hrus[ii].q_sz - q_sf_in[ii] - q_sz_in[ii]) ; // net lateral flux
	
      }
      

      out_flux.apply( hrus, out, dbl_n_sub_step);
      
      // end loop of substeps
    }

    // finish off mass balance at end of step
    for(int ii=0; ii<nhru; ++ii){
      mbv[4] += (hrus[ii].s_sf + hrus[ii].s_rz + hrus[ii].s_uz - hrus[ii].s_rz) * hrus[ii].area; // final state volume
    }
    mbv[5] = mbv[0] + mbv[1] - mbv[2] - mbv[3] - mbv[4];
    for(unsigned int ii=0;  ii<6; ++ii){
      mass_balance(tt,ii) = mbv[ii];
    }
    //Rcpp::Rcout << "finished mass balance" << std::endl;
    
    // copy across output
    //Rcpp::Rcout << "Copying " << outt.size() << " outputs" << std::endl;
    for(unsigned int ii=0;  ii<out.size(); ++ii){
      out_matrix(tt,ii) = out[ii];
    }
    //Rcpp::Rcout << "copied output" << std::endl;
    
    // keep states if required
    //Rcpp::Rcout << "keep states" << std::endl;
    
    if( keep_states[tt] ){
      state_rec(tt) = makeStateList(hrus);
    }
    
    // check user interupt
    Rcpp::checkUserInterrupt(); 
  }

  
  // Rcpp::Rcout << "start copying states" << std::endl;
  // copy back states
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List tmp = mdl[ii];
    Rcpp::NumericVector svec = tmp["states"];
    svec["s_sf"] = hrus[ii].s_sf;
    svec["s_rz"] = hrus[ii].s_rz;
    svec["s_uz"] = hrus[ii].s_uz;
    svec["s_sz"] = hrus[ii].s_sz;
  };

  // end of dt_sim
}

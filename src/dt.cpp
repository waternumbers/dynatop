#include "Rcpp.h"
#include <vector>
#include <algorithm>
#include <execution>
#include <thread>
#include "hru.h"
#include "helpers.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

// ///////////////////////////////////////
// Initialisation
// ///////////////////////////////////////
// [[Rcpp::export]]
void dt_init(Rcpp::List mdl, // hru data frame
	     double const vtol,
	     double const etol,
	     int const max_it,
	     unsigned int const n_thread
	    ){
#if RCPP_PARALLEL_USE_TBB
  #include <tbb/global_control.h>
  unsigned int nt = std::min(n_thread, std::thread::hardware_concurrency()-1);
  Rcpp::Rcout << "Number of threads " << nt << std::endl;
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, nt);
  auto policy = std::execution::par;
#else
  auto policy = std::execution::seq;
#endif


  // set execution policy
  //auto policy = std::execution::seq;
  
  // dimensions and constants
  int nhru = mdl.size(); // number of HRUs

  // storage for inflow fluxes
  std::vector<double> q_sf_in(nhru,0.0);// vector to surface inflow inflow volumes
  std::vector<double> q_sz_in(nhru,0.0);// vector to saturated zone inflow volumes

  // set a time step
  double const Dt(0.0);
  
  // make HRUs
  // Rcpp::Rcout << "making HRUs" << std::endl;
  std::vector<hru> hrus = makeHRUs(mdl,q_sf_in,q_sz_in,vtol,etol,max_it,Dt);
  // Rcpp::Rcout << "hru size is " << nhru << " " << hrus.size() << std::endl;

  // make vector of break points between the bands
  std::vector<int> band_edge{0};
  for(int ii = 1; ii < nhru; ii++){
    if( hrus[ii].band != hrus[ii-1].band ){
      band_edge.push_back( ii );
    }
  }
  band_edge.push_back(hrus.size());
  
  // start loop of hrus
  for(int ii=band_edge.size()-1; ii >0; ii--){ // loop bands
    std::for_each(policy,hrus.begin() + band_edge[ii-1],
		  hrus.begin() + band_edge[ii],
		  []( hru &h ){ h.init(); }
		  );
    // TODO loop to add to downstream
    std::for_each(std::execution::seq,hrus.begin() + band_edge[ii-1],
		  hrus.begin() + band_edge[ii],
		  []( hru &h ){ h.lateral_redistribution(); }
		  );
  }
  // std::for_each( policy, hrus.rbegin(),hrus.rend(),
  // 		 []( hru &h ){ h.init(); } );
		 

  // copy back states
  for(int ii=0; ii<nhru; ++ii){
    // Rcpp::Rcout << "copy hru " << ii << std::endl;
    Rcpp::List tmp = mdl[ii];
    Rcpp::NumericVector svec = tmp["states"];
    svec["s_sf"] = hrus[ii].s_sf / hrus[ii].area;
    svec["s_rz"] = hrus[ii].s_rz / hrus[ii].area;
    svec["s_uz"] = hrus[ii].s_uz / hrus[ii].area;
    svec["s_sz"] = hrus[ii].s_sz / hrus[ii].area;
  };
  //  Rcpp::Rcout << "reacched end" << std::endl;
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
	    int const max_it,
	    unsigned int const n_thread
	    ){
  // Rcpp::Rcout << "Entered function" << std::endl;
#if RCPP_PARALLEL_USE_TBB
  #include <tbb/global_control.h>
  unsigned int nt = std::min(n_thread, std::thread::hardware_concurrency()-1);
  Rcpp::Rcout << "Number of threads " << nt << std::endl;
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, nt);
  //auto policy = std::execution::par;
#else
  unsigned int nt{1}
#endif

  bool isSeq = true;
  if( nt>1 ){ isSeq = false; }
  
  // set execution policy
  //auto policy = std::execution::seq;
  
  // dimensions
  int nhru = mdl.size(); // number of HRUs

  // constant used in simulation
  double const dbl_n_sub_step = (double)n_sub_step;
  double const Dt = timestep / dbl_n_sub_step;


  // create storage for input and output series
  std::vector<double> obs(obs_matrix.ncol(),0.0);// vector to store observed values
  std::vector<double> out(out_matrix.ncol(),0.0);// vector to store observed values
  std::vector<double> mbv(6,0.0);// vector to store mass balance calulations
  
  // storage for inflow fluxes
  std::vector<double> q_sf_in(nhru,0.0);// vector to surface inflow fluxes
  std::vector<double> q_sz_in(nhru,0.0);// vector to saturated zone inflow fluxes
  
  // make HRUs
  std::vector<hru> hrus = makeHRUs(mdl,q_sf_in,q_sz_in,vtol,etol,max_it,Dt);
  //Rcpp::Rcout << "Made HRUs" << std::endl;


  
  // create output flux object
  outFlux out_flux(out_dfn["name_idx"], out_dfn["id_idx"], out_dfn["flux_int"], out_dfn["scale"], dbl_n_sub_step);
  //Rcpp::Rcout << "Made outFlux" << std::endl;

  // make vector of break points between the bands
  std::vector<int> band_edge{0};
  for(int ii = 1; ii < nhru; ii++){
    if( hrus[ii].band != hrus[ii-1].band ){
      band_edge.push_back( ii );
    }
  }
  band_edge.push_back(hrus.size());
  
  // start loop of time steps
  for(int tt = 0; tt < obs_matrix.nrow(); ++tt) {
    //Rcpp::Rcout << "Time step " << tt << std::endl;

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

    // Rcpp::Rcout << hrus[0].s_sf << " " << hrus[0].s_rz << " " << hrus[0].s_uz << " " << hrus[0].s_sz << std::endl;
    
    // compute the mass balance initial storage
    for(int ii=0; ii<nhru; ++ii){
      //      if( hrus[ii].area > 0.0){
	mbv[0] += (hrus[ii].s_sf + hrus[ii].s_rz + hrus[ii].s_uz - hrus[ii].s_sz); // initial state volume
	mbv[1] += hrus[ii].precip; // precip volume
	//}
    }
    mbv[1] = mbv[1] * timestep;
    
    //Rcpp::Rcout << "mbv after initialised of time step inputs " << mbv[0] << std::endl;
    //Rcpp::Rcout << "summed Precip: " << std::accumulate(precip.begin(), precip.end(), 0.0) << std::endl;
    //Rcpp::Rcout << "summed pet: " << std::accumulate(pet.begin(), pet.end(), 0.0) << std::endl;

    
 
    
    // start loop of substeps
    for(int nn = 0; nn < n_sub_step; ++nn){


      std::fill( q_sf_in.begin(), q_sf_in.end(), 0.0) ;
      std::fill( q_sz_in.begin(), q_sz_in.end(), 0.0) ;
      //Rcpp::Rcout << "cleared flux" << std::endl;
  
      // start loop of hrus
      if( isSeq ){
	for(int ii=band_edge.size()-1; ii >0; ii--){ // loop bands
	  std::for_each(std::execution::seq,hrus.begin() + band_edge[ii-1],
			hrus.begin() + band_edge[ii],
			[]( hru &h ){ h.step(); } );
	  // loop to add to downstream
	  std::for_each(std::execution::seq,hrus.begin() + band_edge[ii-1],
			hrus.begin() + band_edge[ii],
			[]( hru &h ){ h.lateral_redistribution(); }
			);
	}
	// this is quicker but not quite binary compatable
	// std::for_each(std::execution::seq,hrus.rbegin(),hrus.rend(),
	// 	      []( hru &h ){ h.step(); h.lateral_redistribution(); }
	// 	      );
      }else{
	for(int ii=band_edge.size()-1; ii >0; ii--){ // loop bands
	  std::for_each(std::execution::par,hrus.begin() + band_edge[ii-1],
			hrus.begin() + band_edge[ii],
			[]( hru &h ){ h.step(); } );
	  // loop to add to downstream
	  std::for_each(std::execution::seq,hrus.begin() + band_edge[ii-1],
			hrus.begin() + band_edge[ii],
			[]( hru &h ){ h.lateral_redistribution(); }
			);
	}
      }
      
      for(int ii= nhru-1; ii >= 0; --ii){
	// mass balance components
	//	if( hrus[ii].area > 0.0){
	  mbv[2] += hrus[ii].aet * Dt ; // actual evapotranspiration
	  mbv[3] += Dt * (hrus[ii].q_sf + hrus[ii].q_sz - hrus[ii].q_sf_in - hrus[ii].q_sz_in) ; // net lateral flux
	  //	}
      }
      

      out_flux.apply( hrus, out) ; //, dbl_n_sub_step);
      
      // end loop of substeps
    }
    //Rcpp::Rcout << hrus[0].s_sf << " " << hrus[0].s_rz << " " << hrus[0].s_uz << " " << hrus[0].s_sz << std::endl;
    //Rcpp::Rcout << hrus[0].r_sf_rz << " " << hrus[0].r_rz_uz << " " << hrus[0].r_uz_sz << std::endl;
    //Rcpp::Rcout << hrus[0].q_sf_in << " " << hrus[0].q_sf << " " << hrus[0].q_sz_in << " " << hrus[0].q_sz << std::endl;
    
    // finish off mass balance at end of step
    for(int ii=0; ii<nhru; ++ii){
      //      if( hrus[ii].area > 0.0){
      mbv[4] += (hrus[ii].s_sf + hrus[ii].s_rz + hrus[ii].s_uz - hrus[ii].s_sz); // final state volume
	//}
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
      //state_rec(tt) = makeStateList(hrus);
      state_rec(tt) = makeStateDataFrame(hrus);
    }
    
    // check user interupt
    Rcpp::checkUserInterrupt(); 
  }

  
  // Rcpp::Rcout << "start copying states" << std::endl;
  // copy back states
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List tmp = mdl[ii];
    Rcpp::NumericVector svec = tmp["states"];
    svec["s_sf"] = hrus[ii].s_sf / hrus[ii].area;
    svec["s_rz"] = hrus[ii].s_rz / hrus[ii].area;
    svec["s_uz"] = hrus[ii].s_uz / hrus[ii].area;
    svec["s_sz"] = hrus[ii].s_sz / hrus[ii].area;
  };

  // end of dt_sim
}

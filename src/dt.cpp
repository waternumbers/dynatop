#include "Rcpp.h"
#include <vector>
#include "hru.h"
#include "lnk.h"
#include "helpers.h"

// ///////////////////////////////////////
// Initialisation
// ///////////////////////////////////////
// [[Rcpp::export]]
void dt_init(Rcpp::DataFrame mdl, // hru data frame
	     double const vtol,
	     double const qtol,
	     int const max_it
	    ){
  
  // dimensions and constants
  int nhru = mdl.nrows(); // number of HRUs
  const double Dt = -99.9;
  std::vector<double> s_rz_0 = mdl["s_rz0"];
  std::vector<double> r_uz_sz_0 = mdl["r_uz_sz0"];
  
  // create storage for all the states and fluxes that might be output
  std::vector<double> s_sf = mdl["s_sf"];
  std::vector<double> s_sp = mdl["s_sp"];
  std::vector<double> s_rz = mdl["s_rz"];
  std::vector<double> s_uz = mdl["s_uz"];
  std::vector<double> s_sz = mdl["s_sz"];

  std::vector<double> mbv(6,0.0);
  std::vector<double> precip(nhru,0.0), pet(nhru,0.0);
  std::vector<double> q_sf(nhru,0.0), q_sz(nhru,0.0);
  std::vector<double> q_sf_in(nhru,0.0), q_sz_in(nhru,0.0);
  std::vector<double> r_sf_sp(nhru,0.0), r_sf_rz(nhru,0.0), r_sp_rz(nhru,0.0), r_rz_uz(nhru,0.0), r_uz_sz(nhru,0.0);
  std::vector<double> e_a(nhru,0.0);

  Rcpp::Rcout << "make HRUS" << std::endl;
  // create hrus - in doing this the data is copied into C++ structures within the HRU class, except those declared above
  std::vector<hru> hrus  = makeHRUs(mdl,
				    s_sf, s_sp, s_rz, s_uz, s_sz,
				    precip, pet,
				    q_sf, q_sz, q_sf_in, q_sz_in,
				    r_sf_sp, r_sf_rz, r_sp_rz, r_rz_uz, r_uz_sz, e_a,
				    Dt, vtol,qtol,max_it);

  // Rcpp::Rcout << "make links" << std::endl;
  // create the post hru evaluation links
  std::vector< slnk > sf_link, sz_link;
  Rcpp::List sf_fd = mdl["sf_flow_direction_cpp"];
  Rcpp::List sz_fd = mdl["sz_flow_direction_cpp"];
  Rcpp::List tmp;
  Rcpp::IntegerVector tmp_int;
  
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List tmp = sf_fd[ii];
    sf_link.push_back( slnk(  tmp["id_idx"], tmp["frc"] ) );
    tmp = sz_fd[ii];
    sz_link.push_back( slnk(  tmp["id_idx"], tmp["frc"] ) );
  }

    // Rcpp::Rcout << "loop HRUS" << std::endl;  
  // start loop of hrus
  for(int ii=0; ii<nhru; ++ii){
    //Rcpp::Rcout << "hru " << ii << std::endl;
    
    hrus[ii].init(s_rz_0[ii], r_uz_sz_0[ii]);
	
    // run post links
    sf_link[ii].apply(q_sf[ii], q_sf_in);
    sz_link[ii].apply(q_sz[ii], q_sz_in);

    //Rcpp::Rcout << "redistributed hru " << ii << std::endl;
    // end loop of hrus
  }

  // Rcpp::Rcout << "start copying states" << std::endl;
  // copy back states
  //for(int ii=0; ii<nhru; ++ii){
  mdl["s_sz"] = s_sz;
  mdl["s_uz"] = s_uz;
  mdl["s_rz"] = s_rz;
  mdl["s_sf"] = s_sf;
  mdl["q_sf"] = q_sf;
  mdl["q_sz"] = q_sz;
    //}
}

// ////////////////////////////////////////
// Simulation
// ///////////////////////////////////////

// [[Rcpp::export]]
void dt_sim(Rcpp::DataFrame mdl, // hru data frame
	    Rcpp::NumericMatrix obs, // external series
	    Rcpp::NumericMatrix mass_balance, // mass balance for each timestep
	    Rcpp::NumericMatrix out, // output series to populate
	    Rcpp::DataFrame out_dfn,
	    std::vector<bool> keep_states,
	    Rcpp::List state_rec,
	    double const timestep,
	    int const n_sub_step,
	    double const vtol,
	    double const qtol,
	    int const max_it
	    ){
  
  // dimensions and constants
  int nhru = mdl.nrows(); // number of HRUs
  //double one = 1.0;
  //double minus_one = -1.0;
  double dbl_n_sub_step = (double)n_sub_step;
  double const Dt = timestep / dbl_n_sub_step;

  // scalings used in links
  // Rcpp::Rcout << "initialisation part 1" << std::endl;
 
  // create storage for all the states and fluxes that migh be output
  std::vector<int> id = mdl["id"];
  std::vector<double> s_sf = mdl["s_sf"];
  std::vector<double> s_rz = mdl["s_rz"];
  std::vector<double> s_uz = mdl["s_uz"];
  std::vector<double> s_sz = mdl["s_sz"];
  std::vector<double> q_sz = mdl["q_sz"];
  std::vector<double> q_sf = mdl["q_sf"];
  std::vector<double> precip(nhru,0.0), pet(nhru,0.0);
  std::vector<double> q_sf_in(nhru,0.0), q_sz_in(nhru,0.0); // inflows at end of step
  std::vector<double> q_sf_prev(nhru,0.0), q_sz_prev(nhru,0.0); // outflows at start to step
  std::vector<double> r_sf_rz(nhru,0.0), r_rz_uz(nhru,0.0), r_uz_sz(nhru,0.0);
  std::vector<double> a_e(nhru,0.0);
  std::vector<double> mbv(6,0.0); // mass balance vector
  std::vector<double> obst(obs.ncol(),0.0);// vector to store observed values
  std::vector<double> outt(out.ncol(),0.0);// vector to store observed values
  std::vector< alnk > precip_link, pet_link; // links to evaluate at start of step
  std::vector< slnk > sf_link, sz_link; // flow redistribution

  // Rcpp::Rcout << "initialisation part 2" << std::endl;
  //std::vector< std::vector<lnk> > post_lnks(nhru); // links to run after evaluating each HRU in substep
  //std::vector<lnk> out_lnks; // links to run after evaluating every HRU in substep
  //std::vector<lnk> end_lnks; // links to run at end of step

  // create hrus - in doing this the data is copied into C++ structures within the HRU class, except those declared above
  // Rcpp::Rcout << "start HRU creation" << std::endl;
  std::vector<hru> hrus  = makeHRUs(mdl,
				    s_sf, s_rz, s_uz, s_sz, q_sz, q_sf,
				    precip, pet,
				    q_sf_in, q_sz_in, q_sf_in_prev, q_sz_in_prev,
				    r_sf_rz, r_rz_uz, r_uz_sz, a_e,
				    Dt, vtol,qtol,max_it);

  // Rcpp::Rcout << "made HRUS" << std::endl;
  
  // create links
  Rcpp::List sf_fd = mdl["sf_flow_direction_cpp"];
  Rcpp::List sz_fd = mdl["sz_flow_direction_cpp"];
  Rcpp::List prcp = mdl["precip_cpp"];
  Rcpp::List pt = mdl["pet_cpp"];
  
  for(int ii=0; ii<nhru; ++ii){
    Rcpp::List tmp = sf_fd[ii];
    sf_link.push_back( slnk(  tmp["id_idx"], tmp["frc"] ) );
    tmp = sz_fd[ii];
    sz_link.push_back( slnk(  tmp["id_idx"], tmp["frc"] ) );
    tmp = prcp[ii];
    precip_link.push_back( alnk( tmp["col_idx"], tmp["area"] ) );
    tmp = pt[ii];
    pet_link.push_back( alnk( tmp["col_idx"], tmp["area"] ) );
  }

  // Rcpp::Rcout << "made links" << std::endl;

  // create output flux object
  outFlux out_flux(out_dfn["name_idx"], out_dfn["id_idx"], out_dfn["flux_int"]);

  // initialise the inflow at the starting time step
  for(int ii=0; ii<nhru; ++ii){
    // run post links
    sf_link[ii].apply(q_sf[ii], q_sf_in);
    sz_link[ii].apply(q_sz[ii], q_sz_in);
  }

  
  // start loop of time steps
  for(int tt = 0; tt < obs.nrow(); ++tt) {
    // Rcpp::Rcout << "Time step " << tt << std::endl;

    // copy obs values
    for(unsigned int ii=0; ii < obst.size(); ++ii){
      obst[ii] = obs(tt,ii) / timestep ; // convert to rate
    }
    //Rcpp::Rcout << "Initialised observed series " << std::endl;
    
    // clear vectors valid for all time step
    std::fill( mbv.begin(), mbv.end(), 0.0) ;
    std::fill( precip.begin(), precip.end(), 0.0) ;
    std::fill( pet.begin(), pet.end(), 0.0) ;
    std::fill( outt.begin(), outt.end(), 0.0) ;
    //Rcpp::Rcout << "Initialised fluxes " << std::endl;
    
    // initialise the precip and pet for each HRU - these are in m^3/s
    for(int ii=0; ii < nhru; ++ii){
      precip_link[ii].apply(obst, precip[ii]);
      pet_link[ii].apply(obst, pet[ii]);
    }

    //Rcpp::Rcout << "After applying precip links " << precip[0] << std::endl;
    
    // compute the mass balance initial storage
    mbv[0] = std::accumulate(s_sf.begin(), s_sf.end(), 0.0);
    mbv[0] += std::accumulate(s_rz.begin(), s_rz.end(), 0.0);
    mbv[0] += std::accumulate(s_uz.begin(), s_uz.end(), 0.0);
    mbv[0] -= std::accumulate(s_sz.begin(), s_sz.end(), 0.0);

    // compute the mass balance precip input
    mbv[1] = std::accumulate(precip.begin(), precip.end(), 0.0);
    mbv[1] = mbv[1] * timestep;

    // Rcpp::Rcout << "mbv after initialised of time step inputs " << mbv[0] << std::endl;
    // Rcpp::Rcout << "summed Precip: " << std::accumulate(precip.begin(), precip.end(), 0.0) << std::endl;
    // Rcpp::Rcout << "summed pet: " << std::accumulate(pet.begin(), pet.end(), 0.0) << std::endl;
    
    // start loop of substeps
    for(int nn = 0; nn < n_sub_step; ++nn){

      //Rcpp::Rcout << "Substep " << nn << std::endl;

      
      // Move inflow record back to past value and set current flow passing records to 0
      // for(int ii=0; ii < nhru; ++ii){
      // 	q_sf_in_prev[ii] = q_sf_in[ii];
      // 	q_sz_in_prev[ii] = q_sz_in[ii];
      // }
      q_sf_in_prev = q_sz_in;
      q_sz_in_prev = q_sz_in;
      std::fill( q_sf_in.begin(), q_sf_in.end(), 0.0) ;
      std::fill( q_sz_in.begin(), q_sz_in.end(), 0.0) ;
      std::fill( a_e.begin(), a_e.end(), 0.0) ;
      //Rcpp::Rcout << "cleared flux" << std::endl;

      // Do mass balance part with the starting fluxes
      mbv[3] += std::accumulate(q_sf.begin(), q_sf.end(), 0.0) /2.0;
      mbv[3] += std::accumulate(q_sz.begin(), q_sz.end(), 0.0) /2.0;
      mbv[3] -= std::accumulate(q_sf_in_prev.begin(), q_sf_in_prev.end(), 0.0) /2.0;
      mbv[3] -= std::accumulate(q_sz_in_prev.begin(), q_sz_in_prev.end(), 0.0) /2.0;
      
      // start loop of hrus
      for(int ii=0; ii<nhru; ++ii){
  	//Rcpp::Rcout << "hru " << ii << " at timestep " << tt << std::endl;

  	hrus[ii].step(); //Dt,vtol,qtol,max_it);

	// run post links
	sf_link[ii].apply(q_sf[ii], q_sf_in);
	sz_link[ii].apply(q_sz[ii], q_sz_in);

  	// end loop of hrus
      }

      // do mass balance calc
      mbv[2] += std::accumulate(a_e.begin(), a_e.end(), 0.0);
      mbv[3] += std::accumulate(q_sf.begin(), q_sf.end(), 0.0) /2.0;
      mbv[3] += std::accumulate(q_sz.begin(), q_sz.end(), 0.0) /2.0;
      mbv[3] -= std::accumulate(q_sf_in.begin(), q_sf_in.end(), 0.0) /2.0;
      mbv[3] -= std::accumulate(q_sz_in.begin(), q_sz_in.end(), 0.0) /2.0;
      
      // Rcpp::Rcout << "q_sf " << std::accumulate(q_sf.begin(), q_sf.end(), 0.0) << std::endl;
      // Rcpp::Rcout << "q_sz " << std::accumulate(q_sz.begin(), q_sz.end(), 0.0)<< std::endl;
      // Rcpp::Rcout << "q_sf_in " << std::accumulate(q_sf_in.begin(), q_sf_in.end(), 0.0) << std::endl;
      // Rcpp::Rcout << "q_sf_in " << std::accumulate(q_sz_in.begin(), q_sz_in.end(), 0.0) << std::endl;

      out_flux.apply( outt, dbl_n_sub_step, precip, pet, a_e, q_sf, q_sf_in, q_sz, q_sz_in,
		      s_sf, s_rz, s_uz, s_sz, r_sf_rz, r_rz_uz, r_uz_sz);
      
      // end loop of substeps
    }

    // finish off mass balance at end of step
    mbv[2] = mbv[2] * Dt;
    mbv[3] = mbv[3] * Dt;
    mbv[4] = std::accumulate(s_sf.begin(), s_sf.end(), 0.0);
    mbv[4] += std::accumulate(s_rz.begin(), s_rz.end(), 0.0);
    mbv[4] += std::accumulate(s_uz.begin(), s_uz.end(), 0.0);
    mbv[4] -= std::accumulate(s_sz.begin(), s_sz.end(), 0.0);
    mbv[5] = mbv[0] + mbv[1] - mbv[2] - mbv[3] - mbv[4];
    for(unsigned int ii=0;  ii<6; ++ii){
      //Rcpp::Rcout << "mbv: " << mbv[ii] << std::endl;
      mass_balance(tt,ii) = mbv[ii];
      //Rcpp::Rcout << "mass_balance: " << mass_balance(tt,ii) << std::endl;
    }
    //Rcpp::Rcout << "finished mass balance" << std::endl;
    
    // copy across output
    //Rcpp::Rcout << "Copying " << outt.size() << " outputs" << std::endl;
    for(unsigned int ii=0;  ii<outt.size(); ++ii){
      out(tt,ii) = outt[ii];
    }
    //Rcpp::Rcout << "copied output" << std::endl;
    
    // keep states if required
    //Rcpp::Rcout << "keep states" << std::endl;
    
    if( keep_states[tt] ){
      state_rec(tt) = Rcpp::DataFrame::create( Rcpp::Named("id") = id,         // simple assign
  					 Rcpp::Named("s_sf") = s_sf,
  					 Rcpp::Named("s_rz") = s_rz,
  					 Rcpp::Named("s_uz") = s_uz,
  					 Rcpp::Named("s_sz") = s_sz,
  					 Rcpp::Named("q_sf") = q_sf,
  					 Rcpp::Named("q_sz") = q_sz);
      //Rcpp::Rcout << "copied states" << std::endl;
    }
    
    // check user interupt
    Rcpp::checkUserInterrupt(); 
  }

  // Rcpp::Rcout << "start copying states" << std::endl;
  // copy back states
  //for(int ii=0; ii<nhru; ++ii){
  mdl["s_sz"] = s_sz;
  mdl["s_uz"] = s_uz;
  mdl["s_rz"] = s_rz;
  mdl["s_sf"] = s_sf;
  mdl["q_sf"] = q_sf;
  mdl["q_sz"] = q_sz;
    //}
}

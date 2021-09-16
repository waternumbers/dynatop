#include "Rcpp.h"
#include "hillslope_hru.h"


// [[Rcpp::export]]
void dt_init(Rcpp::DataFrame hillslope, // hillslope data frame
	     Rcpp::DataFrame channel, // channel data frame
	     Rcpp::DataFrame flow_direction // flow directions data frame
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
  Rcpp::NumericVector m = hillslope["m"]; // transmissivity decay parameter
  Rcpp::NumericVector D = hillslope["D"]; // transmissivity decay parameter
  Rcpp::NumericVector m_2 = hillslope["m_2"]; // second transmissivity decay parameter
  Rcpp::NumericVector omega = hillslope["omega"]; // second transmissivity decay parameter
  Rcpp::IntegerVector opt = hillslope["opt"]; // type of saturated zone
  
  // parts of hillslope needed for initialisation
  Rcpp::NumericVector r_uz_sz_0 = hillslope["r_uz_sz0"];
  Rcpp::NumericVector s_rz_0 = hillslope["s_rz0"];

  // seperate out channel to vector
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
   
  // work out computational timestep - explicit casting of n_sub_step to double
  //double Dt = timestep / (double)n_sub_step;

  // create vectors for storing lateral fluxes, precip and pet
  std::vector<double> q_sf_in(maxid+1,0.0), q_sz_in(maxid+1,0.0);
  std::vector<double> q_sf_out(maxid+1,0.0), q_sz_out(maxid+1,0.0);
  std::vector<double> e_a(maxid+1,0.0);

  // variable use with loops
  int cid; // current id
 
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

  // set redistribution counters to initial value
  int link_cntr = 0; // counter for flow links
  int link_from_id = flow_from[link_cntr];
  
  // loop HSUs
  for(int ii=0; ii<nhillslope; ++ii){
    
    // current id used to reference longer vectors
    cid = id[ii];
    
    // evolve hillslope
    hs_hru[ii].init(s_rz_0[ii], r_uz_sz_0[ii]);

    // resdistribute flows
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
  }
}

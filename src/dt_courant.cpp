#include "Rcpp.h"
#include "hillslope_hru.h"


// [[Rcpp::export]]
void dt_courant(Rcpp::DataFrame hillslope, // hillslope data frame
		Rcpp::NumericMatrix courant, // courant numbers
		double timestep, // time step
		int n_sub_step // number of sub steps
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
  Rcpp::NumericVector D = hillslope["D"]; // maximum storage depth
  Rcpp::NumericVector m_2 = hillslope["m_2"]; // second transmissivity decay parameter
  Rcpp::NumericVector omega = hillslope["omega"]; // second transmissivity decay parameter
  Rcpp::IntegerVector opt = hillslope["opt"]; // type of saturated zone
  
  // work out some dimensions
  int nhillslope = id.size();
   
  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;

  // create dummy variable for passing as storing lateral fluxes, precip and pet
  double dummy_double(0.0);
    
  // variable for output
  std::pair<double,double> cr(0.0,0.0);
  
  // loop to create hillslope class objects
  //std::vector<hillslope> hs_hru;
  for(int ii=0; ii<nhillslope; ++ii){
    hillslope_hru hs = hillslope_hru(id[ii],
				     s_sf[ii], s_rz[ii], s_uz[ii], s_sz[ii],
				     s_bar[ii],   area[ii],   width[ii],
				     dummy_double, dummy_double, // surface zone lateral fluxes
				     dummy_double, dummy_double, // saturated zone lateral fluxes
				     dummy_double, // actual evapotranspiration as a rate [m/s]
				     r_sf_max[ii],   c_sf[ii], // surface store parameters
				     s_rz_max[ii], // root zone store parameters
				     t_d[ii], // unsaturated zone parameters
				     ln_t0[ii], c_sz[ii], m[ii], D[ii], m_2[ii], omega[ii],// saturated zone parameters
				     opt[ii]   ); //type of saturated zone
    cr = hs.courant(Dt);
    courant(ii,0) = cr.first;
    courant(ii,1) = cr.second;
  }
}

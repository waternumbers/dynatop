// #include "gperftools/profiler.h"
#include "Rcpp.h"
#include "dsu.h"

using namespace Rcpp;
// recall: Rcpp deep copies inputs unless they are the correct type!

// [[Rcpp::export]]
void single_hsu_cpp(NumericMatrix state_rec,
		    NumericMatrix ext_rec,
		    NumericMatrix flux_rec,
		    NumericVector properties,
		    double timestep,
		    int n_sub_step){
  // Rcout << "Running init" << std::endl;

  state_rec(4,3) = -999.0;
  
  // Deep copy the first rows of the state and ext to be stores
  //NumericMatrix::Row states = state_rec(0,_);
  NumericVector states = state_rec(0,_);
  NumericVector ext = ext_rec(0,_);
  NumericVector flux = flux_rec(0,_);
 
  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;
  Rcout<< Dt <<std::endl;

  std::vector<flink> tmp;
  
  // initialise the hsu
  hsu h(states[0],states[1],states[2],states[3],
  	ext[0],ext[1],ext[2],ext[3],
  	properties[0],properties[1],properties[2],
  	properties[3],properties[4],
  	properties[5],
  	properties[6],
  	properties[7],properties[8],
  	Dt, tmp
  	);

  // loop data timesteps
  for(int it = 1; it < ext_rec.nrow(); ++it) {
    Rcout << it << std::endl;
    ext = ext_rec(it,_);
    for(int ii = 0; ii < n_sub_step; ++ii){
      h.step();
    }
    state_rec(it,_) = states;
    flux = h.get_flux();
    flux_rec(it,_) = flux;
  }
  Rcout<< states[0] <<std::endl;

}

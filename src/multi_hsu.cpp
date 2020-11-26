// #include "gperftools/profiler.h"
#include "Rcpp.h"
#include "dsu.h"

using namespace Rcpp;
// recall: Rcpp deep copies inputs unless they are the correct type!

// [[Rcpp::export]]
void multi_hsu_cpp(IntegerVector id,
		   NumericMatrix states,
		   NumericMatrix properties,
		   List flow_dir,
		   IntegerMatrix ext_idx,
		   IntegerVector channel_id,
		   NumericVector channel_area,
		   IntegerVector channel_ext_idx,
		   NumericMatrix ext_rec,
		   NumericMatrix channel_inflow,
		   LogicalVector keep_states,
		   List state_rec,
		   double timestep,
		   int n_sub_step){
  // Rcout << "Running init" << std::endl;

  // Deep copy the first rows of ext to be the references store
  NumericVector ext = ext_rec(0,_);
 
  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;
  //Rcout<< Dt <<std::endl;

  // create vectors for storing lateral fluxes - make longer so can keep R id
  int maxid = max(id);
  if(max(channel_id)>maxid){ maxid = max(channel_id); }
  maxid +=1;
  NumericVector q_sf_in(maxid);
  NumericVector q_sz_in(maxid);
  
  // initialise the hsu
  std::vector<hsu> vhsu;

  int ip=0, iep=0, idi=0;
  
  for(int ii=0; ii<states.nrow(); ++ii){
    ip = ext_idx(ii,0);
    iep = ext_idx(ii,1);
    idi = id[ii];
    
    hsu h(states(ii,0),states(ii,1),states(ii,2),states(ii,3),
	  q_sf_in[idi], q_sz_in[idi],ext[ip],ext[iep],
	  properties(ii,0),properties(ii,1),properties(ii,2),
	  properties(ii,3),properties(ii,4),
	  properties(ii,5),
	  properties(ii,6),
	  properties(ii,7),properties(ii,8),
	  Dt
	  );
    vhsu.push_back(h);
  }

  // loop data timesteps
  for(int it = 0; it < ext_rec.nrow(); ++it) {
    Rcout << it << std::endl;
    ext = ext_rec(it,_);
    for(int nn = 0; nn < n_sub_step; ++nn){
      q_sf_in.fill(0.0);
      q_sz_in.fill(0.0);
      for(int ii=0; ii<states.nrow(); ++ii){
	vhsu[ii].step();
	List fd = flow_dir(ii);
	NumericVector fdf = fd("frc");
	IntegerVector fdi = fd("idx");
	std::vector<double> q = vhsu[ii].get_q();
	for(int jj=0; jj < fdi.length(); ++jj){
	  int fdi_ = fdi[jj];
	  double fdf_ = fdf[jj];
	  q_sf_in[fdi_] += q[0]*fdf_;
	  q_sz_in[fdi_] += q[1]*fdf_;
	}
      }
    }
    // populate channel inflow
    for(int ii=0; ii < channel_id.length(); ++ii){
      int id_ = channel_id[ii];
      int edx_ = channel_ext_idx[ii];
      channel_inflow(it,ii) = channel_area[ii]*ext[edx_] + q_sf_in[id_] + q_sz_in[id_];
    }

    // keep states if required
    if( keep_states(it) ){
      state_rec(it) = clone(states);
    }
  }
  
  //Rcout << states(0,0) << std::endl;
}

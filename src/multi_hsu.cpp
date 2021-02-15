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
		   int n_sub_step,
		   bool approx_soln){
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
  std::vector<double> q_sf_in(maxid,0.0), q_sz_in(maxid,0.0);
  std::vector<double> q_sf_out(maxid,0.0), q_sz_out(maxid,0.0);
  //std::vector<double> q_sz_in(maxid);

  // work out number of hillslopes
  int nhillslope = states.nrow();
  
  // convert all flow directions to links - do first to avoid pushback issues?
  List fd;
  NumericVector fdfrc;
  IntegerVector fdidx;
  int idi;
  std::vector< std::vector<flink> > paul;
  std::vector< double > link_from_hsu, link_to_hsu, link_frc;
  //std::vector< std::vector<flowlink> > links(nhillslope);
  for(int ii=0; ii<nhillslope; ++ii){
    idi = id[ii];
    fd = flow_dir(ii);
    fdidx = fd("idx");
    fdfrc = fd("frc");
    std::vector<flink> tmp;
    //std::vector<flowlink> fl(nfl);
    for(int i=0; i<fdidx.length(); ++i){
      int j = fdidx[i];
      flink lnksf(fdfrc[i],q_sf_out[idi],q_sf_in[j]);
      tmp.push_back(lnksf);
      flink lnksz(fdfrc[i],q_sz_out[idi],q_sz_in[j]);
      tmp.push_back(lnksz);

      link_from_hsu.push_back(idi);
      link_to_hsu.push_back(j);
      link_frc.push_back(fdfrc[i]);

      //fl[i].frc = &fdf[i];
      //fl[i].q_sf_in = &q_sf_in[j];
      //fl[i].q_sz_in = &q_sz_in[j];
    }
    paul.push_back(tmp); //[ii] = tmp;
    //links[ii] = fl;
    //paul.push_back(tmp);
  }
  
  
  // initialise the hsu
  std::vector<hsu> vhsu;
  //std::vector< std::tuple< std::vector<int>, std::vector<double> > > f_dir;
  
  int ip=0, iep=0;
  
  for(int ii=0; ii<nhillslope; ++ii){
    //Rcout << "making hsus " << ii << std::endl;
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
	  Dt, paul[ii], q_sf_out[idi], q_sz_out[idi]
	  );
    vhsu.push_back(h);
  }

  int lnk_cnt;
  int nlink = link_from_hsu.size();
  int j,k;
  // loop data timesteps
  for(int it = 0; it < ext_rec.nrow(); ++it) {
    //Rcout << it << std::endl;
    ext = ext_rec(it,_) / timestep;
    //Rcout << "ext is " << ext(0) << " " << ext(1) << std::endl;
    for(int nn = 0; nn < n_sub_step; ++nn){
      std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
      std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);
      //q_sf_in.fill(0.0);
      //q_sz_in.fill(0.0);
      lnk_cnt = 0;
      for(int ii=0; ii<nhillslope; ++ii){
	if(approx_soln){
	  vhsu[ii].astep();
	}else{
	  vhsu[ii].step();
	}

	// std::vector<flink>& lnks = paul[ii];
	// for(uint k=0; k<lnks.size(); ++k){
	//  lnks[k].eval();
	// }

	// idi = id[ii];
	// double q1 = q_sf_out[idi];
	// double q2 = q_sz_out[idi];
	// while(link_from_hsu[lnk_cnt]==idi & lnk_cnt < nlink){
	//   //j = link_from_hsu[lnk_cnt];
	//   k = link_to_hsu[lnk_cnt];
	//   double& f = link_frc[lnk_cnt];
	//   q_sf_in[k] += q1*f;//_sf_out[j]*f;
	//   q_sz_in[k] += q2*f;//q_sz_out[j]*f;
	//   ++lnk_cnt;
	// }
	  
      }

	
    }
    // populate channel inflow
    for(int ii=0; ii < channel_id.length(); ++ii){
      int id_ = channel_id[ii];
      int edx_ = channel_ext_idx[ii];
      channel_inflow(it,ii) = channel_area[ii]*ext[edx_] + q_sf_in[id_] + q_sz_in[id_];
      //channel_inflow(it,ii) = q_sz_in[id_];
    }

    // keep states if required
    if( keep_states(it) ){
      state_rec(it) = clone(states);
    }
  }

  //Rcout << states(0,0) << std::endl;
}


// [[Rcpp::export]]
void multi_hsu_cpp_init(IntegerVector id,
			NumericMatrix states,
			NumericMatrix properties,
			List flow_dir,
			NumericVector s_rz_0,
			NumericVector r_uz_sz_0
			)
{

  // Create dummy external input
  // Values not used but needed for initialisation of hsus
  NumericVector ext(2);

  // Create dummy timestep value
  double Dt = -99.0;
  
  // create vectors for storing lateral fluxes - make longer so can keep R id
  int maxid = max(id);
  //if(max(channel_id)>maxid){ maxid = max(channel_id); }
  maxid +=1;
  //NumericVector q_sf_in(maxid);
  //NumericVector q_sz_in(maxid);
  std::vector<double> q_sf_in(maxid), q_sz_in(maxid);
  std::vector<double> q_sf_out(maxid), q_sz_out(maxid);

  int nhillslope = states.nrow();
    
  // convert all flow directions to links - do first to avoid pushback issues?
  List fd;
  NumericVector fdfrc;
  IntegerVector fdidx;
  std::vector< std::vector<flink> > paul;
  std::vector< std::vector<flowlink> > links(nhillslope);
  int nfl;
  int idi;
  for(int ii=0; ii<states.nrow(); ++ii){
    idi = id[ii];
    fd = flow_dir(ii);
    fdidx = fd("idx");
    fdfrc = fd("frc");
    //nfl = fdi.length();
    std::vector<flink> tmp;
    //std::vector<flowlink> fl(nfl);
    for(int i=0; i<fdidx.length(); ++i){
      int j = fdidx[i];
      flink lnksf(fdfrc[i],q_sf_out[j],q_sf_in[j]);
      tmp.push_back( lnksf );
      flink lnksz(fdfrc[i],q_sz_out[j],q_sz_in[j]);
      tmp.push_back( lnksz );

      //fl[i].frc = &fdf[i];
      //fl[i].q_sf_in = &q_sf_in[j];
      //fl[i].q_sz_in = &q_sz_in[j];
    }
    paul.push_back(tmp);
    //links[ii] = fl;
  }
  
  // initialise the hsu
  std::vector<hsu> vhsu;

  //int idi = -99;

  for(int ii=0; ii<states.nrow(); ++ii){
    idi = id[ii];
    
    hsu h(states(ii,0),states(ii,1),states(ii,2),states(ii,3),
	  q_sf_in[idi], q_sz_in[idi],ext[0],ext[1],
	  properties(ii,0),properties(ii,1),properties(ii,2),
	  properties(ii,3),properties(ii,4),
	  properties(ii,5),
	  properties(ii,6),
	  properties(ii,7),properties(ii,8),
	  Dt, paul[ii], q_sf_out[idi], q_sz_out[idi]
	  );
    vhsu.push_back(h);
  }

  // loop hsu to initialise
  std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
  std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);
  for(int ii=0; ii<states.nrow(); ++ii){
    vhsu[ii].init(s_rz_0[ii],r_uz_sz_0[ii]);
    //List fd = flow_dir(ii);
    // NumericVector fdf = fd("frc");
    // IntegerVector fdi = fd("idx");
    // std::vector<double> q = vhsu[ii].get_q();
    // for(int jj=0; jj < fdi.length(); ++jj){
    //   int fdi_ = fdi[jj];
    //   double fdf_ = fdf[jj];
    //   q_sf_in[fdi_] += q[0]*fdf_;
    //   q_sz_in[fdi_] += q[1]*fdf_;
    // }
  }
}

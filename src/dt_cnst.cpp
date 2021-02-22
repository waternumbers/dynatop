#include "Rcpp.h"

// recall: Rcpp deep copies inputs unless they are the correct type!

// [[Rcpp::export]]
void dt_cnst(std::vector<int> id, // hillslope id
	       Rcpp::NumericMatrix states, // hillslope states
	       std::vector<double> area, // hillslope area - recall area/depth=width
	       std::vector<double> Dx, // hillslope length
	       std::vector<double> beta, // hillslope angle
	       std::vector<double> c_sf,
	       std::vector<double> k_sf,
	       std::vector<double> s_rzmax,
	       std::vector<double> t_d,
	       std::vector<double> m,
	       std::vector<double> ln_t0,
	       std::vector<int> channel_id,
	       std::vector<double> channel_area,
	       std::vector<int> p_idx,
	       std::vector<int> ep_idx,
	       std::vector<int> flow_from,
	       std::vector<int> flow_to,
	       std::vector<double> flow_frc,
	       Rcpp::NumericMatrix obs, // hillslope statesexternal series
	       Rcpp::NumericMatrix channel_inflow, // channel_inflow
	       Rcpp::NumericMatrix mass_balance, // mass balance for each timestep
	       std::vector<bool> keep_states,
	       Rcpp::List state_rec,
	       double timestep,
	       int n_sub_step,
	       bool approx_soln
	       ){
  // work out some dimensions
  int nhillslope = states.nrow();
  int nchannel = channel_id.size();
  int nlink = flow_from.size();

  // work out computational timestep - explicit casting of n_sub_step to double
  double Dt = timestep / (double)n_sub_step;
  
  // create vectors for storing lateral fluxes
  int maxid = std::max( *std::max_element(std::begin(id), std::end(id)),
			*std::max_element(std::begin(channel_id),std::end(channel_id)) );
  maxid +=1; // make longer so can keep R id
  std::vector<double> q_sf_in(maxid), q_sz_in(maxid);
  std::vector<double> q_sf_out(maxid), q_sz_out(maxid);

  // seperate out states to the hillslope
  Rcpp::NumericMatrix::Column l_sf = states.column(0); //( _ , 0);
  Rcpp::NumericMatrix::Column s_rz = states.column(1);//( _ , 1);
  Rcpp::NumericMatrix::Column s_uz = states.column(2);//( _ , 2);
  Rcpp::NumericMatrix::Column l_sz = states.column(3);//( _ , 3);
    
  // compute the property summaries for the hillslope
  std::vector<double> l_szmax(nhillslope),lambda_szmax(nhillslope),cosbeta_m(nhillslope);
  std::vector<double> lambda_sf(nhillslope), w(nhillslope);
  std::vector<double> s_sz(nhillslope), log_l_szmax(nhillslope), alpha(nhillslope);

  for(int i=0;i<nhillslope;++i){
    w[i] = area[i]/Dx[i]; // width
    cosbeta_m[i] = std::cos(beta[i]) /m[i]; //cosin of slope angle divided by m
    alpha[i] = cosbeta_m[i]*Dt/Dx[i]; // scaling factor for celerity calc
    lambda_sf[i] = (c_sf[i]*Dt)/Dx[i]; //constant for surface routing

    // // for exponential profile
    // l_szmax[i] = std::exp(ln_t0[i])*std::sin(beta[i]); // maximum flow from saturated zone
    // log_l_szmax[i] = ln_t0[i] + std::log( std::sin(beta[i]) ); // log of max flow from saturated zone
    // lambda_szmax[i] = cosbeta_m[i]*l_szmax[i]*Dt/Dx[i]; // max courant number 
    // s_sz[i] = (log_l_szmax[i] - std::log(l_sz[i]))/cosbeta_m[i];
    
    // for constant celerity with fixed max depth
    l_szmax[i] = 0.5*c_sf[i]; // maximum flow from saturated zone
    s_sz[i] = 0.5 - (l_sz[i] / c_sf[i]);
    
  }

  //Rcpp::Rcout << "initialised hillslope properties" << std::endl;
  
  // variable use with loop
  double p, ep;
  int cid,ip,iep;
  double r_sf_rz,r_rz_uz,r_uz_sz;
  double max_uz,r_uz_pos,r_uz_0, gamma, tuz, et, l_sz_tmp;
  double l_sf_in, l_sz_in, cin;
  std::vector<double> pbodge(12,-99.0);
  for(int it = 0; it < obs.nrow(); ++it) {
    //Rcpp::Rcout << it << std::endl;

    // work out the mass in each store
    mass_balance(it,0) = 0.0;
    for(int ii=0; ii<nhillslope; ++ii){
      // Rcpp::Rcout << l_sf[ii] <<std::endl;
      // Rcpp::Rcout << s_rz[ii] <<std::endl;
      // Rcpp::Rcout << s_uz[ii] <<std::endl;
      // Rcpp::Rcout << l_sz[ii] <<std::endl;
      // Rcpp::Rcout << l_sf[ii]/c_sf[ii] <<std::endl;
      // Rcpp::Rcout << s_sz[ii] <<std::endl;
      mass_balance(it,0) += area[ii]*((l_sf[ii]/c_sf[ii]) + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }
    mass_balance(it,1) = 0.0;
    mass_balance(it,2) = 0.0;
    mass_balance(it,3) = 0.0;
    
    // initialise channel inflow
    for(int ii=0; ii < nchannel; ++ii){
      channel_inflow(it,ii) = 0.0;
    }
    
    for(int nn = 0; nn < n_sub_step; ++nn){
      //Rcpp::Rcout << it << " " << nn << std::endl;
      std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
      std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);
      
      int link_cntr = 0; // counter for flow links
      int link_from_id = flow_from[link_cntr];
      
      for(int ii=0; ii<nhillslope; ++ii){


	// current id used to reference longer vectors
	cid = id[ii];

	// inputs
	ip = p_idx[cid];
	p = obs(it,ip) / timestep;
	iep = ep_idx[cid];
	ep = obs(it,iep) / timestep;

	// standardise inflows by width
	l_sf_in = q_sf_in[cid] / w[ii];
	l_sz_in = q_sz_in[cid] / w[ii];

	// max downward flux times (superscript pos direction)
	r_sf_rz = std::min( k_sf[ii] , (l_sf[ii] + lambda_sf[ii]*l_sf_in)/(Dx[ii]*lambda_sf[ii]) );
	r_rz_uz = std::max( 0.0 , ((s_rz[ii]-s_rzmax[ii])/Dt) + p + r_sf_rz - ep );
	// max volume in unsat from drain down
	max_uz = s_uz[ii] + Dt*r_rz_uz;
	
	// explicit estimates of r_uz_sz
	// use approximate solution  based on previous storage
	r_uz_sz = std::min(max_uz / (t_d[ii]*s_sz[ii] + Dt),1.0/t_d[ii]);

	// store the debugging states at the start
	pbodge[0] = s_uz[ii];
	pbodge[1] = r_uz_sz ;
	pbodge[2] = s_sz[ii];
	pbodge[3] = l_sz[ii];
	pbodge[4] = l_szmax[ii];
	// flow between stores after downwards pass
	//Rcpp::Rcout << "r_sf_rz " << r_sf_rz << std::endl;
	//Rcpp::Rcout << "r_rz_uz " << r_rz_uz << std::endl;
	//Rcpp::Rcout << "r_uz_sz " << r_uz_sz << std::endl;

	// solve l_sz assuming constant celerity
	gamma = lambda_sf[ii]; // take constant so that mass balance is ensured
	// initial estimate
	l_sz_tmp = (l_sz[ii] + gamma*(l_sz_in + Dx[ii]*r_uz_sz)) / (1.0+gamma);
	

	pbodge[5] = gamma;
	pbodge[6] = l_sz_tmp;
	pbodge[7] = gamma*l_sz_in;
	pbodge[8] = (1.0+gamma)*l_szmax[ii] ;//(1.0+gamma)*l_szmax[ii];
	pbodge[9] = 0.0;
	pbodge[12] = gamma*Dx[ii];
	
	// handle saturation
	if( l_sz_tmp > l_szmax[ii] ){
	  pbodge[9] = 1.0;
	  // this calc is wrong for some reason....
	  //manually using values returns the correct solution
	  r_uz_sz = ( ((1.0+gamma)*l_szmax[ii]) - l_sz[ii] - (gamma*l_sz_in)) / (gamma*Dx[ii]); //pbodge[12];
	  
	  l_sz_tmp = l_szmax[ii];
	}

	pbodge[10] = r_uz_sz;
	pbodge[11] = l_sz_tmp;
	// // explicit solution to the saturated zone
	// gamma = (cosbeta_m[ii]*l_sz[ii])*Dt/Dx[ii]; // lambda with explicit celerity
	// double l_sz_tmp = (l_sz[ii] + gamma*(l_sz_in + Dx[ii]*r_uz_sz)) / (1.0+gamma);
	// // handle saturation
	// if( l_sz_tmp > l_szmax[ii] ){
	//   r_uz_sz = ((1.0+gamma)*l_szmax[ii] - l_sz[ii] - gamma*l_sf_in )/Dx[ii];
	//   l_sz_tmp = l_szmax[ii];
	// }

	// // semi-implict solution
	// gamma = alpha[ii]*(l_sz_in + Dx[ii]*r_uz_sz) - 1;
	// l_sz[ii] = ( gamma + std::sqrt( std::pow(gamma,2) + 4.0*alpha[ii]*l_sz[ii]) )/(2.0*alpha[ii]);
	// // handle saturation
	// if( l_sz[ii] > l_szmax[ii] ){
	//   r_uz_sz = ( (1+lambda_szmax[ii])*l_szmax[ii] - l_sz[ii] - lambda_szmax[ii]*l_sz_in ) /(lambda_szmax[ii]*Dx[ii]);
	//   l_sz[ii] = l_szmax[ii];
	// }

	l_sz[ii] = l_sz_tmp;
	// for constant celerity
	s_sz[ii] = 0.5 - (l_sz[ii] / c_sf[ii]);
	
	//s_sz[ii] = (log_l_szmax[ii] - std::log(l_sz[ii]))/cosbeta_m[ii];
	//Rcpp::Rcout << "s_sz is " << s_sz[ii] << std::endl;
	
	
	// solve unsaturated zone
	tuz = std::min( s_sz[ii], s_uz[ii] + Dt*(r_rz_uz -r_uz_sz) );

	if(tuz < -1e-16){
	  Rcpp::Rcout << "iteration " << it << std::endl;
	  Rcpp::Rcout << "ii " << ii << std::endl;
	  Rcpp::Rcout << "id " << cid << std::endl;

	  Rcpp::Rcout << "Initial states " << std::endl;
	  Rcpp::Rcout << "s_uz " << pbodge[0] << std::endl;
	  Rcpp::Rcout << "r_uz_sz " << pbodge[1] << std::endl;
	  Rcpp::Rcout << "s_sz " << pbodge[2] << std::endl;
	  Rcpp::Rcout << "l_sz " << pbodge[3] << std::endl;
	  Rcpp::Rcout << "l_szmax " << pbodge[4] << std::endl;
	  
	  Rcpp::Rcout << "Initial Computed values " << std::endl;
	  Rcpp::Rcout << "gamma " << pbodge[5] << std::endl;
	  Rcpp::Rcout << "l_sz_tmp " << pbodge[6] << std::endl;
	  Rcpp::Rcout << "gamma * l_sz_in " << pbodge[7] << std::endl;
	  Rcpp::Rcout << "(1.0+gamma) * l_sz_max " << pbodge[8] << std::endl;
	  Rcpp::Rcout << "gamma * Dx " << pbodge[12] << std::endl;
	  
	  Rcpp::Rcout << "After adjustment " << std::endl;
	  Rcpp::Rcout << "Was saturation adjustment called? " << pbodge[9] << std::endl;
	  Rcpp::Rcout << "r_uz_sz " << pbodge[10] << std::endl;
	  Rcpp::Rcout << "l_sz_tmp " <<pbodge[11] << std::endl;

	  Rcpp::Rcout << "Final states " << std::endl;
	  Rcpp::Rcout << "s_uz " << tuz << std::endl;
	  Rcpp::Rcout << "l_sz " << l_sz[ii] << std::endl;
	  Rcpp::Rcout << "l_sz " << s_sz[ii] << std::endl;
	  
	}
	
	//Rcpp::Rcout << "s_uz tuz " << s_uz[ii] << " " << tuz << std::endl;
	r_rz_uz = (tuz - s_uz[ii])/Dt + r_uz_sz;
	s_uz[ii] = tuz;

	//Rcpp::Rcout << "solved uz" << std::endl;
	//Rcpp::Rcout << "corrected r_rz_uz " << r_rz_uz << std::endl;

	// solve root zone
	tuz = s_rz[ii];
	s_rz[ii] = std::min( s_rzmax[ii] , (s_rz[ii] + Dt*(p + r_sf_rz - r_rz_uz)) /
			     (1.0+((ep*Dt)/s_rzmax[ii])) );
	r_sf_rz = ( (1.0+((ep*Dt)/s_rzmax[ii]))*s_rz[ii] - tuz - Dt*p + Dt*r_rz_uz ) /Dt ;
	
	// solve the surface
	l_sf[ii] = (l_sf[ii] + lambda_sf[ii]*(l_sf_in - Dx[ii]*r_sf_rz)) / (1 + lambda_sf[ii]);
	
	// add incoming and outgoing vertical fluxes to mass balance
	et = ep*(s_rz[ii]/s_rzmax[ii]);
	
	mass_balance(it,1) -= et*area[ii]*Dt;
	mass_balance(it,2) += p*area[ii]*Dt;
	
	//Rcpp::Rcout << "solved r_sf_rz " << r_sf_rz << std::endl;
	//Rcpp::Rcout << "solved r_rz_uz " << r_rz_uz << std::endl;
	//Rcpp::Rcout << "solved r_uz_sz " << r_uz_sz << std::endl;
	// tranfer on the outflow
	//Rcpp::Rcout << "flows for cid " << cid << std::endl;
	while( (link_from_id == cid) & (link_cntr<nlink) ){
	  int j = flow_to[link_cntr];
	  double f = flow_frc[link_cntr];
	  //Rcpp::Rcout << "from " << link_from_id << std::endl;
	  //Rcpp::Rcout << "to " << j << std::endl;
	  //Rcpp::Rcout << "frac " << f << std::endl;
	  //Rcpp::Rcout << "width " << w[ii] << std::endl;
	  //Rcpp::Rcout << "cntr " << link_cntr << std::endl;
	  //Rcpp::Rcout << "value " << f*w[ii]*l_sz[ii] << std::endl;
	  q_sf_in[j] += f*w[ii]*l_sf[ii];
	  q_sz_in[j] += f*w[ii]*l_sz[ii]; //l_sz[ii];
	  link_cntr +=1;
	  if(link_cntr < nlink) {
	    link_from_id = flow_from[link_cntr];
	  }
	}
	
	// end of hillslope
      }

      // loop channels for the sub step
      for(int ii=0; ii < nchannel; ++ii){
	cid = channel_id[ii];
	ip = p_idx[cid];
	p = obs(it,ip)/Dt;
	cin = (channel_area[ii]*p + q_sf_in[cid] + q_sz_in[cid])*Dt;
	//Rcpp::Rcout << "channel id " << cid << std::endl;
	//Rcpp::Rcout << "channel ip " << ip << std::endl;
	//Rcpp::Rcout << "q_sf_in " << q_sf_in[cid] << std::endl;
	//Rcpp::Rcout << "q_sz_in " <<  q_sz_in[cid] << std::endl;
	//Rcpp::Rcout << "obs(it,ip) " <<  obs(it,ip) << std::endl;
	mass_balance(it,2) += p*channel_area[ii]*Dt;  
	channel_inflow(it,ii) += cin;
	mass_balance(it,3) -= cin;
      }
    }
    
    // convert channel inflow to rate
    for(int ii=0; ii < nchannel; ++ii){
      channel_inflow(it,ii) = channel_inflow(it,ii)/timestep;
    }
    
    if( keep_states[it] ){
      state_rec(it) = Rcpp::clone(states);
    }

    Rcpp::checkUserInterrupt(); 
  }
}


// // [[Rcpp::export]]
// void multi_hsu_cpp_init(IntegerVector id,
// 			NumericMatrix states,
// 			NumericMatrix properties,
// 			List flow_dir,
// 			NumericVector s_rz_0,
// 			NumericVector r_uz_sz_0
// 			)
// {

//   // Create dummy external input
//   // Values not used but needed for initialisation of hsus
//   NumericVector ext(2);

//   // Create dummy timestep value
//   double Dt = -99.0;
  
//   // create vectors for storing lateral fluxes - make longer so can keep R id
//   int maxid = max(id);
//   //if(max(channel_id)>maxid){ maxid = max(channel_id); }
//   maxid +=1;
//   //NumericVector q_sf_in(maxid);
//   //NumericVector q_sz_in(maxid);
//   std::vector<double> q_sf_in(maxid), q_sz_in(maxid);
//   std::vector<double> q_sf_out(maxid), q_sz_out(maxid);

//   int nhillslope = states.nrow();
    
//   // convert all flow directions to links - do first to avoid pushback issues?
//   List fd;
//   NumericVector fdfrc;
//   IntegerVector fdidx;
//   std::vector< std::vector<flink> > paul;
//   std::vector< std::vector<flowlink> > links(nhillslope);
//   int nfl;
//   int idi;
//   for(int ii=0; ii<states.nrow(); ++ii){
//     idi = id[ii];
//     fd = flow_dir(ii);
//     fdidx = fd("idx");
//     fdfrc = fd("frc");
//     //nfl = fdi.length();
//     std::vector<flink> tmp;
//     //std::vector<flowlink> fl(nfl);
//     for(int i=0; i<fdidx.length(); ++i){
//       int j = fdidx[i];
//       flink lnksf(fdfrc[i],q_sf_out[j],q_sf_in[j]);
//       tmp.push_back( lnksf );
//       flink lnksz(fdfrc[i],q_sz_out[j],q_sz_in[j]);
//       tmp.push_back( lnksz );

//       //fl[i].frc = &fdf[i];
//       //fl[i].q_sf_in = &q_sf_in[j];
//       //fl[i].q_sz_in = &q_sz_in[j];
//     }
//     paul.push_back(tmp);
//     //links[ii] = fl;
//   }
  
//   // initialise the hsu
//   std::vector<hsu> vhsu;

//   //int idi = -99;

//   for(int ii=0; ii<states.nrow(); ++ii){
//     idi = id[ii];
    
//     hsu h(states(ii,0),states(ii,1),states(ii,2),states(ii,3),
// 	  q_sf_in[idi], q_sz_in[idi],ext[0],ext[1],
// 	  properties(ii,0),properties(ii,1),properties(ii,2),
// 	  properties(ii,3),properties(ii,4),
// 	  properties(ii,5),
// 	  properties(ii,6),
// 	  properties(ii,7),properties(ii,8),
// 	  Dt, paul[ii], q_sf_out[idi], q_sz_out[idi]
// 	  );
//     vhsu.push_back(h);
//   }

//   // loop hsu to initialise
//   std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
//   std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);
//   for(int ii=0; ii<states.nrow(); ++ii){
//     vhsu[ii].init(s_rz_0[ii],r_uz_sz_0[ii]);
//     //List fd = flow_dir(ii);
//     // NumericVector fdf = fd("frc");
//     // IntegerVector fdi = fd("idx");
//     // std::vector<double> q = vhsu[ii].get_q();
//     // for(int jj=0; jj < fdi.length(); ++jj){
//     //   int fdi_ = fdi[jj];
//     //   double fdf_ = fdf[jj];
//     //   q_sf_in[fdi_] += q[0]*fdf_;
//     //   q_sz_in[fdi_] += q[1]*fdf_;
//     // }
//   }
// }

#include "Rcpp.h"

// recall: Rcpp deep copies inputs unless they are the correct type!

// [[Rcpp::export]]
void dt_init(std::vector<int> id, // hillslope id
	     Rcpp::NumericMatrix states, // hillslope states
	     std::vector<double> area, // hillslope area
	     std::vector<double> Dx, // hillslope length
	     std::vector<double> beta, // hillslope angle
	     std::vector<double> c_sf,
	     std::vector<double> s_rzmax,
	     std::vector<double> t_d,
	     std::vector<double> m,
	     std::vector<double> ln_t0,
	     std::vector<int> channel_id, // need to allow full passing
	     std::vector<int> flow_from,
	     std::vector<int> flow_to,
	     std::vector<double> flow_frc,
	     std::vector<double> s_rz_0, // initial root zone depth as fraction
	     std::vector<double> r_uz_sz_0 // initial root zone depth as fraction	     
	     ){
  // work out some dimensions
  int nhillslope = states.nrow();
  int nlink = flow_from.size();
   
  // create vectors for storing lateral fluxes
  int maxid = std::max( *std::max_element(std::begin(id), std::end(id)),
			*std::max_element(std::begin(channel_id),std::end(channel_id)) );
  //maxid +=1; // make longer so can keep R id
  std::vector<double> q_sf_in(maxid,0.0), q_sz_in(maxid,0.0);
  std::vector<double> q_sf_out(maxid,0.0), q_sz_out(maxid,0.0);
  
  // seperate out states to the hillslope
  Rcpp::NumericMatrix::Column l_sf = states.column(0); //( _ , 0);
  Rcpp::NumericMatrix::Column s_rz = states.column(1);//( _ , 1);
  Rcpp::NumericMatrix::Column s_uz = states.column(2);//( _ , 2);
  Rcpp::NumericMatrix::Column l_sz = states.column(3);//( _ , 3);
  
  // compute the property summaries for the hillslope
  std::vector<double> l_szmax(nhillslope,-99.0), cosbeta_m(nhillslope,-99.0),
    log_l_szmax(nhillslope,-99.0);
  
  for(int i=0;i<nhillslope;++i){
    // // for exponential profile
    l_szmax[i] = std::exp(ln_t0[i])*std::sin(beta[i]);
    log_l_szmax[i] = ln_t0[i] + std::log( std::sin(beta[i]) );
    cosbeta_m[i] = std::cos(beta[i]) /m[i];

    // for constant celerity with fixed max depth
    //l_szmax[i] = 0.5 * c_sf[i];
  }
  
  //Rcpp::Rcout << "initialised hillslope properties" << std::endl;

  // variable use with loop
  int cid;
  double s_sz;
  double l_sz_in;
  double w;
  double r_uz_sz;
    
  int link_cntr = 0; // counter for flow links
  int link_from_id = flow_from[link_cntr];
  
  for(int ii=0; ii<nhillslope; ++ii){
    // current id used to reference longer vectors
    cid = id[ii];

    // compute width
    w = area[ii]/Dx[ii];
    
    // standardise inflows by width
    l_sz_in = q_sz_in[cid] / w;

    // apply the initial value to the surface and root zones
    l_sf[ii] = 0.0;
    s_rz[ii] = s_rzmax[ii]*s_rz_0[ii];

        
    // out flux under steady state
    l_sz[ii] = std::min(l_szmax[ii], l_sz_in + Dx[ii]*r_uz_sz_0[ii]);

    // solve to find actual r_uz_sz
    r_uz_sz = (l_sz[ii] - l_sz_in)/Dx[ii];

    
    // for exponential
    s_sz = std::max(0.0, (log_l_szmax[ii] - std::log(l_sz[ii]))/cosbeta_m[ii]);
    // for constant celerity with max depth of 0.5
    //s_sz = 0.5 - (l_sz[ii] / c_sf[ii]);
    
    s_uz[ii] = std::min( s_sz, r_uz_sz*t_d[ii]*s_sz );
    //cpp::Rcout << l_sz[ii] << " " << s_sz << " " << s_uz[ii] << std::endl;

    // tranfer on the outflow    
    while( (link_from_id == cid) & (link_cntr<nlink) ){
      int& j = flow_to[link_cntr];
      double& f = flow_frc[link_cntr];
      q_sf_in[j] += f*w*l_sf[ii];
      q_sz_in[j] += f*w*l_sz[ii];
      link_cntr +=1;
      if(link_cntr < nlink) {
  	link_from_id = flow_from[link_cntr];
      }
    }
  }
}


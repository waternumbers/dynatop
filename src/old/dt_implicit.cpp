#include "Rcpp.h"
#include <boost/integer/common_factor.hpp>  
#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>
// recall: Rcpp deep copies inputs unless they are the correct type!


// structure for implicit solution
template <class T>
struct fsz {
  fsz(double const& muz_, double const& l_,
      double const& lin_,
      double const& dx_, double const& dt_,
      double const& td_, double const& cbm_,
      double const& llmax_, double const& lmax_) :
    muz(muz_), l(l_), lin(lin_),
    dx(dx_), dt(dt_),
    td(td_), cbm(cbm_), llmax(llmax_), lmax(lmax_)
  { /* Constructor just stores value a to find root of. */ }
  double operator()(double const& x){
    double sz = (llmax - std::log(x))/cbm;
    double ruz = std::min( muz/(td*sz + dt), 1/td);
    double lambda = cbm*x*dt/dx;
    //Rcpp::Rcout << "x, lambda " << x << " " << lambda << std::endl;
    return x - std::min(lmax,(l + lambda*(lin + dx*ruz))/(1.0+lambda));
  }
private:
  double muz, l, lin, dx,dt, td, cbm,llmax,lmax; // to be 'cube_rooted'.
};


// [[Rcpp::export]]
void dt_implicit(std::vector<int> id, // hillslope id
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
		 int n_sub_step
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
    l_szmax[i] = std::exp(ln_t0[i])*std::sin(beta[i]); // maximum flow from saturated zone
    log_l_szmax[i] = ln_t0[i] + std::log( std::sin(beta[i]) ); // log of max flow from saturated zone
    lambda_szmax[i] = cosbeta_m[i]*l_szmax[i]*Dt/Dx[i]; // max courant number 
    s_sz[i] = (log_l_szmax[i] - std::log(l_sz[i]))/cosbeta_m[i];
    
    // for constant celerity with fixed max depth
    //l_szmax[i] = 0.5*c_sf[i]; // maximum flow from saturated zone
    //s_sz[i] = 0.5 - (l_sz[i] / c_sf[i]);
    
  }


  // boost optimisation parameters
  const boost::uintmax_t opt_maxit = 1000;
  int digits = std::numeric_limits<double>::digits;
  int get_digits = (digits * 3) /4;
  boost::math::tools::eps_tolerance<double> tol(get_digits);
  boost::uintmax_t opt_it = opt_maxit;
  std::pair<double, double> opt_res; // solution for output
  



  
  //Rcpp::Rcout << "initialised hillslope properties" << std::endl;
  
  // variable use with loop
  double p, ep;
  int cid,ip,iep;
  double r_sf_rz,r_rz_uz,r_uz_sz;
  double max_uz, tuz, et, r_uz_sz_pos, r_uz_sz_0;
  double l_sf_in, l_sz_in, cin;
  std::vector<double> mbv(4,0.0); // mass balance vector
  std::vector<double> ch_in(nchannel,0.0); // channel inflow vector
 
  for(int it = 0; it < obs.nrow(); ++it) {
    //Rcpp::Rcout << it << std::endl;
    std::fill(mbv.begin(), mbv.end(), 0.0);
    // work out the mass in each store
    for(int ii=0; ii<nhillslope; ++ii){
      mbv[0] += area[ii]*((l_sf[ii]/c_sf[ii]) + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    }
    // // work out the mass in each store
    // mass_balance(it,0) = 0.0;
    // for(int ii=0; ii<nhillslope; ++ii){
    //   // Rcpp::Rcout << l_sf[ii] <<std::endl;
    //   // Rcpp::Rcout << s_rz[ii] <<std::endl;
    //   // Rcpp::Rcout << s_uz[ii] <<std::endl;
    //   // Rcpp::Rcout << l_sz[ii] <<std::endl;
    //   // Rcpp::Rcout << l_sf[ii]/c_sf[ii] <<std::endl;
    //   // Rcpp::Rcout << s_sz[ii] <<std::endl;
    //   mass_balance(it,0) += area[ii]*((l_sf[ii]/c_sf[ii]) + s_rz[ii] + s_uz[ii] - s_sz[ii]);
    // }
    // mass_balance(it,1) = 0.0;
    // mass_balance(it,2) = 0.0;
    // mass_balance(it,3) = 0.0;
    // initialise channel inflow
    std::fill(ch_in.begin(), ch_in.end(), 0.0);
    // for(int ii=0; ii < nchannel; ++ii){
    //   channel_inflow(it,ii) = 0.0;
    // }
    
    for(int nn = 0; nn < n_sub_step; ++nn){
      //Rcpp::Rcout << it << " " << nn << std::endl;
      std::fill(q_sf_in.begin(), q_sf_in.end(), 0.0);
      std::fill(q_sz_in.begin(), q_sz_in.end(), 0.0);

      
      int link_cntr = 0; // counter for flow links
      int link_from_id = flow_from[link_cntr];
      
      for(int ii=0; ii<nhillslope; ++ii){

	//double& paul_w = w[ii];
	//double& paul_Dx = Dx[ii];
	
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

	// see if we need to optimise
	r_uz_sz_pos = std::min( max_uz/Dt, 1/t_d[ii]);
	r_uz_sz_0 = ((1+lambda_szmax[ii])*l_szmax[ii] - l_sz[ii] - lambda_szmax[ii]*l_sz_in)/(Dx[ii]*lambda_szmax[ii]);
	if( r_uz_sz_pos >= r_uz_sz_0 ){
	  l_sz[ii] = l_szmax[ii];
	  s_sz[ii] = 0.0;
	  r_uz_sz = r_uz_sz_0;
	}else{
	  // set up for optimisation
	  opt_it = opt_maxit;
	  
	  fsz<double> oc(max_uz, l_sz[ii],l_sz_in,
			 Dx[ii],Dt,
			 t_d[ii],cosbeta_m[ii],
			 log_l_szmax[ii], l_szmax[ii]);
	  
	  double fcr = oc(l_sz[ii]);
	  //Rcpp::Rcout << "fcr is " << fcr << std::endl;
	  if(fcr!=0.0){
	    // initial range estimates
	    double upr = l_sz[ii], lwr = 0.0;
	    if(fcr < 0.0){
	      lwr = l_sz[ii];
	      upr = l_szmax[ii];
	    }
	    if( oc(lwr)*oc(upr) > 0.0 ){
	      Rcpp::Rcout << "limits " << it << " " << ii << " " << lwr << " " << upr << " " << l_szmax[ii] << " " << fcr << std::endl;
	      Rcpp::Rcout << "values " << oc(lwr) << " " << oc(upr) << std::endl;
	    }
	    opt_res = boost::math::tools::toms748_solve(oc, lwr, upr, tol,opt_it);
	    if(opt_it >= opt_maxit)
	      { // commented out to pass CRAN tests 
		//std::cout << "Unable to locate solution in chosen iterations:"
		//  " Current best guess is between " << r.first << " and " << r.second << std::endl;
	      }
	    l_sz[ii] = (opt_res.second + opt_res.first)/2.0;
	    // handle rounding errors
	    l_sz[ii] = std::min(l_szmax[ii],l_sz[ii]);
	    l_sz[ii] = std::max(l_sz[ii],0.0);
	  }
	  s_sz[ii] = (log_l_szmax[ii] - std::log(l_sz[ii]))/cosbeta_m[ii];
	  r_uz_sz = std::min( max_uz/(t_d[ii]*s_sz[ii] + Dt), 1/t_d[ii]);
	}	
	
	// solve unsaturated zone
	tuz = std::min( s_sz[ii], s_uz[ii] + Dt*(r_rz_uz -r_uz_sz) );

	// if(tuz < -1e-16){
	//   Rcpp::Rcout << "iteration " << it << std::endl;
	//   Rcpp::Rcout << "ii " << ii << std::endl;
	//   Rcpp::Rcout << "id " << cid << std::endl;

	//   Rcpp::Rcout << "Initial states " << std::endl;
	//   Rcpp::Rcout << "s_uz " << pbodge[0] << std::endl;
	//   Rcpp::Rcout << "r_uz_sz " << pbodge[1] << std::endl;
	//   Rcpp::Rcout << "s_sz " << pbodge[2] << std::endl;
	//   Rcpp::Rcout << "l_sz " << pbodge[3] << std::endl;
	//   Rcpp::Rcout << "l_szmax " << pbodge[4] << std::endl;
	  
	//   Rcpp::Rcout << "Initial Computed values " << std::endl;
	//   Rcpp::Rcout << "gamma " << pbodge[5] << std::endl;
	//   Rcpp::Rcout << "l_sz_tmp " << pbodge[6] << std::endl;
	//   Rcpp::Rcout << "gamma * l_sz_in " << pbodge[7] << std::endl;
	//   Rcpp::Rcout << "(1.0+gamma) * l_sz_max " << pbodge[8] << std::endl;
	//   Rcpp::Rcout << "gamma * Dx " << pbodge[12] << std::endl;
	  
	//   Rcpp::Rcout << "After adjustment " << std::endl;
	//   Rcpp::Rcout << "Was saturation adjustment called? " << pbodge[9] << std::endl;
	//   Rcpp::Rcout << "r_uz_sz " << pbodge[10] << std::endl;
	//   Rcpp::Rcout << "l_sz_tmp " <<pbodge[11] << std::endl;

	//   Rcpp::Rcout << "Final states " << std::endl;
	//   Rcpp::Rcout << "s_uz " << tuz << std::endl;
	//   Rcpp::Rcout << "l_sz " << l_sz[ii] << std::endl;
	//   Rcpp::Rcout << "l_sz " << s_sz[ii] << std::endl;
	  
	// }
	
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
	mbv[1] -= et*area[ii]*Dt;
	mbv[2] += p*area[ii]*Dt;

	// tranfer on the outflow
	double q_sf_out = w[ii]*l_sf[ii];
	double q_sz_out = w[ii]*l_sz[ii];
	while( (link_from_id == cid) & (link_cntr<nlink) ){
	  int j = flow_to[link_cntr];
	  double f = flow_frc[link_cntr];
	  q_sf_in[j] += f*q_sf_out; //w[ii]*l_sf[ii];
	  q_sz_in[j] += f*q_sz_out; //w[ii]*l_sz[ii];
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
      	mbv[2] += p*channel_area[ii]*Dt; 
      	ch_in[ii] += cin;
      	mbv[3] -= cin;
      }
    }

    for(int ii=0; ii < 4; ++ii){
      mass_balance(it,ii) = mbv[ii];
    }
    
    // convert channel inflow to rate
    for(int ii=0; ii < nchannel; ++ii){
      channel_inflow(it,ii) = ch_in[ii]/timestep; //channel_inflow(it,ii)/timestep;
    }
    
    if( keep_states[it] ){
      state_rec(it) = Rcpp::clone(states);
    }

    Rcpp::checkUserInterrupt(); 
  }
}

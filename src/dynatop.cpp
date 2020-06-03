#include "Rcpp.h"
using namespace Rcpp;

double fode(double a, double b, double x0,double t){
  // Rcout << "a is: " << a << std::endl;
  // Rcout << "b is: " << b << std::endl;
  // Rcout << "x0 is: " << x0 << std::endl;
  // Rcout << "t is: " << t << std::endl;
  double ebt = exp(-b*t);
  double kappa = (1-ebt)/b;
  if( b == 0.0 ){
    kappa = t;
  }
  // Rcout << "b is: " << b << std::endl;
  // Rcout << "ebt is: " << ebt << std::endl;
  // Rcout << "kappa is: " << kappa << std::endl;
  double x = x0*ebt + a*kappa;
  return(x);
}


// [[Rcpp::export]]
void hs_sim_cpp(NumericMatrix channel_inflow, NumericMatrix mass_errors,
		List state_record,
		const NumericMatrix obs, const LogicalVector keep_states,
		bool mass_check,
		List hillslope, const List channel, const List ts,
		const List sz_dir, const List sf_dir, const List sqnc){
  
  // seperate out the hillslope
  // states
  NumericVector s_sf = hillslope("s_sf"), s_rz = hillslope("s_rz"),
    s_uz = hillslope("s_uz"), s_sz = hillslope("s_sz"), l_sz = hillslope("l_sz");
  // parameters & properties
  IntegerVector id = hillslope("id"), pcol = hillslope("precip"),
    ecol = hillslope("pet");
  //id = id - 1; pcol = pcol - 1; ecol = ecol - 1; // adjustment of 0 start in indexing
  
  const NumericVector t_sf = hillslope("t_sf"),q_sfmax = hillslope("q_sfmax"),
    s_rzmax = hillslope("s_rzmax"), t_d = hillslope("t_d"), 
    s_bar = hillslope("s_bar"), m = hillslope("m"), ln_t0 = hillslope("ln_t0"),
    delta_x = hillslope("delta_x"), area = hillslope("area");

  // seperate out the channel
  // properties
  const IntegerVector cid = channel("id"), cpcol = channel("precip");
  //cid = cid - 1; cpcol = cpcol - 1; // adjustment of 0 start in indexing
  
  const NumericVector carea = channel("area");

  // seperate out the time variables
  const int n_sub_step = ts("n_sub_step");
  const double step = ts("step"), sub_step = ts("sub_step");

  // Rcout << "sub_step" << sub_step << std::endl;
  
  // seperate out the sequnces
  const IntegerVector sqnc_sf = sqnc("sf"), sqnc_sz = sqnc("sz");
  // sqnc_sf = sqnc_sf - 1; sqnc_sz = sqnc_sz - 1;  // adjustment of 0 start in indexing
  
  // dimensions
  const int nit = obs.nrow(), nhsu = id.size(),
    nch = cid.size(), max_id = max(id);
		     
  // local constants used in the evaluation
  const NumericVector beta = atan( s_bar );
  const NumericVector cosbeta_m = cos(beta) / m;
  const NumericVector sinbeta = sin(beta);
  const NumericVector l_szmax = exp( ln_t0 )* sinbeta / delta_x;

  // declare storages for lateral flows
  NumericVector il_sf_in(max_id), il_sz_in(max_id);

  // declare storages for inputs
  //NumericVector p(nhsu), e_p(nhsu), cp(nch);// inputs
  //NumericVector orw(obs.ncol());
  
  // declare variables used in loops
  int ii, ii_id, iidx;
  List dir;
  IntegerVector idx;
  NumericVector frc;
  
    
  // declare variables used in surface
  double il_sf;

  // declare variables used in rootzone
  double p, e_p, q_sf_rz, ie_t, iq_rz_sf, iq_rz_uz;
  
  // declare variables used in unsaturated zone
  double iq_uz_sz, iq_uz_sf;
    
  // declare variables used in saturated zone
  double il_sz, bq_uz_sz, lbq_uz_sz, l_szin, lbar, lambda, lp, iq_sz_sf;

  // declare data frame of states
  // these are pointers so update and need to be cloned to store
  DataFrame df_state = DataFrame::create( Named("id") = id,
					  Named("id") = s_sf,
					  Named("id") = s_rz,
					  Named("id") = s_uz,
					  Named("id") = s_sz,
					  Named("id") = l_sz);

  
  // loop time steps
  for(int it =0; it < nit; it++) {
    
    // initialise the mass_errors and channel_inflow
    if( mass_check){
      mass_errors(it,0) = sum( (s_sf + s_rz + s_uz - s_sz)*area ); // initial states
      mass_errors(it,2) = 0.0; // precipitation input - added to in code
      mass_errors(it,3) = 0.0; // e_t - added to in code
    }
    
    // loop channel elements to initialise the flow to channel
    for(int ii =0; ii < nch; ii++){
      ii_id = cid(ii) - 1;
      p = obs(it,cpcol(ii)-1);
      channel_inflow(it,ii) = p;
      if( mass_check ){
	mass_errors(it,2) += p*carea(ii);
      }
    }
 
    // loop the sub steps - ss does nothing
    for(int ss=0; ss < n_sub_step; ss++){
      
      // declare storages for lateral flows as zero
      il_sf_in.fill(0.0);
      il_sz_in.fill(0.0);

      //Rcout << "Start Surface" <<std::endl;
      
      // loop hsus for the surface store
      for(int sq = 0; sq < nhsu; sq++){
	// ii is location of HSU in vector
	ii = sqnc_sf(sq) - 1;
	// ii_id is the HSU id
	ii_id = id(ii) - 1;

	// initialise flow estiamte
	il_sf = s_sf(ii) + il_sf_in(ii_id);	
	// solve ODE
	s_sf(ii) = fode(il_sf_in(ii_id) / sub_step,
			sub_step / t_sf(ii),
			s_sf(ii),
			sub_step);	
	// finalise outflow estimate
	il_sf = il_sf - s_sf(ii);
	// pass downslope
	dir = sf_dir(ii);
	idx = dir("idx");
	frc = dir("frc");
	for(int jj =0; jj < idx.size(); jj++){
	  iidx = idx(jj) - 1;
	  il_sf_in(iidx) = il_sf_in(iidx) + frc(jj)*il_sf;
	}
      }

      // loop hsus for the subsurface
      for(int sq = 0; sq < nhsu; sq++){
	// ii is location of HSU in vector
	ii = sqnc_sz(sq) - 1;
	// ii_id is the HSU id
	ii_id = id(ii) - 1;


	
	// Step2: solve the rootzone
	// work out inputs
	p = obs(it,pcol(ii)-1) / step;
	e_p = obs(it,ecol(ii)-1) / step;
	// determine inflow for surface
	q_sf_rz = std::min( s_sf(ii) / sub_step, q_sfmax(ii) );
	// correct surface
	s_sf(ii) = s_sf(ii) - (q_sf_rz*sub_step);
	// initialise estimate of e_t
	ie_t =  s_rz(ii) + ((p+q_sf_rz)*sub_step);
	// update storage estimate
	s_rz(ii) = fode(p+q_sf_rz,
			e_p / s_rzmax(ii),
			s_rz(ii),sub_step);
	// finalise estimate of evap-transpiration
	ie_t = ie_t - s_rz(ii);
	// add p and actual e_t into mass check
	if( mass_check ){
	  mass_errors(it,2) = mass_errors(it,2) + (p*sub_step*area(ii));
	  mass_errors(it,3) = mass_errors(it,3) + (ie_t*area(ii));
	}
	// work out flow - presume goes to unsat
	iq_rz_uz = std::max( s_rz(ii) - s_rzmax(ii), 0.0);
	iq_rz_sf = 0.0;
	// correct storage
	s_rz(ii) = s_rz(ii) - iq_rz_uz;
	// all for saturation with flow to surface
	if( s_sz(ii) <= 0 ){
	  iq_rz_sf = iq_rz_uz;
	  iq_rz_uz = 0.0;
	}
	
	// Step 3: Unsaturated zone
	// initialise outflow
	iq_uz_sz = s_uz(ii) + iq_rz_uz;
	// solve for new state
	s_uz(ii) = fode( iq_rz_uz/sub_step,
			 1.0 / (t_d(ii) * s_sz(ii)),
			 s_uz(ii),sub_step );
	// finalise outflow
	iq_uz_sz = iq_uz_sz - s_uz(ii);
	iq_uz_sf = 0.0;

	// Step 4: Solve saturated zone
	// this is omega=theta=1 solution with
	// explicit velocity and average inflow                    
	
	il_sz = (sub_step/2.0)*l_sz(ii); // initialise integral calc
	bq_uz_sz = iq_uz_sz / sub_step; // inflow rate from unsaturated zone
	lbq_uz_sz = bq_uz_sz * sinbeta(ii) ;// flow for uz as rate allowing for angle
	
	// compute the inflow limited by saturation
	l_szin = std::min( il_sz_in(ii_id) / sub_step, l_szmax(ii) );
	
	// compute flow for velocity calculation
	lbar = std::min( (2.0/3.0)*(l_szin+lbq_uz_sz) + l_sz(ii)/3.0,
			 l_szmax(ii) );
	
	// compute lambda
	// using  c = l*Delta_x*cos(beta)/m
	// lambda = c*Delta_t/Delta_x
	lambda = lbar*sub_step*cosbeta_m(ii);
	lp = 1.0;
	if( lambda<=0.0 ){
	  lp = 0.0;
	}
	
	// solve for estimate of outflow
	l_sz(ii) = std::min( (lp*l_sz(ii) + lambda*(l_szin + bq_uz_sz)) / (1.0+lambda),
			     l_szmax(ii) );
	
	// update integral of outflow
	il_sz = il_sz + (sub_step/2.0)*l_sz(ii);
	
	// update volumes in hillslope
	s_sz(ii) = s_sz(ii) + il_sz - il_sz_in(ii_id) - iq_uz_sz;
	
	// pass lateral flux downslope
	dir = sz_dir(ii);
	idx = dir("idx");
	frc = dir("frc");
	for(int jj =0; jj < idx.size(); jj++){
	  iidx = idx(jj) - 1;
	  il_sz_in(iidx) = il_sz_in(iidx) + frc(jj)*il_sz;
	}
	
	// step 5 - correct the stores for saturation flows
	// work out flow to surface
	if( s_sz(ii) <= 0.0 ){
	  iq_sz_sf = -s_sz(ii);
	  iq_uz_sf = s_uz(ii);
	}else{
	  iq_sz_sf = 0.0;
	  iq_uz_sf = 0.0;
	}
	// apply to stores
	s_sz(ii) = s_sz(ii) + iq_sz_sf;
	s_uz(ii) = s_uz(ii) - iq_uz_sf;
	s_sf(ii) = s_sf(ii) + iq_rz_sf + iq_sz_sf + iq_uz_sf;
      }
      
      // loop channel elements to  update volume of flow to channel
      for(int ii =0; ii < nch; ii++){
	ii_id = cid(ii) - 1;
	channel_inflow(it,ii) = channel_inflow(it,ii) +
	  il_sf_in(ii_id) + il_sz_in(ii_id);
      }      
      
    }
    
    // finalise the mass errors calc
    if( mass_check ){
      mass_errors(it,1) = sum( (s_sf + s_rz + s_uz - s_sz)*area ); // final states
      mass_errors(it,4) = sum( channel_inflow(it,_) * carea); // channel inflow
      mass_errors(it,5) = mass_errors(it,0) - mass_errors(it,1) +
	mass_errors(it,2) - mass_errors(it,3) - mass_errors(it,4);
    }
    
    // convert channel inflow to m3/s
    channel_inflow(it,_) = channel_inflow(it,_) * carea / step;
    
    // store states if requested
    if( keep_states(it) ){
      state_record(it) = clone(df_state);
    }
  }
}

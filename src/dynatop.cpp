#include "Rcpp.h"
using namespace Rcpp;

double fode(double a,double b,double x0,double t){
  double ebt = exp(-b*t);
  double kappa = t;
  if( b > 0.0 ){
    kappa = (1-ebt)/b;
  }
  double x = x0*ebt + a*kappa;
  return x;
}

// [[Rcpp::export]]
void hs_sim_cpp(const NumericMatrix obs, NumericMatrix channel_inflow, NumericMatrix mass_errors,
		List hillslope, const List channel, const List ts,
		const List sz_dir, const List sf_dir){
  
  // seperate out the hillslope
  // states
  NumericVector s_sf = hillslope("s_sf"), s_rz = hillslope("s_rz"),
    s_uz = hillslope("s_uz"), s_sz = hillslope("s_sz"), l_sz = hillslope("l_sz");
  // parameters & properties
  const IntegerVector id = hillslope("id"), pcol = hillslope("precip"),
    ecol = hillslope("pet");
  const NumericVector t_sf = hillslope("t_sf"),q_sfmax = hillslope("q_sfmax"),
    s_rzmax = hillslope("s_rzmax"), t_d = hillslope("t_d"), 
    s_bar = hillslope("s_bar"), m = hillslope("m"), ln_t0 = hillslope("ln_t0"),
    delta_x = hillslope("delta_x");

  // seperate out the channel
  // properties
  const IntegerVector cid = channel("id");
  const NumericVector carea = channel("area");

  // seperate out the time variables
  const int n_sub_step = ts("n_sub_step");
  const double step = ts("step"), sub_step = ts("sub_step");
  
  // dimensions
  const int nit = obs.nrow(), nhsu = id.size(),
    nch = cid.size(), max_id = max(id);

  // Rcout << "Finished Declarations" << std::endl;
  // Rcout << max_id << std::endl;
		     
  // local constants used in the evaluation
  const NumericVector beta = atan( s_bar );
  const NumericVector cosbeta_m = cos(beta) / m;
  const NumericVector sinbeta = sin(beta);
  const NumericVector l_szmax = exp( ln_t0 )* sinbeta / delta_x;
  
  // Rcout << "Computed constants" << std::endl;
  
  // loop time steps
  for(int it =0; it < nit; it++) {

    // initialise channel inflow to zero
    for(int ii =0; ii < nch; ii++){
      channel_inflow(it,ii) = 0.0;
    }

    
    // loop the sub steps
    for(int ss=0; ss < n_sub_step; ss++){
      // declare storages for lateral flows
      NumericVector il_sf_in(max_id), il_sz_in(max_id);
      
      // loop hsus for the surface store
      for(int ii = 0; ii < nhsu; ii++){
	// solve for outflow
	int ii_id = id(ii) - 1;

	// Rcout << "ii is: " << ii << std::endl
	//      << "ii_id is: " << ii_id << std::endl;

	double a = il_sf_in(ii_id) / sub_step;
	double b = sub_step / t_sf(ii);
	double il_sf = fode(a,b,s_sf(ii),sub_step);
	// mass balance for new state
	s_sf(ii) += il_sf_in(ii_id) - il_sf;
	// pass downslope
	List dir = sf_dir(ii);
	IntegerVector idx = dir("idx");
	NumericVector frc = dir("frc");
	int nd = idx.size();
	for(int jj = 0; jj < nd; jj++){
	  int i = idx(jj) - 1;
	  double f = frc(jj);
	  // Rcout << "i is: " << i << std::endl;
	  il_sf_in( i ) += f*il_sf;
	}				 
      }
      // Rcout << "looped surface" << std::endl;
      // loop hsus for the subsurface
      for(int ii = 0; ii < nhsu; ii++){

	// set external inputs
	int pc = pcol(ii), ec = ecol(ii);
	double p = obs(it,pc)/step, e_p = obs(it,ec)/step;
	
	// Step 2: solve the root zone for hillslope elements
	// compute inflow from surface
	double q_sf_rz = s_sf(ii) / sub_step;
	if ( q_sf_rz > q_sfmax(ii) ){
	  q_sf_rz = q_sfmax(ii);
	}
	// correct surface state
	s_sf(ii) -= q_sf_rz*sub_step;
	// solve ODE
	double a = p+q_sf_rz;
	double b = e_p / s_rzmax(ii);
	double ts_rz = fode(a,b,s_rz(ii),sub_step);
	
	       
	// work out actual evapotranspiration by mass balance
	// ie_t <- hillslope$s_rz + p*ts$sub_step + iq_sf_rz - ts_rz
	      
	//       if(mass_check){
	// 	mass_errors[it,'e_t'] =mass_errors[it,'e_t'] +
	// 			       sum(ie_t*hillslope$area)
	// 			       }
	// new storage value
	s_rz(ii) = ts_rz;
	if( ts_rz > s_rzmax(ii) ){ s_rz(ii) = s_rzmax(ii); }
	
        // split root zone flow
	double iq_rz_uz = ts_rz - s_rz(ii); // presume goes to unsaturated zone
	double iq_rz_sf = 0.0;
	if ( s_sz(ii) <= 0.0 ){
	  iq_rz_sf = iq_rz_uz;
	  iq_rz_uz = 0.0;
	}
	                                    
	// Step 3: Unsaturated zone
	// solve ode
	double ts_uz = fode( iq_rz_uz/sub_step,
			     1 / (t_d(ii) * s_sz(ii)),
			     s_uz(ii),sub_step );
	// work out outflow
	double iq_uz_sz = s_uz(ii) + iq_rz_uz - ts_uz;
	s_uz(ii) = ts_uz;

	// Step 4: Solve saturated zone
	// this is omega=theta=1 solution with
	// explicit velocity and average inflow                    
	int ii_id = id(ii) - 1;
	double il_sz = (sub_step/2.0)*l_sz(ii); // initialise integral calc
	double bq_uz_sz = iq_uz_sz / sub_step;
	double lbq_uz_sz = bq_uz_sz * sinbeta(ii) ; // flow for uz as rate allowing for angle
	
	double l_szin = il_sz_in(ii_id) / sub_step;
	if (l_szin > l_szmax(ii) ){ l_szin = l_szmax(ii); }
	
	// compute flow for velocity calculation
	double lbar = (2/3)*(l_szin+lbq_uz_sz) + l_sz(ii)/3.0;

	// compute lambda
	// using  c = l*Delta_x*cos(beta)/m
	// lambda = c*Delta_t/Delta_x
	double lambda = lbar*sub_step*cosbeta_m(ii);
	double lp = 1.0;
	if(lambda==0.0){lp = 0;}

	// solve for estimate of outflow
	double tl_sz = (lp*l_sz(ii) + lambda*(l_szin + bq_uz_sz)) / (1+lambda);
	if( tl_sz > l_szmax(ii) ){
	  tl_sz = l_szmax(ii);
	}
	l_sz(ii) = tl_sz;
	
	// integral of outflow
	il_sz += (sub_step/2.0)*l_sz(ii);
	// update volumes in hillslope
	s_sz(ii) = s_sz(ii) + il_sz - il_sz_in(ii_id) - iq_uz_sz;
	double iq_sz_sf = 0.0;
	if( s_sz(ii) <= 0.0 ){
	  iq_sz_sf = -s_sz(ii);
	  s_sz(ii) = 0.0;
	}
	// pass lateral flux downslope
	List dir = sz_dir(ii);
	IntegerVector idx = dir("idx");
	NumericVector frc = dir("frc");
	int nd = idx.size();
	for(int jj = 0; jj < nd; jj++){
	  int i = idx(jj) - 1 ;
	  double f = frc(jj);
	  il_sz_in( i ) += f*il_sz;
	}
	           
	// step 5 - correct the stores for saturation flows
	double iq_uz_sf = 0.0;
	if( s_sz(ii) <= 0.0 ){
	  iq_uz_sf = s_uz(ii);
	}
	s_uz(ii) -= iq_uz_sf;
	s_sf(ii) += iq_rz_sf + iq_sz_sf + iq_uz_sf;
      }
      
      // loop channel elements to  update volume of flow to channel
      for(int ii =0; ii < nch; ii++){
	int ii_id = cid(ii);
	channel_inflow(it,ii) += il_sf_in(ii_id) + il_sz_in(ii_id);
      }

    }
    // Rcout << "looped subsurface" << std::endl;
    // convert channel inflow to m3/s
    for(int ii =0; ii < nch; ii++){
      channel_inflow(it,ii) = channel_inflow(it,ii) * carea(ii) / step;
    }
    Rcout << channel_inflow(it,0) << std::endl;
  }
}

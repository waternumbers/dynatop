// #include "gperftools/profiler.h"
#include "Rcpp.h"
#include "hsu.h"
#include "csu.h"

using namespace Rcpp;
// Thing to do
// 2. Profile and optimise if possible
// 3. check matches r code on MS sims


// [[Rcpp::export]]
void hs_init_cpp(DataFrame hillslope, const List sqnc, const double q0){
  // Rcout << "Running init" << std::endl;
  
  // unpack the saturated zone sequence
  std::vector<int> sqnc_sz = as<std::vector<int>>(sqnc("sz"));
  // unpack the required elements of the hillslope
  IntegerVector id = hillslope("id_index"), p_idx = hillslope("precip_index"),
    e_idx = hillslope("pet_index");
  NumericVector area = hillslope("area"), delta_x = hillslope("delta_x"),
    s_bar = hillslope("s_bar");
  NumericVector t_sf = hillslope("t_sf_value"),
    q_sfmax = hillslope("q_sfmax_value"), s_rzmax = hillslope("s_rzmax_value"),
    s_rz0 = hillslope("s_rz0_value"),
    t_d = hillslope("t_d_value"), m = hillslope("m_value"),
    ln_t0 = hillslope("ln_t0_value");
  NumericVector s_sf = hillslope("s_sf"), s_rz = hillslope("s_rz"),
    s_uz = hillslope("s_uz"), s_sz = hillslope("s_sz"), l_sz = hillslope("l_sz");
  List sf_dir = hillslope("sf_dir_index"), sz_dir = hillslope("sz_dir_index") ;

  // give a time value - not needed except for initialisation of HSU
  double ts = -999.0;
  // Loop to set up the HSUs
  List sf,sz;
  std::vector<hsu> vhsu;
  for(int ii =0; ii<id.size(); ++ii){
    // compute the surface redistribution
    sf = sf_dir(ii);
    sz = sz_dir(ii);
  
    hsu tmp = hsu( id(ii), p_idx(ii), e_idx(ii),
  		   area(ii), s_bar(ii), delta_x(ii),
  		   as<std::vector<double>>(sf("frc")),
  		   as<std::vector<int>>(sf("idx")),
  		   as<std::vector<double>>(sz("frc")),
  		   as<std::vector<int>>(sz("idx")),
  		   t_sf(ii), q_sfmax(ii), s_rzmax(ii),
  		   s_rz0(ii),
  		   t_d(ii), m(ii), ln_t0(ii),
  		   s_sf(ii), s_rz(ii), s_uz(ii), s_sz(ii),
  		   l_sz(ii),
  		   ts);
    vhsu.push_back(tmp);
  }
  // initialise
  int len_redist_vec = max(id)+1; // using rcpp max
  std::vector<double> l_sz_in(len_redist_vec,0.0);
  for(uint sq = 0; sq < sqnc_sz.size(); ++sq){
    int ii = sqnc_sz[sq];
    //Rcout << sq << " " << ii << std::endl;
    vhsu[ii].initialise(q0,l_sz_in);
  }
    
  // I would have though that altering the numericvectors chnaged the data frame but...
  hillslope("s_sf") = s_sf;
  hillslope("s_rz") = s_rz;
  hillslope("s_uz") = s_uz;
  hillslope("s_sz") = s_sz;
  hillslope("l_sz") = l_sz;
}


// [[Rcpp::export]]
void hs_sim_cpp(DataFrame hillslope, const DataFrame channel,
		const List sqnc, const NumericMatrix obs, const List ts,
		NumericMatrix channel_inflow,
		bool mass_check, NumericMatrix mass_errors,
		LogicalVector keep_states, List state_record){
  
  //ProfilerStart("dynaprof.log");
  
  // unpack the sequences
  std::vector<int> sqnc_sz = as<std::vector<int>>(sqnc("sz"));
  std::vector<int> sqnc_sf = as<std::vector<int>>(sqnc("sz"));
  
  // unpack the required elements of the hillslope
  IntegerVector id = hillslope("id_index"),p_idx = hillslope("precip_index"),
    e_idx = hillslope("pet_index");
  NumericVector area = hillslope("area"), delta_x = hillslope("delta_x"),
    s_bar = hillslope("s_bar");
  NumericVector t_sf = hillslope("t_sf_value"),
    q_sfmax = hillslope("q_sfmax_value"), s_rzmax = hillslope("s_rzmax_value"),
    s_rz0 = hillslope("s_rz0_value"),
    t_d = hillslope("t_d_value"), m = hillslope("m_value"),
    ln_t0 = hillslope("ln_t0_value");
  NumericVector s_sf = hillslope("s_sf"), s_rz = hillslope("s_rz"),
    s_uz = hillslope("s_uz"), s_sz = hillslope("s_sz"), l_sz = hillslope("l_sz");
  List sf_dir = hillslope("sf_dir_index"), sz_dir = hillslope("sz_dir_index") ;

  // unpack time values
  int n_sub_step = ts("n_sub_step");
  double step = ts("step"), sub_step = ts("sub_step");

  // unpack the channel information
  IntegerVector channel_id = channel("id_index"), channel_p_idx = channel("precip_index");
  NumericVector channel_area = channel("area");

  // find out highest id - presumed a hillslope
    
  // Loop to set up the hillslope HSUs
  std::vector<hsu> vhsu;
  List sf,sz;
  for(int ii =0; ii<id.size(); ++ii){
    // compute the surface redistribution
    sf = sf_dir(ii);
    sz = sz_dir(ii);
    
    hsu tmp = hsu( id(ii), p_idx(ii), e_idx(ii),
		   area(ii), s_bar(ii), delta_x(ii),
		   as<std::vector<double>>(sf("frc")),
		   as<std::vector<int>>(sf("idx")),
		   as<std::vector<double>>(sz("frc")),
		   as<std::vector<int>>(sz("idx")),
		   t_sf(ii), q_sfmax(ii), s_rzmax(ii),
		   s_rz0(ii),
		   t_d(ii), m(ii), ln_t0(ii),
		   s_sf(ii), s_rz(ii), s_uz(ii), s_sz(ii),
		   l_sz(ii),
		   sub_step);
    vhsu.push_back(tmp);
  }

  // loop to form channel HSUs
  std::vector<csu> vcsu;
  for(int ii =0; ii<channel_id.size(); ++ii){
    csu tmp = csu(channel_id(ii),
		  channel_p_idx(ii),
		  channel_area(ii), sub_step);
    vcsu.push_back(tmp);
  }
  
  // declare data frame of states
  // these are pointers so update and need to be cloned to store
  DataFrame df_state = DataFrame::create( Named("id") = id,
					  Named("s_sf") = s_sf,
					  Named("s_rz") = s_rz,
					  Named("s_uz") = s_uz,
					  Named("s_sz") = s_sz,
					  Named("l_sz") = l_sz);
  
  // lateral flux records - initialise to 0.0
  int len_redist_vec = max(id) +1; // using rcpp max
  std::vector<double> il_sf_in(len_redist_vec,0.0), il_sz_in(len_redist_vec,0.0);

  // observation vector
  std::vector<double> obs_vec(obs.ncol(),0.0);

  // channel inflow volume vector
  std::vector<double> iq_ch(channel_id.size(),0.0);

  // mass balance vector
  std::vector<double> mass_bal(6,0.0);

  // Loop time
  NumericVector tmp;
  for(int it =0; it < obs.nrow(); ++it) {

    // check for interupt
    if(it % 100 == 0){ Rcpp::checkUserInterrupt(); }
    
    // set the observed vector as a rate
    tmp = obs(it,_) / step;
    obs_vec = as<std::vector<double>>(tmp);
    
    // initialise the mass balance?
    if( mass_check ){
      std::fill(mass_bal.begin(), mass_bal.end(), 0.0);
      // store initial mass if required
      for(uint sq = 0; sq < vhsu.size(); sq++){
	mass_bal[1] += vhsu[sq].get_vol();
      }
    }
    
    // initialise the channel inflow volume to zero
    for(uint sq = 0; sq < vcsu.size(); sq++){
      vcsu[sq].initialise();
    }
    
    // Loop the sub steps
    for(int ss=0; ss < n_sub_step; ss++){

      // clear the lateral fluxes
      std::fill(il_sf_in.begin(), il_sf_in.end(), 0.0);
      std::fill(il_sz_in.begin(), il_sz_in.end(), 0.0);

      // single loop
      for(uint sq = 0; sq < sqnc_sz.size(); sq++){
      	int ii = sqnc_sz[sq];
	double ipa(0.0), iea(0.0);
	vhsu[ii].evolve(obs_vec,il_sf_in,il_sz_in,ipa,iea);
	
      	// get precip and ie_t values for mass balance
      	if(mass_check){
      	  mass_bal[1] += ipa;
      	  mass_bal[2] += iea;
      	}
      }

      // Add to the channel volume
      for(uint sq = 0; sq < vcsu.size(); sq++){
	double ipa(0.0);
	vcsu[sq].add_inflow(obs_vec,il_sf_in,il_sz_in,ipa);
	if(mass_check){
	  mass_bal[1] += ipa;
	}
      }
      
      
    }
    
    // compute the channel inflow
    for(uint sq = 0; sq < vcsu.size(); sq++){
      channel_inflow(it,sq) = vcsu[sq].get_inflow();
    }

    // work out final mass in hsus
    if( mass_check ){
      for(uint sq = 0; sq < vhsu.size(); sq++){
	mass_bal[3] += vhsu[sq].get_vol();
      }
      for(uint sq = 0; sq < vcsu.size(); sq++){
	mass_bal[4] += vcsu[sq].get_vol();
      }
      mass_bal[5] = mass_bal[0] + mass_bal[1] - mass_bal[2] - mass_bal[3] - mass_bal[4];
      for(uint sq = 0; sq < mass_bal.size(); sq++){
	mass_errors(it,sq) = mass_bal[sq];
      }
    }

    // store states if requested
    if( keep_states(it) ){
      state_record(it) = clone(df_state);
    }
  }    
  //ProfilerStop();
} 

#include "gperftools/profiler.h"
#include "Rcpp.h"
using namespace Rcpp;

// Thing to do
// 1. Alter indexing of redistribution so starts at 0
// 2. Profile and optimise if possible
// 3. check matches r code

// hillslope class
// note the code takes pointers to the inputs used
class hillslope_hsu {
private:
  // properties of the HSU - passed in as references
  const int id, p_idx, e_idx;
  const double area, s_bar, delta_x;

  // redistribution variables - passed in as references
  const std::vector<double> sf_frc;
  const std::vector<int> sf_idx;
  const std::vector<double> sz_frc;
  const std::vector<int> sz_idx;
  
  // parameter values - passed in as references
  const double t_sf, q_sfmax, s_rzmax, s_rz0, t_d, m, ln_t0;

  // states - passed as references
  double &s_sf, &s_rz, &s_uz, &s_sz, &l_sz;
  
  // time step - passed in as reference
  const double timestep;

  // local variables computed from properties
  double l_szmax, sinbeta, cosbeta_m, iq_sfmax;
  
  // flux vectors
  double iq_sf_rz, iq_rz_uz, iq_uz_sz;

public:
  // constructor
  hillslope_hsu(const int& id_, const int& p_idx_, const int& e_idx_,
		const double& area_ , const double& s_bar_ , const double& delta_x_ ,
    		const std::vector<double>& sf_frc_, const std::vector<int>& sf_idx_,
		const std::vector<double>& sz_frc_, const std::vector<int>& sz_idx_,
		const double& t_sf_, const double& q_sfmax_, const double& s_rzmax_,
		const double s_rz0_, const double& t_d_,
		const double& m_, const double& ln_t0_,
		double& s_sf_, double& s_rz_, double& s_uz_,
		double& s_sz_, double& l_sz_,
		const double& timestep_):
    // set properties
    id{id_}, p_idx{p_idx_}, e_idx{e_idx_},
    area{area_}, s_bar{s_bar_}, delta_x{delta_x_},
    sf_frc{sf_frc_}, sf_idx{sf_idx_},
    sz_frc{sz_frc_}, sz_idx{sz_idx_},
    t_sf{t_sf_}, q_sfmax{q_sfmax_},
    s_rzmax{s_rzmax_}, s_rz0{s_rz0_},
    t_d{t_d_}, m{m_}, ln_t0{ln_t0_},
    s_sf{s_sf_},s_rz{s_rz_},s_uz{s_uz_},s_sz{s_sz_},l_sz{l_sz_},
    timestep{timestep_}
  {
    // initialise local (non)constants
    // local constants used in the evaluation
    double beta = std::atan(s_bar);
    sinbeta = std::sin(beta);
    cosbeta_m = std::cos(beta) / m;
    l_szmax = std::exp(ln_t0) * sinbeta / delta_x;
    iq_sfmax = (q_sfmax)*(timestep);
  }

  int get_id(){
    return id;
  }

  double get_vol(){
    return area * (s_sf + s_rz + s_uz - s_sz);
  }

  void fode(double& x, double& a, double& b){
    // Rcout << "a is: " << a << std::endl;
    // Rcout << "b is: " << b << std::endl;
    // Rcout << "x0 is: " << x0 << std::endl;
    // Rcout << "t is: " << t << std::endl;
    double ebt = exp(-b*(timestep));
    double kappa = (1-ebt)/b;
    if( b == 0.0 ){
      kappa = timestep;
    }
    // Rcout << "b is: " << b << std::endl;
    // Rcout << "ebt is: " << ebt << std::endl;
    // Rcout << "kappa is: " << kappa << std::endl;
    x = x*ebt + a*kappa;
  }

  void initialise(const double& q_uz_sz, std::vector<double>& il_sz_rec){
    // set surface to 0
    s_sf = 0.0;
    // root zone based on fractions filled
    s_rz = s_rzmax*s_rz0;
    // initial flow to saturated zone
    // solve steady state saturated zone
    double il_sf_in = il_sz_rec[id] / area;
    l_sz = std::min( l_szmax, (il_sf_in/timestep) + q_uz_sz );
    double il_sz = l_sz*area;
    // Pass downslope
    // check flow vector
    if (!sz_idx.empty()) {
      // pass flow on
      for(uint jj =0; jj < sz_idx.size(); ++jj){
    	int ii = sz_idx[jj];
    	il_sz_rec[ii] +=  sz_frc[jj]*area*il_sz;
      }
    }
    // compute deficit based on flow
    s_sz = std::max( 0.0, ( std::log(l_szmax) - std::log(l_sz))/cosbeta_m );
    // solve for unsaturated zone
    s_uz = t_d*q_uz_sz*s_sz;
  }
    
  void sf(std::vector<double>& il_sf_rec){
    //std::cout << id_idx << std::endl;
    // get inflow per unit area
    double il_sf_in = il_sf_rec[id] / area;
    // initialise outflow
    double il_sf = s_sf + il_sf_in;
    // solve for storage
    double a = il_sf_in/timestep;
    double b = timestep/t_sf;
    fode(s_sf, a,b); //il_sf_in/timestep, timestep/t_sf);
    // finalise lateral outflow
    il_sf -= s_sf;
    // check flow vector
    if (!sf_idx.empty()) {
      // pass flow on
      for(uint jj =0; jj < sf_idx.size(); ++jj){
    	int ii = sf_idx[jj];
    	il_sf_rec[ii] +=  sf_frc[jj]*area*il_sf;
      }
    }
    // solve for verticle flow
    iq_sf_rz = std::min(s_sf,iq_sfmax);
    s_sf -= iq_sf_rz;
  }
  
  void rz(std::vector<double>& obs, double& ipa, double& iea){
    //Rcout << "Called rz" <<std::endl;
    double ip = obs[p_idx]*timestep ;
    ipa = ip*area;
    
    // initialise estimate of e_t
    double ie =  s_rz + ip + iq_sf_rz ;
    // update storage estimate
    double a = (ip+iq_sf_rz)/timestep;
    double b = obs[e_idx] / s_rzmax;
    fode(s_rz, a, b);
    // finalise estimate of evap-transpiration
    ie -= s_rz;
    iea = ie*area;
    
    // work out flow - presume goes to unsat
    iq_rz_uz = std::max( s_rz - s_rzmax, 0.0);
    // correct storage
    s_rz -= iq_rz_uz;
    // if saturated then flow to surface
    if( s_sz <= 0.0 ){
      s_sf += iq_rz_uz; // send water to surface in effect iq_rz_sf
      iq_rz_uz = 0.0;
    }
  }
  
  void uz(){
    //Rcout << "Called uz" <<std::endl;
    iq_uz_sz = s_uz + iq_rz_uz;
    // solve for new state
    double a = iq_rz_uz/timestep;
    double b = 1.0 / (t_d * s_sz);
    fode( s_uz, a,b); //iq_rz_uz/timestep, 1.0 / (t_d * s_sz));
    // finalise outflow
    iq_uz_sz -= s_uz;
    // iq_uz_sf = 0.0;
  }
  
  void sz(std::vector<double>& il_sz_rec){
    // total inflow
    double il_sz_in = il_sz_rec[id] / area;
    
    double il_sz = l_sz; // initialise integral calc
    
    double q_uz_sz = iq_uz_sz / timestep; // inflow rate from unsaturated zone
    //double lbq_uz_sz = bq_uz_sz * sinbeta ;// flow for uz as rate allowing for angle
	
    // compute the inflow limited by saturation
    double l_sz_in = std::min( il_sz_in / timestep, l_szmax );
	
    // compute flow for velocity calculation
    double lbar = std::min( (2.0/3.0)*(l_sz_in+(q_uz_sz*sinbeta)) + l_sz/3.0,
			    l_szmax );
    
    // compute lambda
    double lambda = lbar*timestep*cosbeta_m;
    double lp = 1.0;
    if( lambda<=0.0 ){ lp = 0.0; }
    
    // solve for estimate of outflow
    l_sz = std::min( (lp*l_sz + lambda*(l_sz_in + q_uz_sz)) / (1.0+lambda),
		     l_szmax );
    
    // update integral of outflow
    il_sz += l_sz;
    il_sz = (timestep/2.0)*il_sz;
    
    // update volumes in hillslope
    s_sz += il_sz - il_sz_in - iq_uz_sz;
	
    // pass lateral flux downslope
    // check flow vector
    if (!sz_idx.empty()) {
      // pass flow on
      for(uint jj =0; jj < sz_idx.size(); ++jj){
    	int ii = sz_idx[jj];
    	il_sz_rec[ii] +=  sz_frc[jj]*area*il_sz;
      }
    }
    
    // pass flow back up vertically
    if( s_sz <= 0.0 ){
      s_sf += s_uz - s_sz; // ass iq_sz_sf and iq_uz_sf
      s_sz = 0.0;
      s_uz = 0.0;
    }

  }
};

// channel_hsu
// really very simple but here to be developed
class channel_hsu {
private:
  // properties of the HSU - passed in as references
  const int id, p_idx; 
  const double area;
  
  // time step - passed in as reference
  const double timestep;

  // local variables computed from properties
  double total_volume;
  int n_add;

public:
  // constructor
  channel_hsu(const int& id_, const int& p_idx_,
		const double& area_, const double& timestep_):
    // set properties
    id{id_}, p_idx{p_idx_},
    area{area_}, timestep{timestep_}  {
      // initialise local (non)constants
      // local constants used in the evaluation
      total_volume = 0.0;
      n_add = 0.0;
    }

  int get_id(){
    return id;
  }

  double get_vol(){
    return total_volume;
  }

  double get_inflow(){
    return total_volume/(n_add*timestep);
  }

  void initialise(){
    total_volume = 0.0;
    n_add = 0;
  }

  void add_inflow(std::vector<double>& obs, std::vector<double>& il_sf_rec,
		  std::vector<double>& il_sz_rec, double& ipa){
    ipa = area*timestep*obs[p_idx];
    total_volume += ipa + il_sf_rec[id] + il_sz_rec[id];
    n_add += 1;
  }
};

  
// [[Rcpp::export]]
void hs_init_cpp(DataFrame hillslope, const List sqnc, const double q0){
  Rcout << "Running init" << std::endl;
  
  // unpack the saturated zone sequence
  std::vector<int> sqnc_sz = as<std::vector<int>>(sqnc("sz"));
  // unpack the required elements of the hillslope
  IntegerVector id = hillslope("id"),p_idx = hillslope("precip_index"),
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
  List sf_dir = hillslope("sf_dir"), sz_dir = hillslope("sz_dir") ;
  // give a time value - not needed except for initialisation of HSU
  double ts = -999.0;
  // Loop to set up the HSUs
  List sf,sz;
  int max_id = max(id); // using rcpp max
  std::vector<hillslope_hsu> hsu;
  for(int ii =0; ii<id.size(); ++ii){
    // compute the surface redistribution
    sf = sf_dir(ii);
    sz = sz_dir(ii);
    
    hillslope_hsu tmp = hillslope_hsu( id(ii), p_idx(ii), e_idx(ii),
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
    hsu.push_back(tmp);
  }
  // initialise
  std::vector<double> il_sz_in(max_id,0.0);
  for(uint sq = 0; sq < sqnc_sz.size(); ++sq){
    int ii = sqnc_sz[sq];
    hsu[ii].initialise(q0,il_sz_in);
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
  IntegerVector id = hillslope("id"),p_idx = hillslope("precip_index"),
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
  List sf_dir = hillslope("sf_dir"), sz_dir = hillslope("sz_dir") ;

  // unpack time values
  int n_sub_step = ts("n_sub_step");
  double step = ts("step"), sub_step = ts("sub_step");

  // unpack the channel information
  IntegerVector channel_id = channel("id"), channle_p_idx = channel("precip_index");
  NumericVector channel_area = channel("area");

  // find out highest id - presumed a hillslope
  int max_id = max(id); // using rcpp max
  
  // Loop to set up the hillslope HSUs
  std::vector<hillslope_hsu> hsu;
  List sf,sz;
  for(int ii =0; ii<id.size(); ++ii){
    // compute the surface redistribution
    sf = sf_dir(ii);
    sz = sz_dir(ii);
    
    hillslope_hsu tmp = hillslope_hsu( id(ii), p_idx(ii), e_idx(ii),
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
    hsu.push_back(tmp);
  }

  // loop to form channel HSUs
  std::vector<channel_hsu> csu;
  for(int ii =0; ii<channel_id.size(); ++ii){
    channel_hsu tmp = channel_hsu( id(ii), p_idx(ii), area(ii), sub_step);
    csu.push_back(tmp);
  }

  // declare data frame of states
  // these are pointers so update and need to be cloned to store
  DataFrame df_state = DataFrame::create( Named("id") = id,
					  Named("id") = s_sf,
					  Named("id") = s_rz,
					  Named("id") = s_uz,
					  Named("id") = s_sz,
					  Named("id") = l_sz);
  
  // lateral flux records - initialise to 0.0
  std::vector<double> il_sf_in(max_id,0.0), il_sz_in(max_id,0.0);

  // observation vector
  std::vector<double> obs_vec(obs.ncol(),0.0);

  // channel inflow volume vector
  std::vector<double> iq_ch(channel_id.size(),0.0);

  // mass balance vector
  std::vector<double> mass_bal(6,0.0);

  // Loop time
  NumericVector tmp;
  for(int it =0; it < obs.nrow(); ++it) {

    // set the observed vector as a rate
    tmp = obs(it,_) / step;
    obs_vec = as<std::vector<double>>(tmp);
    
    // initialise the mass balance?
    if( mass_check ){
      std::fill(mass_bal.begin(), mass_bal.end(), 0.0);
      // store initial mass if required
      for(uint sq = 0; sq < hsu.size(); sq++){
	mass_bal[1] += hsu[sq].get_vol();
      }
    }
    
    // initialise the channel inflow volume to zero
    for(uint sq = 0; sq < channel_id.size(); sq++){
      csu[sq].initialise();
    }
    
    // Loop the sub steps
    for(int ss=0; ss < n_sub_step; ss++){

      // clear the lateral fluxes
      std::fill(il_sf_in.begin(), il_sf_in.end(), 0.0);
      std::fill(il_sz_in.begin(), il_sz_in.end(), 0.0);

      // surface loop
      for(uint sq = 0; sq < sqnc_sf.size(); sq++){
	int ii = sqnc_sf[sq];
	//Rcout << ii << std::endl;
	hsu[ii].sf(il_sf_in);
      }
    
      // subsurface loop
      for(uint sq = 0; sq < sqnc_sz.size(); sq++){
      	int ii = sqnc_sz[sq];
	double ipa(0.0), iea(0.0);
      	hsu[ii].rz(obs_vec,ipa,iea);
      	hsu[ii].uz();
      	hsu[ii].sz(il_sz_in);
	
      	// get precip and ie_t values for mass balance
      	if(mass_check){
      	  mass_bal[1] += ipa;
      	  mass_bal[2] += iea;
      	}
      }

      // Add to the channel volume
      for(uint sq = 0; sq < channel_id.size(); sq++){
	double ipa(0.0);
	csu[sq].add_inflow(obs_vec,il_sf_in,il_sz_in,ipa);
	if(mass_check){
	  mass_bal[1] += ipa;
	}
      }
      
      
    }
    
    // compute the channel inflow
    for(uint sq = 0; sq < csu.size(); sq++){
      channel_inflow(it,sq) = csu[sq].get_inflow();
    }

    // work out final mass in hsus
    if( mass_check ){
      for(uint sq = 0; sq < hsu.size(); sq++){
	mass_bal[3] += hsu[sq].get_vol();
      }
      for(uint sq = 0; sq < csu.size(); sq++){
	mass_bal[4] += csu[sq].get_vol();
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

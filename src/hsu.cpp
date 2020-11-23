#include "hsu.h"
// hillslope class
// note the code takes references to the inputs used

// constructor
hsu::hsu(const int& id_, const int& p_idx_, const int& e_idx_,
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
  // empty and size flags for indexes
  e_sf = sf_idx.empty();
    n_sf = sf_idx.size();
    e_sz = sz_idx.empty();
    n_sz = sz_idx.size();
}

int hsu::get_id(){
    return id;
  }

double hsu::get_vol(){
  return area * (s_sf + s_rz + s_uz - s_sz);
}

void hsu::fode(double& x, double& a, double& b){
  // x = (x + a*timestep) / (1 + b*timestep);
  if( b== 0.0 ){
    x = x + a*timestep;
  }else{
    double ebt = std::exp(-b*(timestep));
    x = x*ebt + a*(1-ebt)/b;
  }
}

void hsu::initialise(const double& q_uz_sz, std::vector<double>& l_sz_rec){
  // set surface to 0
  s_sf = 0.0;
  // root zone based on fractions filled
  s_rz = std::min(1.0,s_rz0)*s_rzmax;
  // initial flow to saturated zone
  // solve steady state saturated zone
  double l_sf_in = std::min( l_szmax,l_sz_rec[id] / area);
  l_sz = std::min( l_szmax, l_sf_in + q_uz_sz );
  // Pass downslope
  // check flow vector
  if (!e_sz) {
    // pass flow on
    for(uint jj =0; jj < n_sz; ++jj){
      int ii = sz_idx[jj];
      l_sz_rec[ii] +=  sz_frc[jj]*area*l_sz;
    }
  }
  // compute deficit based on flow
  s_sz = std::max( 0.0, ( std::log(l_szmax) - std::log(l_sz))/cosbeta_m );
  // solve for unsaturated zone
  s_uz = t_d*q_uz_sz*s_sz;
}

void hsu::sf(std::vector<double>& il_sf_rec){
  // get inflow per unit area
  double il_sf_in = il_sf_rec[id] / area;
  // initialise outflow
  double il_sf = s_sf + il_sf_in;
  // solve for storage
  double a = il_sf_in/timestep;
  double b = timestep/t_sf;
  fode(s_sf, a,b);
  // finalise lateral outflow
  il_sf -= s_sf;
  // check flow vector
  if (!e_sf) {
    int ii;
    //pass flow on
    for(uint jj =0; jj < n_sf; ++jj){
      ii = sf_idx[jj];
      il_sf_rec[ii] +=  sf_frc[jj]*area*il_sf;
    }
  }
  // solve for verticle flow
  iq_sf_rz = std::min(s_sf,iq_sfmax);
  s_sf -= iq_sf_rz;
}

void hsu::rz(std::vector<double>& obs, double& ipa, double& iea){
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
  // if( s_sz <= 0.0 ){
  //   s_sf += iq_rz_uz; // send water to surface in effect iq_rz_sf
  //   iq_rz_uz = 0.0;
  // }
}

void hsu::uz(){
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

void hsu::sz(std::vector<double>& il_sz_rec){
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
  if (!e_sz) {
    double fa = area*il_sz;

    // pass flow on
    int ii;
    double q;
    for(uint jj =0; jj < n_sz; ++jj){
      //ii = sz_idx[jj];
      //q = il_sz_rec[jj];
      //q = 
      il_sz_rec[sz_idx[jj]] +=  sz_frc[jj]*fa; //area*il_sz;
      //il_sz_rec[sz_idx[jj]] = il_sz_rec[sz_idx[jj]] + sz_frc[jj]*fa; 
    }
  }
  
  // pass flow back up vertically
  if( s_sz <= s_uz ){
    s_sf += s_uz - s_sz; // as iq_sz_sf and iq_uz_sf
    s_sz = 0.0;
    s_uz = 0.0;
  }
  
}



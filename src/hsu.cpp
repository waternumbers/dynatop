#include "hsu.h"
#include <iostream>
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
  width = area/delta_x;
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

// compute storage from l_sz
double hsu::fsz(double& x){
  return (std::log(width*sinbeta) + ln_t0 - std::log(x))/cosbeta_m ;
}

// compute kinematic velocity form l_sz
double hsu::fc(double& x){
  return cosbeta_m*x/width;
}

void hsu::initialise(const double& q_uz_sz, std::vector<double>& l_sz_rec){
  // set surface to 0
  s_sf = 0.0;
  // root zone based on fractions filled
  s_rz = std::min(1.0,s_rz0)*s_rzmax;
  // initial flow to saturated zone
  // solve steady state saturated zone
  double l_sz_in = std::min( l_szmax,l_sz_rec[id] / area);
  l_sz = std::min( l_szmax, l_sz_in + q_uz_sz );
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
  s_sz = std::max( 0.0, fsz(l_sz)); //( std::log(l_szmax) - std::log(l_sz))/cosbeta_m );
  // solve for unsaturated zone
  s_uz = t_d*q_uz_sz*s_sz;
  if(s_uz > s_sz){ s_uz = s_sz; }
}

void hsu::evolve( std::vector<double>& obs,
		  std::vector<double>& il_sf_rec,
		  std::vector<double>& il_sz_rec,
		  double& ipa, double& iea){
  // get inflows
  double il_sf_in = il_sf_rec[id] / area;
  double il_sz_in = il_sf_rec[id] / area;
  double ip = obs[p_idx] ; //*timestep ;
  double iep = obs[e_idx] ; //*timestep ;
  ipa = ip*area;

  // solve down
  // surface
  // these are constant and could move to initialisation
  double gamma_sf = std::exp(-timestep/t_sf);
  double lambda_sf = (t_sf/timestep)*(1-gamma_sf);
  double iq_sf_sz = ((gamma_sf/lambda_sf)*s_sf)+il_sf_in ;
  if( iq_sf_sz > (timestep*q_sfmax) ){
    iq_sf_sz = timestep*q_sfmax;
  }
  //double iq_sf_sz = std::min(timestep*q_sfmax, ((gamma_sf/lambda_sf)*s_sf)+il_sf_in)

  // root zone
  double gamma_rz = 1.0;
  double lambda_rz = 1.0;
  if( iep > 0.0 ){
    gamma_rz = std::exp(-iep/s_rzmax) ;
    lambda_rz = (s_rzmax/iep)*(1-gamma_rz) ;
  }
  double trz = gamma_rz * s_rz + lambda_rz*(ip+iq_sf_rz) ;
  double iq_rz_uz = std::max(0.0, (trz-s_rzmax)/lambda_rz) ;

  // sat parameters
  double chat = (cosbeta_m*delta_x/area)*l_sz ; //velocity
  chat = chat*timestep/delta_x ;
  double gamma_sz = 1.0;
  double lambda_sz = 1.0;
  if( chat > 0.0 ){
    gamma_sz = std::exp(-chat);
    lambda_sz = (1-gamma_sz)/chat ;
  }
  
  // unsaturated zone
  double iq_uz_sz = s_uz + iq_rz_uz ;
  if( l_sz < l_szmax ){
    // may be free draining
    double gamma_uz = exp(-timestep/(t_d*s_sz));
    double lambda_uz = ( (t_d*s_sz)/timestep )*(1-gamma_uz);
    double tuz = gamma_uz*s_uz + lambda_uz*iq_rz_uz ;
    double tiq = s_uz + iq_rz_uz - tuz ;
    double tlz = gamma_sz*l_sz + lambda_sz*(il_sz_in+tiq) ;
    double tsz = fsz(tlz);
    if( tsz > tuz ){
      iq_uz_sz = tiq;
    }
  }

  // see if it saturates
  double iq_sat = ((l_szmax - gamma_sz*l_sz)/lambda_sz) - il_sz_in ;
  if( iq_uz_sz > iq_sat ){
    // adjust fluxes
    iq_uz_sz = iq_sat;
    iq_rz_uz = iq_uz_sz - s_uz;
    iq_sf_rz = std::min(iq_sf_rz,
		   ((s_rzmax - gamma_rz*s_rz)/lambda_rz) - iep + iq_rz_uz) ;
  }

  //std::cout << iq_sf_rz << " " << iq_rz_uz << " " << iq_uz_sz << std::endl;
  
  // solve with final flows
  // saturated zone
  l_sz = gamma_sz*l_sz + lambda_sz*(il_sz_in+iq_uz_sz) ;
  s_sz = fsz(l_sz) ;
  double il_sz = timestep*l_sz ; // this is a tempory bodge to check compiling
  // unsaturated zone
  s_uz = s_uz + iq_rz_uz - iq_uz_sz;
  // root zone inc iet 
  double iet = s_rz + iep + iq_sf_rz - iq_rz_uz;
  s_rz = gamma_rz*s_rz + lambda_rz*(iep + iq_sf_rz - iq_rz_uz);
  iet = iet - s_rz;
  iea = iet*area;
  // surface inc outflow
  double il_sf = s_sf + il_sf_in - iq_sf_rz;
  s_sf = gamma_sf*s_sf + lambda_sf*(il_sf_in - iq_sf_rz);
  il_sf  = il_sf - s_sf;

  //std::cout << s_sf << " " << s_rz << " " << s_uz << " " << s_sz << " " << l_sz << std::endl;
  
  // redistribute the outflow
  if (!e_sz) {
    double fsz = area*il_sz;
    double fsf = area*il_sf;
    // pass flow on
    for(uint jj =0; jj < n_sz; ++jj){
      il_sz_rec[sz_idx[jj]] +=  sz_frc[jj]*fsz;
      il_sf_rec[sz_idx[jj]] +=  sz_frc[jj]*fsf;
    }
  }
}


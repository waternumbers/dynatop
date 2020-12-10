// a test for passing data to the hillslope class
#include <vector>
#include <cmath>
#include <iostream>
#include <boost/integer/common_factor.hpp>  
#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>

class hsu {
private:

  //local state variables - referenced to above state vector [m]
  double &l_sf, &s_rz, &s_uz, &l_sz;

  // lateral flux values [m3/s]
  double &q_sf_in, &q_sz_in;
  double &p, &ep;

  // //local property variables
  double &w, &Dx, &beta;
  double &t_sf, &k_sf;
  double &s_rzmax;
  double &t_d;
  double &m, &ln_t0;
  const double &Dt;
  
  // // other summaries of properties used internally across multiple functions
  double l_szmax=-99.0, log_l_szmax=-99.0, lambda_szmax=-99.0, lambda_sf=-99.0, cosbeta_m=-99.0;
  double l_sz_in=-99.0, l_sf_in=-99.0;

  // // internal fluxes and maximum values of stores
  double max_uz = 0.0, r_sf_rz=0.0, r_rz_uz=0.0, r_uz_sz=0.0, et=0.0;
  
  // parameters of the root finding
  const boost::uintmax_t maxit = 1000;
  int digits = std::numeric_limits<double>::digits;
  int get_digits = (digits * 3) /4;

  // private function
  double flambda_sz(double& l);
  double fsz(double& l);
  double fopt(double l);
public:
  // constructor
  /* hsu(std::vector<double>& states, */
  /*     std::vector<double>& ext, */
  /*     std::vector<double>& prop, */
  /*     double& timestep_); */
  hsu(double& l_sf_, double& s_rz_, double& s_uz_, double& l_sz_,
      double& q_sf_in_, double& q_sz_in_,double& p_, double& ep_,
      double& w_, double& Dx_, double& beta_,
      double& t_sf_, double& k_sf_,
      double& s_rzmax_,
      double& t_d_,
      double& m_, double& ln_t0_,
      double& timestep_);
  void step();
  std::vector<double> get_flux();
  std::vector<double> get_q();
  void init(double& s_rz_0, double& r_uz_sz_0);
};

#include <vector>
#include <cmath>
// hillslope class
// note the code takes pointers to the inputs used
class hsu {
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
  double l_szmax, sinbeta, cosbeta_m, iq_sfmax, width;
  
  // flux vectors
  double iq_sf_rz, iq_rz_uz, iq_uz_sz;

  // dimensions
  uint n_sf, n_sz;
  bool e_sf, e_sz;

public:
  // constructor
  hsu(const int& id_, const int& p_idx_, const int& e_idx_,
      const double& area_ , const double& s_bar_ , const double& delta_x_ ,
      const std::vector<double>& sf_frc_, const std::vector<int>& sf_idx_,
      const std::vector<double>& sz_frc_, const std::vector<int>& sz_idx_,
      const double& t_sf_, const double& q_sfmax_, const double& s_rzmax_,
      const double s_rz0_, const double& t_d_,
      const double& m_, const double& ln_t0_,
      double& s_sf_, double& s_rz_, double& s_uz_,
      double& s_sz_, double& l_sz_,
      const double& timestep_);
  int get_id();
  double get_vol();
  double fc(double& x);
  double fsz(double& x);
  
  void initialise(const double& q_uz_sz, std::vector<double>& l_sz_rec);  
  void evolve(std::vector<double>& obs,
	      std::vector<double>& il_sf_rec,
	      std::vector<double>& il_sz_rec,
	      double& ipa, double& iea);
};

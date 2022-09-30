#ifndef HRU
#define HRU

#include <vector>
#include <utility>
#include <cmath>
#include <memory>
#include <limits>
#include "Rcpp.h" // this is included just to get warning messages out!
#include "sf.h"
#include "sz.h"

// Class for the Hillslope HRU
class hru {
  std::unique_ptr<sfc> sf;
  std::unique_ptr<szc> sz;

  
  double const width, Dx;
  std::vector<double> const sf_param;
  double const s_rzmax, t_d;
  std::vector<double> const sz_param;
  
  std::vector<int> const precip_lnk_id;
  std::vector<double> const precip_lnk_frc;
  std::vector<int> const pet_lnk_id;
  std::vector<double> const pet_lnk_frc;
  std::vector<int> const sf_lnk_id;
  std::vector<double> const sf_lnk_frc;
  std::vector<int> const sz_lnk_id;
  std::vector<double> const sz_lnk_frc;

  
  // double v_sf_in, v_sz_in;
  void  lateral_redistribution(std::vector<double>&, std::vector<double>&);
  double fsz(double&, double&, double&, double&, double const&);
  double fsf(double&, double&, double&, double&, double const&);
  
public:
  // variables initialised
  int const id;
  double s_sf, s_rz, s_uz, s_sz, q_sf, q_sz;
  double const area;
  

  double q_sf_in, q_sz_in;
  double r_sf_rz, r_rz_uz, r_uz_sz;
  double precip, pet, aet;
  
  // initialisation
  hru(int const,
      std::vector<double>,
      std::vector<double> const,
      int const, std::vector<double> const,
      std::vector<double> const,
      std::vector<double> const,
      int const, std::vector<double> const,
      std::vector<int> const, std::vector<double> const,
      std::vector<int> const, std::vector<double> const,
      std::vector<int> const, std::vector<double> const,
      std::vector<int> const, std::vector<double> const
      );

  void init(std::vector<double>&, std::vector<double>&, double, double, double const&, double const&, int const&);
  void update_met(std::vector<double>&);
  void step(std::vector<double>&, std::vector<double>&, double const&, double const&, int const&, double const&);
};

#endif

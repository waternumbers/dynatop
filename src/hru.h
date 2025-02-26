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

  //std::vector<double> states;
  //std::vector<double> const properties;

  std::unique_ptr<sfc> sf;
  std::vector<double> const rz_param;
  std::vector<double> const uz_param;
  std::unique_ptr<szc> sz;
  
  std::vector<int> const precip_lnk_id;
  std::vector<double> const precip_lnk_frc;
  std::vector<int> const pet_lnk_id;
  std::vector<double> const pet_lnk_frc;
  std::vector<int> const sf_lnk_id;
  std::vector<double> const sf_lnk_frc;
  std::vector<int> const sz_lnk_id;
  std::vector<double> const sz_lnk_frc;

  //double const &s_rzmax = rz_param[0];
  //double const &t_d = uz_param[0];
  
public:
  // variables initialised
  int const id;

  double s_sf, s_rz, s_uz, s_sz;
  double const area;
  
  //double &s_sf = states[0];
  //double &s_rz = states[1];
  //double &s_uz = states[2];
  //double &s_sz = states[3];
  //double const &area{ properties[0] }; // area if the area of the HRU

  double q_sf, q_sz;
  double q_sf_in, q_sz_in;
  double v_sf_rz, v_rz_uz, v_uz_sz;
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

  void init(std::vector<double>&, std::vector<double>&, double, double, int const&);
  void update_met(std::vector<double>&);
  void  lateral_redistribution(std::vector<double>&, std::vector<double>&);
  void step(std::vector<double>&, std::vector<double>&, int const&, double const&);
};

#endif

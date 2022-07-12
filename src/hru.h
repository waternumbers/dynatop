#ifndef HRU
#define HRU

#include <vector>
#include <utility>
#include <cmath>
#include <memory>
#include "Rcpp.h" // this is included just to get warning messages out!
#include "q_sf.h"
#include "q_sz.h"

// Class for the Hillslope HRU
class hru {
  // these are references to values stored with in the class
  double &s_sf, &s_rz, &s_uz, &s_sz; // states
  double &q_sf, &q_sz, &q_sf_in, &q_sz_in; // lateral fluxes
  double &r_sf_sp, &r_sf_rz, &r_sp_rz, &r_rz_uz, &r_uz_sz; // between store fluxes
  double const s_rzmax; // since multiple options for s_rz not implimented
  double const r_uzmax; // since multiple options for s_uz not implimented
  double const area;
  double const &Dt, &ztol, &etol;
  int const &max_it;
  std::unique_ptr<qsf> sf_func;
  std::unique_ptr<qsf> sp_func;
  std::unique_ptr<qsz> sz_func;
public:
  // initialisation
  hru(double& ,double&, double&, double&, 
      double& ,double&, double&, double&, 
      double& ,double&, double&,
      double& ,double&, double&, double&, double& ,double&,
      double const, double const, double const,
      int const, std::vector<double> const,
      int const, std::vector<double> const,
      int const, std::vector<double> const,
      int const, std::vector<double> const,
      int const, std::vector<double> const,
      double const&, double const&, int const&);
  
  std::pair<double, double> courant(double& Dt);
  void init(double& s_rz_0, double& r_uz_sz_0);
  void step();
  double fz(double&);
};

#endif

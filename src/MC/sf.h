#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double Dx;
public:
  // initialisation
  sfc();
  double kappa{-999.9}, eta{-999.9};
  virtual void update(double const&); // multiply value by storage to give outflow
};

// two section constant velocity
class sfc_cnst: public sfc {
protected:
  double v_raf, q_raf, v_sf;
public:
  sfc_cnst(std::vector<double> const&, std::vector<double> const&);
  void update(double const&);
};

// kinematic with RAF
class sfc_kin: public sfc {
  double v_raf, q_raf, rho, width;
public:
  sfc_kin(std::vector<double> const&, std::vector<double> const&);
  void update(double const&);
};

#endif

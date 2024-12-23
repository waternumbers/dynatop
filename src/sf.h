#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double v_1, s_1, v_2, rho;
public:
  // initialisation
  sfc();
  virtual double fk(double const&); // multiply value by storage to give outflow
};

// two section constant velocity
class sfc_cnst: public sfc {
public:
  sfc_cnst(std::vector<double> const&, double const&, double const&);
};

// kinematic with RAF
class sfc_kin: public sfc {
public:
  sfc_kin(std::vector<double> const&, double const&, double const&, double const&);
  double fk(double const&);
};

#endif

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
  virtual double fs(double const&); // multiply value by storage to give outflow
  virtual double fv(double const&); // multiply value by storage to give outflow
};

// two section constant velocity
class sfc_cnst: public sfc {
protected:
  double v_raf, s_raf, v_sf;
public:
  sfc_cnst(std::vector<double> const&, std::vector<double> const&);
  double fs(double const&);
  double fv(double const&);
};

// // kinematic with RAF
// class sfc_kin: public sfc {
//   double v_raf, q_raf, rho, width;
// public:
//   sfc_kin(std::vector<double> const&, std::vector<double> const&);
//   double fs(double const&);
//   double fq(double const&);
// };

#endif

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
  double kappa;
  double eta;
  virtual void update(double const&); // multiply value by storage to give outflow
};

// two section constant velocity
class sfc_cnst: public sfc {
  double v_1, q_1, v_2;
public:
  //double kappa;
  //double eta;
  sfc_cnst(std::vector<double> const&, std::vector<double> const&);
  void update(double const&);
};

// kinematic with RAF
class sfc_kin: public sfc {
  double v_1, q_1, rho, width;
public:
  //double kappa;
  //double eta;
  
  sfc_kin(std::vector<double> const&, std::vector<double> const&);
  void update(double const&);
};

#endif

#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
public:
   // initialisation
  sfc();
  virtual double fq(double&, double&); // outflow given storage and inflow
};

// constant celerity kinematic wave
class sfc_cnstCD: public sfc {
protected:
  double kappa, eta;
public:
  sfc_cnstCD(std::vector<double> const&, double const&);
  double fq(double&, double&);
};

#endif

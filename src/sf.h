#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double kappa0, eta0;
public:
  // initialisation
  sfc();
  virtual void fke(double&, double&, double const&);
};

// constant celerity kinematic wave
class sfc_cnstC: public sfc {
 public:
  sfc_cnstC(std::vector<double> const&, double const& ,double const&);
};
// constant celerity and diffusivity
class sfc_cnstCD: public sfc {
 public:
  sfc_cnstCD(std::vector<double> const&, double const& ,double const&);
};

#endif


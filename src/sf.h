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
  sfc(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual void init(double&, double&, double&, double&);
  virtual double max_down(double&, double&, double const&);
  virtual void solve (double&, double&, double&, double&, double const&);
};

// constant celerity kinematic wave
class sfc_cnst: public sfc {
 protected:
  double kappa, eta;
 public:
  sfc_cnst(std::vector<double> const&, double const& ,double const&);
  void init(double&, double&, double&, double&);
  double max_down(double&, double&, double const&);
  void solve(double&, double&, double&, double&, double const&); 
};

#endif


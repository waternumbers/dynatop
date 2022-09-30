#ifndef SZ
#define SZ

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class for the saturated zone flow
class szc {
protected:
public:
  double q_szmax, sz_max;
  // initialisation
  szc(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual double fq(double&,double&);
};

// exponential
class szc_exp: public szc {
 protected:
  //double const &t0, &m;
  double psi, eta, kappa;  
 public:
  szc_exp(std::vector<double> const, double const ,double const, double const, double const);
  double fq(double&,double&);
};


#endif


#ifndef SZ
#define SZ

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class for the saturated zone flow
class szc {
protected:
  double psi;
public:
  double q_szmax;
  // initialisation
  szc(); 
  virtual double fs(double const&); // compute storage given representative flow
  virtual double fq(double const&); // compute representative flow given storage
};

// exponential
class szc_exp: public szc {
 protected:
 public:
  szc_exp(std::vector<double> const&, std::vector<double> const&);
  double fs(double const&);
  double fq(double const&);
};

// bounded exponential
class szc_bexp: public szc {
 protected:
  double kappa, omega;
 public:
  szc_bexp(std::vector<double> const&, std::vector<double> const&);
  double fs(double const&);
  double fq(double const&);
};

// constant velocity
class szc_cnst: public szc {
protected:
  double kappa, omega;
public:
  szc_cnst(std::vector<double> const&, std::vector<double> const&);
  double fs(double const&);
  double fq(double const&);
};

// double exponential
class szc_dexp: public szc {
protected:
  double kappa, omega;
public:
  szc_dexp(std::vector<double> const&, std::vector<double> const&);
  double fs(double const&);
  double fq(double const&);
};


#endif


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
  // initialisation
  szc();
  double q_szmax;
  double eta{0.5}, kappa{-999.9}, h{-999.9};
  virtual void update(double const&); // update
};

// exponential
class szc_exp: public szc {
 protected:
  double psi, Dx, width;
 public:
  szc_exp(std::vector<double> const&, std::vector<double> const&);
  void update(double const&); // update
};

// // bounded exponential
// class szc_bexp: public szc {
//  protected:
//  public:
//   szc_bexp(std::vector<double> const&, std::vector<double> const&);
//   double fts(double const&);
//   double ftq(double const&);
// };

// // constant velocity
// class szc_cnst: public szc {
// protected:
// public:
//   szc_cnst(std::vector<double> const&, std::vector<double> const&);
//   double ftq(double const&);
//   double fts(double const&);
// };

// double exponential
class szc_dexp: public szc {
protected:
  double psi, psi2, omega, Dx, width;
public:
  szc_dexp(std::vector<double> const&, std::vector<double> const&);
  void update(double const&); // update
};

#endif


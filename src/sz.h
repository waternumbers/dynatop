#ifndef SZ
#define SZ

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class for the saturated zone flow
class szc {
protected:
  double psi, omega, kappa, q_szmax;
public:
  // initialisation
  szc(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual double fq(double const&, double const&); // compute outflow given storage and lateral inflow
  virtual double fs(double const&, double const&); // compute storage given lateral inflow and outflow
  virtual double ftq(double const&); // compute transmissivity flow given storage
  virtual double fts(double const&); // compute transmissivity storage given transmissivity flow
  //virtual void update(double&, double&, double&, double const&, double const&, int const&); // update
};

// exponential
class szc_exp: public szc {
 protected:
 public:
  szc_exp(std::vector<double> const&, std::vector<double> const&);
  double fts(double const&);
  double ftq(double const&);
};

// bounded exponential
class szc_bexp: public szc {
 protected:
 public:
  szc_bexp(std::vector<double> const&, std::vector<double> const&);
  double fts(double const&);
  double ftq(double const&);
};

// constant velocity
class szc_cnst: public szc {
protected:
public:
  szc_cnst(std::vector<double> const&, std::vector<double> const&);
  double ftq(double const&);
  double fts(double const&);
};

// double exponential
class szc_dexp: public szc {
protected:
public:
  szc_dexp(std::vector<double> const&, std::vector<double> const&);
  double ftq(double const&);
  double fts(double const&);
};




#endif


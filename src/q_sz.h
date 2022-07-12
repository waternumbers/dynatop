#ifndef QSZ
#define QSZ

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class for the saturated zone flow
class qsz {
 public:
  // initialisation
  qsz(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual double fq(double&);
  virtual void init(double&, double&, double&, double&);
};

// exponential
class qsz_exp: public qsz {
 protected:
  double qmax, kappa;
 public:
  qsz_exp(std::vector<double> const&, double const& ,double const&, double const&);
  double fq(double&);
  void init(double&, double&, double&, double&);
};

// bounded exponential
class qsz_bexp: public qsz {
 protected:
  double qmax, kappa, omega, D;
 public:
  qsz_bexp(std::vector<double> const&, double const& ,double const&, double const&);
  double fq(double&);
  void init(double&, double&, double&, double&);
};

// double exponential
class qsz_dexp: public qsz {
 protected:
  double qmax, omega, kappa1, kappa2;
 public:
  qsz_dexp(std::vector<double> const&, double const& ,double const&, double const&);
  double fq(double&);
  void init(double&, double&, double&, double&);
};

// constant velocity
class qsz_cnst: public qsz {
 protected:
  double qmax, kappa, D;
 public:
  qsz_cnst(std::vector<double> const&, double const& ,double const&, double const&);
  double fq(double&);
  void init(double&, double&, double&, double&);
};

#endif


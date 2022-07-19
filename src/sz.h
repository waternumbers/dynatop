#ifndef SZ
#define SZ

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class for the saturated zone flow
class szc {
 public:
  // initialisation
  szc(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual double fD();
  virtual double fq_szmax();
  virtual double fc(double&);
  virtual void init(double&, double&, double&, double&, double&);
};

// bounded exponential
class szc_bexp: public szc {
 protected:
  double qmax, kappa, omega, D;
 public:
  szc_bexp(std::vector<double> const&, double const& ,double const&, double const&);
  double fD();
  double fq_szmax();
  double fc(double&);
  void init(double&, double&, double&, double&, double&);
};

/* // double exponential */
/* class qsz_dexp: public qsz { */
/*  protected: */
/*   double qmax, omega, kappa1, kappa2; */
/*  public: */
/*   qsz_dexp(std::vector<double> const&, double const& ,double const&, double const&); */
/*   double fq(double&); */
/*   void init(double&, double&, double&, double&); */
/* }; */

// constant velocity
class szc_cnst: public szc {
 protected:
  double qmax, celerity, D;
 public:
  szc_cnst(std::vector<double> const&, double const& ,double const&, double const&);
  double fD();
  double fc(double&);
  double fq_szmax();
  void init(double&, double&, double&, double&);
};

#endif


#ifndef SZ
#define SZ

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class for the saturated zone flow
class szc {
protected:
  double Dx;
public:
  double q_szmax, D;
  // initialisation
  szc(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual double fc(double const&);
  virtual void fke(double&, double&, double const&);
};


// bounded exponential
class szc_bexp: public szc {
 protected:
  double psi, omega;
 public:
  szc_bexp(std::vector<double> const&, double const& ,double const&, double const&);
  double fc(double const&);
  void fke(double&, double&, double const&);
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
  double kappa0, eta0, celerity;
 public:
  szc_cnst(std::vector<double> const&, double const& , double const&);
  double fc(double const&);
  void fke(double&, double&, double const&);
};

#endif


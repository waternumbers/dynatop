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
  double q_szmax;
   // initialisation
  szc(); //std::vector<double> const&, double const& ,double const&, double const&);
  virtual double fq(double&,double&);
  virtual double fa(double&);
};

// exponential
class szc_exp: public szc {
 protected:
  //double const &t0, &m;
  double Dx; 
  double psi;
 public:
  szc_exp(std::vector<double> const, double const ,double const, double const);
  double fq(double&,double&);
  double fa(double&);
};


#endif


#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double Dx, eta;
public:
  // initialisation
  sfc();
  virtual double fq(double&, double&);
  virtual double finit(double&, double&);
};

// constant celerity kinematic wave
class sfc_cnstC: public sfc {
protected:
  double celerity;
public:
  sfc_cnstC(std::vector<double> const&, double const& ,double const&);
  double fq(double&, double&);
  double finit(double&, double&);
};


// constant celerity and diffusivity
class sfc_cnstCD: public sfc {
protected:
  double celerity;
public:
  sfc_cnstCD(std::vector<double> const&, double const& ,double const&);
  double fq(double&, double&);
  double finit(double&, double&);
};

#endif


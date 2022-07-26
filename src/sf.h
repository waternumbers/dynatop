#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
public:
  double eta;
  // initialisation
  sfc();
  virtual double fq(double&, double&); // outflow given storage and upslope crosssectional area
  virtual double fa(double&); // x-sec are given discharge
};

// constant celerity kinematic wave
class sfc_cnstC: public sfc {
protected:
  double const Dx, celerity;
public:
  sfc_cnstC(std::vector<double> const&, double const&);
  double fq(double&, double&); // outflow given storage and upslope crosssectional area
  double fa(double&); // x-sec are given discharge
};


// constant celerity and diffusivity
class sfc_cnstCD: public sfc {
protected:
  double const  &Dx, celerity;
public:
  sfc_cnstCD(std::vector<double> const&, double const&);
  double fq(double&, double&); // outflow given storage and upslope crosssectional area
  double fa(double&); // x-sec are given discharge
};

#endif


#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
public:
   // initialisation
  sfc();
  virtual double fq(double&, double&); // outflow given storage and inflow
  virtual double fs(double&, double&); // storage given outflow and inflow
};

// constant celerity & diffusivity with RAF
class sfc_cnst: public sfc {
protected:
  double kappa, eta, s_raf, t_raf, q_raf;
  
public:
  sfc_cnst(std::vector<double> const&, std::vector<double> const&);
  double fq(double&, double&);
  double fs(double&, double&);
};

// kinematic with RAF
class sfc_kin: public sfc {
protected:
  double kappa, omega, s_raf, t_raf, q_raf;
  
public:
  sfc_kin(std::vector<double> const&, std::vector<double> const&);
  double fq(double&, double&);
  double fs(double&, double&);
};

// compound channel with RAF
class sfc_comp: public sfc {
protected:
  double kappa_1, eta_1, s_1, kappa_2, eta_2, q_1_max;
  
public:
  sfc_comp(std::vector<double> const&, std::vector<double> const&);
  double fq(double&, double&);
  double fs(double&, double&);
};

#endif

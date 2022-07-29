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
  virtual double max_vert(double&, double&, double const&); // max possible down flow from surface given initial state, inflow and timestep
  virtual void solve(double&, double&, double&, double&, double const&); // solve for storage and outflow given inflows and timestep
  virtual void init(double&, double&, double&, double&);
};

// constant celerity kinematic wave
class sfc_cnstC: public sfc {
protected:
  double const Dx, celerity;
  double const eta = 0.5;
public:
  sfc_cnstC(std::vector<double> const&, double const&);
  double max_vert(double&, double&, double const&); // max possible down flow from surface given initial state, inflow and timestep
  void solve(double&, double&, double&, double&, double const&); // solve for storage and outflow given inflows and timestep
  void init(double&, double&, double&, double&);
};


// constant celerity and diffusivity
class sfc_cnstCD: public sfc {
protected:
  double const  Dx, celerity;
  double eta;
public:
  sfc_cnstCD(std::vector<double> const&, double const&);
  double max_vert(double&, double&, double const&); // max possible down flow from surface given initial state, inflow and timestep
  void solve(double&, double&, double&, double&, double const&); // solve for storage and outflow given inflows and timestep
  void init(double&, double&, double&, double&);
};

// constant celerity with raf
class sfc_cnstC_raf: public sfc {
protected:
  double const  Dx, celerity, sraf, traf;
  double eta = 0.5;
public:
  sfc_cnstC_raf(std::vector<double> const&, double const&);
  double max_vert(double&, double&, double const&); // max possible down flow from surface given initial state, inflow and timestep
  void solve(double&, double&, double&, double&, double const&); // solve for storage and outflow given inflows and timestep
  void init(double&, double&, double&, double&);
};

#endif


#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double kappa_1, eta_1, q_1, kappa_2, eta_2;
public:
  double kappa, eta; // these are used
   // initialisation
  sfc();
  virtual void update(double const&);
};

// constant celerity & diffusivity with RAF
class sfc_cnst: public sfc {
public:
  sfc_cnst(std::vector<double> const&, std::vector<double> const&);
};

// kinematic with RAF
class sfc_kin: public sfc {
public:
  sfc_kin(std::vector<double> const&, std::vector<double> const&);
  void update(double const&);
};

// compound channel with RAF
class sfc_comp: public sfc {  
public:
  sfc_comp(std::vector<double> const&, std::vector<double> const&);
};

// arbitary area flow relationship
class sfc_arb_kin: public sfc {
  std::vector<double> s_val, q_val;
public:
  sfc_arb_kin(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  void update(double const&);
};

// raf with power law
class sfc_power_law: public sfc {
private:
  double Dx;
public:
  sfc_power_law(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  void update(double const&);
};

// MCT
class sfc_mct: public sfc {
private:
  double grd, Dx, n, ca, sa, B0;
public:
  sfc_mct(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  void update(double const&);
};

#endif

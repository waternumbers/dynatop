#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double T_1, T_2, S_1;
public:
   // initialisation
  sfc();
  virtual double fT(double const&);
  virtual double fS(double const&);
};

// constant velocity with RAF
class sfc_cnst: public sfc {
public:
  sfc_cnst(std::vector<double> const&, std::vector<double> const&);
};

// kinematic with RAF
class sfc_kin: public sfc {
  double eta;
public:
  sfc_kin(std::vector<double> const&, std::vector<double> const&);
  double fT(double const&);
  double fS(double const&);
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
  double fT(double const&);
  double fS(double const&);
};

// raf with power law
class sfc_power_law: public sfc {
private:
  double kappa, eta;
public:
  sfc_power_law(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  double fT(double const&);
  double fS(double const&);
};

// MCT
class sfc_mct: public sfc {
private:
  double grd, Dx, n, ca, sa, B0;
public:
  sfc_mct(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  double fT(double const&);
  double fS(double const&);
};

// MCT with double rectangle channel
class sfc_mct_rect: public sfc {
private:
  double Dx, B0, n, ca, sa, A_crit, beta, y_crit;
  double solve_depth(double const&);
public:
  sfc_mct_rect(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  double fT(double const&);
  double fS(double const&);
};

#endif

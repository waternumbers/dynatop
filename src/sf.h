#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double kappa_1, eta_1, s_1, kappa_2, eta_2;
public:
   // initialisation
  sfc();
  virtual double fq(double const&); // outflow given storage
  virtual double fs(double const&, double const&); // storage given outflow and inflow
  virtual void update(double&, double&, double const&, double const&,
		      double const&, double const&, int const&);
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
  double fq(double const&);
  double fs(double const&, double const&);
};

// power law with raf
class sfc_power_law: public sfc {
public:
  sfc_power_law(std::vector<double> const&, std::vector<double> const&);
  double fq(double const&);
  double fs(double const&, double const&);
};
// MCT
class sfc_mct: public sfc {
  double Cs, Ds;
  double grd, n, Dx, ca, sa, B0;
public:
  sfc_mct(std::vector<double> const&, std::vector<double> const&);
  double fq(double const&);
  double fs(double const&, double const&);
  void update(double&, double&, double const&, double const&,
    double const&, double const&, int const&);
  void internal_update(double const&);
};
// MCT
class sfc_arb_mct: public sfc {
  double Cs, Ds;
  double Dx, grd;
  std::vector<double> a_val,q_val,B_val,c_val;
public:
  sfc_arb_mct(std::vector<double> const&, std::vector<double> const&);
  double fq(double const&);
  double fs(double const&, double const&);
  void update(double&, double&, double const&, double const&,
    double const&, double const&, int const&);
  void internal_update(double const&);
};
#endif

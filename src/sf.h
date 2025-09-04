#ifndef SF
#define SF

#include <vector>
#include <utility>
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// generic class
class sfc {
protected:
  double k_1, k_2, S_1;
public:
   // initialisation
  sfc();
  virtual std::pair<double,double> fq(double const&);
  virtual double fs(double const&);
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
  std::pair<double,double> fq(double const&);
  double fs(double const&);
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
  std::pair<double,double> fq(double const&);
  double fs(double const&);
};

// raf with power law
class sfc_power_law: public sfc {
private:
  double kappa, eta;
public:
  sfc_power_law(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  std::pair<double,double> fq(double const&);
  double fs(double const&);
};

// MCT
class sfc_mct: public sfc {
private:
  double grd, Dx, n, ca, sa, B0;
  double Ay(double const&), Py(double const&), Qy(double const&), dQ_dy(double const&);
  //auto Ay, Py, Qy, dQ_dy;
public:
  sfc_mct(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  std::pair<double,double> fq(double const&);
  double fs(double const&);
};

// MCT with double rectangle channel
class sfc_mct_rect: public sfc {
private:
  double Dx, B0, n, ca, sa, s_crit, beta, y_crit;
  double solve_storage(double const&);
  double Ay(double const&), Py(double const&), Qy(double const&), dQ_dy(double const&);
public:
  sfc_mct_rect(std::vector<double> const&, std::vector<double> const&); //, std::vector<double> const&);
  std::pair<double,double> fq(double const&);
  double fs(double const&);
};

#endif

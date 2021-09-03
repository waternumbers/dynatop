#ifndef HILLSLOPE_H
#define HILLSLOPE_H

#include <iostream>
#include <cmath>
#include <boost/math/tools/roots.hpp>
#include <boost/bind/bind.hpp>
#include <boost/bind/placeholders.hpp>

// Class for the Hillslope HRU
class hillslope_hru {
  // these are references to values stored with in the class
  double& s_sf; double& s_rz; double& s_uz; double& s_sz; // states
  double const& s_bar; double const& area; double const& width; // physical properties
  double& q_sf_in; double& q_sf_out; double& q_sz_in; double& q_sz_out;// fluxes in and out of surface and saturated zones
  double& e_a; // actual evapotranspiration as a rate [m/s]
  double const& r_sf_max; double const& c_sf; // surface store parameters
  double const& s_rz_max; // root zone store parameters
  double const& t_d; // unsaturated zone parameters
  double const& ln_t0; double const& m; double const& D; // saturated zone parameters
  int const& type_sz; // type of saturated zone
  double Dx, beta, l_sz_max,cosbeta_m; // values computed during initialisation
  double r_rz_uz, l_sz_in, Dts; // values used during optimisation of sz
public:
  hillslope_hru(double& ,double&, double&, double&,
		double const&, double const&, double const&,
		double&, double&, // surface zone lateral fluxes
		double&, double&, // saturated zone lateral fluxes
		double&, // actual evapotranspiration as a rate [m/s]
		double const&, double const&, // surface store parameters
		double const&, // root zone store parameters
		double const&, // unsaturated zone parameters
		double const&, double const&, double const&, // saturated zone parameters
		int const& // type of sautruated zone
	    );
  std::pair<double, double> courant(double& Dt);
  void init(double& s_rz_0, double& r_uz_sz_0);
  void implicit_step(double& pet, double& precip, double& Dt, int& max_it);
  double fsz(double&, double&); //,double&,double&,double&);
  struct TerminationCondition {
    bool operator()(double min, double max);
  };
  struct Fsz {
    double const& suz; double const& ssz;
    double const& lszmax; double const& cosbetam;
    double const& td; double const d; double const m;
    double const& dt; double const& dx;
    double const& lszin; double const& rrzuz;
    int const& typesz;
  public:
    Fsz(double const& suz_, double const& ssz_,
	double const& lszmax_, double const& cosbetam_,
	double const& td_, double const d_, double const m_,
	double const& dt_, double const& dx_,
	double const& lszin_, double const& rrzuz_,
	int const& typesz);
    double operator()(double const& x);
  };
};

#endif

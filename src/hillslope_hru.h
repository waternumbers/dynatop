#ifndef HILLSLOPE_H
#define HILLSLOPE_H

// #include <iostream> // this is needed for debugging messages
#include <cmath>
#include "Rcpp.h" // this is included just to get warning messages out!

// Class for the Hillslope HRU
class hillslope_hru {
  // these are references to values stored with in the class
  int& id; // hru id
  double& s_sf; double& s_rz; double& s_uz; double& s_sz; // states
  double const& s_bar; double const& area; double const& width; // physical properties
  double& q_sf_in; double& q_sf_out; double& q_sz_in; double& q_sz_out;// fluxes in and out of surface and saturated zones
  double& e_a; // actual evapotranspiration as a rate [m/s]
  double const& r_sf_max; double const& c_sf; // surface store parameters
  double const& s_rz_max; // root zone store parameters
  double const& t_d; // unsaturated zone parameters
  double const& ln_t0; double const& c_sz; double const& m; double const& D; double const& m_2; double const& omega;// saturated zone parameters
  int const& opt; // type of saturated zone
  double Dx, beta, l_sz_max,cosbeta_m,cosbeta_m_2, r_uz_sz_max; // values computed during initialisation
  double v_rz_uz, l_sz_in, Dt_Dx; // values used during optimisation of sz
  // v_rz_uz = Dt*r_rz_uz - used to stop rescaling by Dt in code
  

public:
  hillslope_hru(int&,
		double& ,double&, double&, double&,
		double const&, double const&, double const&,
		double&, double&, // surface zone lateral fluxes
		double&, double&, // saturated zone lateral fluxes
		double&, // actual evapotranspiration as a rate [m/s]
		double const&, double const&, // surface store parameters
		double const&, // root zone store parameters
		double const&, // unsaturated zone parameters
		double const&, double const&, double const&, double const&, double const&, double const&,// saturated zone parameters
		int const& // type of sautruated zone
	    );
  std::pair<double, double> courant(double& Dt);
  void init(double& s_rz_0, double& r_uz_sz_0, double& tol, int& max_it);
  void implicit_step(double& pet, double& precip, double& Dt, double& tol, int& max_it);
  double fsz(double&, double&); //,double&,double&,double&);
  double flz(double&); //,double&,double&,double&);
};

#endif

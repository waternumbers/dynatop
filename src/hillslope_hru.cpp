#include "hillslope_hru.h"

hillslope_hru::hillslope_hru(double& s_sf_,double& s_rz_,double& s_uz_,double& s_sz_,
			     double const& s_bar_, double const& area_, double const& width_,
			     double& q_sf_in_, double& q_sf_out_, // surface zone lateral fluxes
			     double& q_sz_in_, double& q_sz_out_, // saturated zone lateral fluxes
			     double& e_a_, // actual evapotranspiration as a rate [m/s]
			     double const& r_sf_max_, double const& c_sf_, // surface store parameters
			     double const& s_rz_max_, // root zone store parameters
			     double const& t_d_, // unsaturated zone parameters
			     double const& ln_t0_, double const& m_, double const& D_, // saturated zone parameter
			     int const& type_sz_
			     ):
  s_sf(s_sf_), s_rz(s_rz_), s_uz(s_uz_), s_sz(s_sz_),
  s_bar(s_bar_), area(area_), width(width_),
  q_sf_in(q_sf_in_), q_sf_out(q_sf_out_),
  q_sz_in(q_sz_in_), q_sz_out(q_sz_out_),
  e_a(e_a_),
  r_sf_max(r_sf_max_), c_sf(c_sf_),
  s_rz_max(s_rz_max_),
  t_d(t_d_),
  ln_t0(ln_t0_), m(m_), D(D_),
  type_sz(type_sz_)
{
  // compute summary for saturated zone
  beta = std::atan(s_bar_);
  l_sz_max = std::exp(ln_t0_)*std::sin(beta);
  cosbeta_m = std::cos(beta) / m;
  Dx = area_/width;
  r_uz_sz_max = 1/t_d;
}

std::pair<double, double> hillslope_hru::courant(double& Dt){
  std::pair<double, double> cr(-99.0,-99.0);
  cr.first = c_sf*Dt/Dx;
  switch (type_sz) {
  case 1: // exponential transmissivity
    cr.second = std::exp(ln_t0)*std::sin(2.0*beta)/(2.0*m*Dx);
    break;
  case 2: // constant celerity
    cr.second = m*Dt/Dx;
    break;
  case 3: // bounded exponential
    cr.second = std::exp(ln_t0)*std::sin(2.0*beta)/(2.0*m*Dx);
    break;
  }
  
  return cr;
}

void hillslope_hru::init(double& s_rz_0, double& r_uz_sz_0){
  double l_sz_in, l_sz;
  double r_uz_sz;
  l_sz_in = q_sz_in / width; // standardise inflows by width
  s_sf = 0.0; // initialise surface store
  s_rz = s_rz_max*s_rz_0; // initialise root zone store
  l_sz = std::min(l_sz_max, l_sz_in + Dx*r_uz_sz_0); // outflow flux under steady state
  r_uz_sz = (l_sz - l_sz_in)/Dx; // solve to find actual r_uz_sz
  // compute saturated zone storage deficit
  s_sz = std::max(0.0, (std::log(l_sz_max) - std::log(l_sz))/cosbeta_m);
  s_uz = std::min( s_sz, r_uz_sz*t_d*s_sz ); // compute unsaturated zone storage
}

void hillslope_hru::implicit_step(double& pet, double& precip, double& Dt, int& max_it)
{
  double l_sz;
  double r_sf_rz;
 
  std::pair<double, double> opt_res; // solution for output of optimiser
  
  // standardise inflow to the saturated zone by width
  l_sz_in = q_sz_in / width;

  // compute first downward flux estimate from surface ans root zone
  // these are given as \check{r} in documentation	
  r_sf_rz = std::min( r_sf_max , (s_sf + Dt*q_sf_in/area)/Dt );
  r_rz_uz = std::max(0.0 ,
		     (s_rz + Dt*(precip + r_sf_rz - pet) - s_rz_max)/Dt);

  // initialise the bands for the search
  double lwr(0), upr(D);

  // test for saturation
  double x(0.0); // declared here so can be passed into fsz to test for saturation
  if( fsz(x,Dt) >= 0.0 ){
    // then saturated
    l_sz = l_sz_max;
    r_rz_uz = ( s_sz + (Dt/Dx)*(l_sz - l_sz_in) - s_uz )/Dt;
    s_sz = 0.0;
  }else{
    // not saturated test to see if wetting or drying
    if( fsz(s_sz,Dt) >= 0.0 ){
      // then wetting - solution between current s_sz and 0
      lwr = 0.0;
      upr = s_sz;
    }else{
      // drying
      lwr = s_sz;
      // drying - need to work lower depth to bracket
      upr = 2.0*s_sz + 0.01;
      double fupr = fsz(upr,Dt);
      while( (fupr < 0.0) & (upr < D)){
	upr += upr;
	fupr = fsz(upr,Dt);
      }
    }

    int it(0);
    double e;
    while( (it < max_it) and ((upr-lwr)>1e-16) ){
      x = (upr+lwr)/2.0;
      e = fsz(x,Dt);
      if(e<=0){ lwr=x; }
      if(e>=0){ upr=x; }
      it +=1;
    }
    if(it >max_it){
      Rcpp::Rcout << it << std::endl;
    }
    
    s_sz = (upr+lwr)/2.0;
    r_rz_uz = std::min( r_rz_uz, (s_sz + Dt/t_d - s_uz)/Dt );
    l_sz = flz(s_sz);
  }

  // solve unsaturated zone
  s_uz = ( (t_d*s_sz)/(t_d*s_sz + Dt) ) * (s_uz+Dt*r_rz_uz);
  // compute revised r_sf_rz
  r_sf_rz = std::min( r_sf_rz,
		      (s_rz_max - s_rz - Dt*(precip-pet-r_rz_uz))/Dt
		      );
  // solve for root zone
  s_rz = ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) ) * ( s_rz_max / (s_rz_max + Dt*pet) );
  // solve for surface
  s_sf = ( Dx / (Dx + Dt*c_sf) ) * ( s_sf + Dt*( (q_sf_in/area) - r_sf_rz ));
  // compute actual pet loss
  e_a = pet*(s_rz/s_rz_max);

  // compute the outflow
  q_sf_out = width*c_sf*s_sf;
  q_sz_out = width*l_sz;

}

double hillslope_hru::flz(double& x){ // compute saturated zone outflow
  //double l = l_sz_max*std::exp(-x*cosbeta_m);
  double l=-99.0;
  switch (type_sz) {
  case 1: // exponential transmissivity
    l = l_sz_max*std::exp(-x*cosbeta_m);
    break;
  case 2: // constant celerity
    l = std::max( m*(D-x) , 0.0 );
    break;
  case 3: // bounded exponential
    l = l_sz_max * ( std::exp(-x*cosbeta_m) - std::exp( -D *cosbeta_m ) );
    break;
  }
  return l;
}

double hillslope_hru::fsz(double& x, double& Dt){ // compute for saturated zone function
  //double l = l_sz_max*std::exp(-x*cosbeta_m);
  double l = flz(x);
  double r = std::min( r_uz_sz_max , (s_uz + Dt*r_rz_uz)/(t_d*x + Dt) );
  //double r = (s_uz + Dt*r_rz_uz)/(t_d*x + Dt);
  //if( r > r_uz_sz_max ){ r = r_uz_sz_max; }
  
  return x - s_sz - Dt*((l - l_sz_in)/Dx - r);
}


void hillslope_hru::explicit_step(double& pet, double& precip, double& Dt, int& max_it)
{
  double l_sz;
  double r_sf_rz;
 
  std::pair<double, double> opt_res; // solution for output of optimiser
  
  // standardise inflow to the saturated zone by width
  l_sz_in = q_sz_in / width;

  // compute first downward flux estimate from surface ans root zone
  // these are given as \check{r} in documentation	
  r_sf_rz = std::min( r_sf_max , (s_sf + Dt*q_sf_in/area)/Dt );
  r_rz_uz = std::max(0.0 ,
		     (s_rz + Dt*(precip + r_sf_rz - pet) - s_rz_max)/Dt);

  // compute inital estimate fo the outflow
  l_sz = flz(s_sz);

  // parts of the quadratic
  double qc = s_sz - s_sz +Dt*(l_sz-l_sz_in-r_rz_uz);
  double qb = t_d*(s_sz +Dt*(l_sz-l_sz_in)) + Dt;
  double& qa = t_d;

  // test for saturation
  if( qc <= 0.0 ){
    // is saturated so limit r_rz_uz
    r_rz_uz = ( s_sz + (Dt/Dx)*(l_sz - l_sz_in) - s_uz )/Dt;
    s_sz = 0.0;
  }else{
    // test if reaches D
    if( 
  // initialise the bands for the search
  double lwr(0), upr(D);

  // test for saturation
  double x(0.0); // declared here so can be passed into fsz to test for saturation
  if( fsz(x,Dt) >= 0.0 ){
    // then saturated
    l_sz = l_sz_max;
    r_rz_uz = ( s_sz + (Dt/Dx)*(l_sz - l_sz_in) - s_uz )/Dt;
    s_sz = 0.0;
  }else{
    // not saturated test to see if wetting or drying
    if( fsz(s_sz,Dt) >= 0.0 ){
      // then wetting - solution between current s_sz and 0
      lwr = 0.0;
      upr = s_sz;
    }else{
      // drying
      lwr = s_sz;
      // drying - need to work lower depth to bracket
      upr = 2.0*s_sz + 0.01;
      double fupr = fsz(upr,Dt);
      while( (fupr < 0.0) & (upr < D)){
	upr += upr;
	fupr = fsz(upr,Dt);
      }
    }

    int it(0);
    double e;
    while( (it < max_it) and ((upr-lwr)>1e-16) ){
      x = (upr+lwr)/2.0;
      e = fsz(x,Dt);
      if(e<=0){ lwr=x; }
      if(e>=0){ upr=x; }
      it +=1;
    }
    if(it >max_it){
      Rcpp::Rcout << it << std::endl;
    }
    
    s_sz = (upr+lwr)/2.0;
    r_rz_uz = std::min( r_rz_uz, (s_sz + Dt/t_d - s_uz)/Dt );
    l_sz = l_sz_max*std::exp(-s_sz*cosbeta_m);
  }

  // solve unsaturated zone
  s_uz = ( (t_d*s_sz)/(t_d*s_sz + Dt) ) * (s_uz+Dt*r_rz_uz);
  // compute revised r_sf_rz
  r_sf_rz = std::min( r_sf_rz,
		      (s_rz_max - s_rz - Dt*(precip-pet-r_rz_uz))/Dt
		      );
  // solve for root zone
  s_rz = ( s_rz + Dt*(precip + r_sf_rz - r_rz_uz) ) * ( s_rz_max / (s_rz_max + Dt*pet) );
  // solve for surface
  s_sf = ( Dx / (Dx + Dt*c_sf) ) * ( s_sf + Dt*( (q_sf_in/area) - r_sf_rz ));
  // compute actual pet loss
  e_a = pet*(s_rz/s_rz_max);

  // compute the outflow
  q_sf_out = width*c_sf*s_sf;
  q_sz_out = width*l_sz;

}

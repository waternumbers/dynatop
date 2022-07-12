#include "hillslope_hru.h"

hillslope_hru::hillslope_hru(int& id_, double& s_sf_,double& s_rz_,double& s_uz_,double& s_sz_,
			     double const& s_bar_, double const& area_, double const& width_,
			     double& q_sf_in_, double& q_sf_out_, // surface zone lateral fluxes
			     double& q_sz_in_, double& q_sz_out_, // saturated zone lateral fluxes
			     double& e_a_, // actual evapotranspiration as a rate [m/s]
			     double const& s_raf_, double const& t_raf_, // runoff attentuation feature
			     double const& r_sf_max_, double const& c_sf_, // surface store parameters
			     double const& s_rz_max_, // root zone store parameters
			     double const& t_d_, // unsaturated zone parameters
			     double const& ln_t0_, double const& c_sz_, double const& m_, double const& D_, double const& m_2_, double const& omega_,// saturated zone parameter
			     int const& opt_
			     ):
  id(id_), s_sf(s_sf_), s_rz(s_rz_), s_uz(s_uz_), s_sz(s_sz_),
  //s_bar(s_bar_), area(area_), width(width_),
  width(width_),
  q_sf_in(q_sf_in_), q_sf_out(q_sf_out_),
  q_sz_in(q_sz_in_), q_sz_out(q_sz_out_),
  e_a(e_a_),
  s_raf(s_raf_), t_raf(t_raf_),
  r_sf_max(r_sf_max_), c_sf(c_sf_),
  s_rz_max(s_rz_max_),
  t_d(t_d_),
  ln_t0(ln_t0_), c_sz(c_sz_), m(m_), D(D_), m_2(m_2_), omega(omega_),
  opt(opt_)
{
  // summary from topgraphic data
  beta = std::atan(s_bar_);
  Dx = area_/width;
  // summary from unsaturated zone
  r_uz_sz_max = 1/t_d;
  // compute summary for saturated zone
  cosbeta_m = std::cos(beta) / m;
  cosbeta_m_2 = std::cos(beta) / m_2;
  switch (opt) {
  case 1: // exponential transmissivity
    l_sz_max = std::exp(ln_t0_)*std::sin(beta);   
    break;
  case 2: // constant celerity
    l_sz_max = c_sz*D;
    break;
  case 3: // bounded exponential
    l_sz_max = std::exp(ln_t0_)*std::sin(beta)*( 1 - std::exp( -D*cosbeta_m ) );
    break;
  case 4: // double exponential
    l_sz_max = std::exp(ln_t0_)*std::sin(beta);
    break;
  }
}

std::pair<double, double> hillslope_hru::courant(double& Dt){
  std::pair<double, double> cr(-99.0,-99.0);
  cr.first = c_sf*Dt/Dx;
  switch (opt) {
  case 1: // exponential transmissivity
    cr.second = Dt*std::exp(ln_t0)*std::sin(2.0*beta)/(2.0*m*Dx);
    break;
  case 2: // constant celerity
    cr.second = c_sz*Dt/Dx;
    break;
  case 3: // bounded exponential
    cr.second = Dt*std::exp(ln_t0)*std::sin(2.0*beta)/(2.0*m*Dx);
    break;
  case 4: // double exponential
    cr.second = Dt*(std::exp(ln_t0)*std::sin(2.0*beta)/(2*Dx))*( omega/m + (1.0-omega)/m_2 );
  }
  
  return cr;
}

void hillslope_hru::init(double& s_rz_0, double& r_uz_sz_0, double& tol, int& max_it){
  double l_sz_in, l_sz;
  double r_uz_sz;
  l_sz_in = q_sz_in / width; // standardise inflows by width
  s_sf = 0.0; // initialise surface store
  s_rz = s_rz_max*s_rz_0; // initialise root zone store
  l_sz = std::min(l_sz_max, l_sz_in + Dx*r_uz_sz_0); // outflow flux under steady state
  r_uz_sz = (l_sz - l_sz_in)/Dx; // solve to find actual r_uz_sz
  // compute saturated zone storage deficit
  switch (opt) {
  case 1: // exponential transmissivity
    s_sz = std::max(0.0, (std::log(l_sz_max) - std::log(l_sz))/cosbeta_m);
    break;
  case 2: // constant celerity
    s_sz = (c_sz*D - l_sz)/c_sz;
    break;
  case 3: // bounded exponential
    s_sz = std::max(0.0,
		    -std::log( (l_sz/(std::exp(ln_t0)*std::sin(beta))) +
			       std::exp( -D*cosbeta_m ) ) / cosbeta_m );
    break;
  case 4: // double exponential
    // this needs a numeric search
    // test for zero flow which is impossible but then set a large s_sz
    if( l_sz <= 0.0 ){
      Rcpp::warning("ID: %i. No lateral flux in initialisation - setting depth to 100");
      s_sz = 100.0;
    }else{
      // test for saturation
      double x(0.0); // declared here so can be passed into fsz to test for saturation
      double lwr(0), upr(0); // upper and lower limits of search
      if( (flz(x)-l_sz) <= 0.0 ){
	// then saturated
	s_sz = 0.0;
      }else{
	// not saturated then determine lwr bound
	lwr = 0.0;
	// need to work upper depth to bracket
	upr = 0.01;
	while( (flz(upr)-l_sz) >= 0.0 ){
	  upr += upr;
	}
	upr = std::min(upr,D); // since loop may result in upr > D
      }
      // bisection
      int it(0);
      double e;
      while( (it <= max_it) and ((upr-lwr)>tol) ){
	x = (upr+lwr)/2.0;
	e = flz(x) - l_sz;
	if(e>=0.0){ lwr=x; }
	if(e<=0.0){ upr=x; }
	it +=1;
      }
      if(it >max_it){
	Rcpp::warning("ID: %i. No root found within %i iterations of initialisation. Difference between bounds is %d",
		      id, max_it, upr-lwr);
      }
      s_sz = (upr+lwr)/2.0;
    }
    break;
  }
  s_uz = std::min( s_sz, r_uz_sz*t_d*s_sz ); // compute unsaturated zone storage

  // compute the outflow
  q_sf_out = width*c_sf*s_sf;
  q_sz_out = width*l_sz;
}

void hillslope_hru::implicit_step(double& pet, double& precip, double& Dt, double& tol, int& max_it, double& ftol)
{
  // NOTE: Throughout the vertical fluxes r_*_* are computed as volumes v_*_* = Dt*r_*_*
  // This stops constant rescaling and speeds up the code
  
  double l_sz, l_sf_in;
  //double r_sf_rz;
  double v_sf_rz; //, v_rz_uz;
  
  // standardise inflows by width
  l_sz_in = q_sz_in / width;
  l_sf_in = q_sf_in / width;

  // set ratio of Dt_Dx used
  Dt_Dx = Dt/Dx;

  // single HRU mass balance for development
  //double mass_ballance = s_sf + s_rz + s_uz - s_sz + Dt*precip + Dt_Dx * (l_sz_in + l_sf_in);
  

  // compute first downward flux estimate from surface and root zone
  // these are given as \check{r} in documentation
  // These are computed a volumes - that is Dt*\check{r}
  v_sf_rz = std::min( Dt*r_sf_max , s_sf + Dt_Dx*l_sf_in );
  v_rz_uz = std::max(0.0 ,
		     (s_rz + Dt*(precip - pet) + v_sf_rz - s_rz_max));

  // initialise the bands for the search
  double lwr(0.0), upr(D);

    
  // compute flow estimate
  double z(0.0); // declared here so can be passed into fsz to test for saturation
  if( fsz(z,Dt) >= 0.0 ){
    // then saturated
    l_sz = l_sz_max;
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
	lwr = upr;
	upr += upr;
	fupr = fsz(upr,Dt);
      }
      upr = std::min(upr,D); // since loop may result in upr > D
    }

    int it(0);
    double e(2*ftol);
    while( (it <= max_it) and ( ((upr-lwr)>tol) or (std::abs(e)>=ftol) ) ){
      z = (upr+lwr)/2.0;
      e = fsz(z,Dt);
      if(e<=0.0){ lwr=z; }
      if(e>=0.0){ upr=z; }
      it +=1;
    }
    if(it >max_it){
      Rcpp::warning("ID: %i. No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		    id, max_it, upr-lwr,e);
    }
    z = (upr+lwr)/2.0;
    l_sz = flz(z);
  }

  // solve for saturated zone storage and v_rz_uz
  z = 0.0;
  if( hsz(z,l_sz,Dt) >= 0.0 ){
    // then saturated
    v_rz_uz = ( s_sz + Dt_Dx*(l_sz - l_sz_in) - s_uz );
    s_sz = 0.0;
  }else{
    z = D;
    if( hsz(z,l_sz,Dt) >= 0.0 ){
      // set z to be z_C in vignette
      z = std::max(0.0,s_uz + v_rz_uz - (Dt/t_d));
      if( hsz(z,l_sz,Dt) >= 0.0 ){
	z = s_sz + (Dt/Dx)*(l_sz - l_sz_in) - (Dt/t_d);
      }else{
	double bb = - ( s_sz + Dt_Dx*(l_sz - l_sz_in) - (Dt/t_d));
	double cc = (Dt/t_d) * (s_uz + v_rz_uz - s_sz - Dt_Dx*(l_sz - l_sz_in));
	z = ( -bb + std::sqrt( std::pow(bb,2.0) - 4*cc ) )/2.0;
      }
    }else{
      // have reached D
      z = D ;
      l_sz = ( D - s_sz + Dt_Dx*l_sz_in + Dt*std::min((1/t_d),(s_uz + v_rz_uz) / (t_d*D +Dt)) )/ Dt_Dx ;
    }
    
    s_sz = z;
    v_rz_uz = std::min( v_rz_uz, (s_sz + Dt/t_d - s_uz) );
  }
  
  // solve unsaturated zone
  s_uz = ( (t_d*s_sz)/(t_d*s_sz + Dt) ) * (s_uz+v_rz_uz);
  
  // compute revised r_sf_rz
  v_sf_rz = std::min( v_sf_rz,
		      s_rz_max - s_rz - Dt*(precip-pet)+v_rz_uz
		      );
  // solve for root zone
  s_rz = ( s_rz + Dt*precip + v_sf_rz - v_rz_uz ) * ( s_rz_max / (s_rz_max + Dt*pet) );

  // solve for surface

  //s_sf = ( Dx / (Dx + Dt*c_sf) ) * ( s_sf + Dt_Dx*l_sf_in - v_sf_rz );

  
  double s_sf_tmp = s_sf + Dt_Dx*l_sf_in - v_sf_rz; // tempory middle calc
  s_sf = s_sf_tmp / (1 + (Dt/t_raf)) ; // ) * s_sf_tmp; // valid if new value of s_sf is less then s_raf
  if ( (s_sf > s_raf) | (s_raf==0.0) ){  
    s_sf = ( Dx / (Dx + Dt*c_sf) ) * ( s_sf_tmp  + ((Dt_Dx*c_sf)-(Dt/t_raf))*s_raf );
  }
  

  // compute actual pet loss
  e_a = pet*(s_rz/s_rz_max);

  // // single HRU mass balance for development
  // mass_ballance = mass_ballance - Dt_Dx*(l_sz + c_sf*s_sf) - Dt*e_a;
  // mass_ballance = mass_ballance - (s_sf + s_rz + s_uz - s_sz);
  // if( std::abs(mass_ballance) > 1e-6){
  //     Rcpp::Rcout << "At end: " << mass_ballance << " : " << mth << std::endl;
  //     Rcpp::Rcout << s_sf  << " : " << s_rz << " : " << s_uz << " : " << s_sz << " : " << l_sz << std::endl;
  // }
  
  // compute the outflow
  q_sf_out = width*(c_sf*std::max(0.0,s_sf-s_raf) + (Dx/t_raf)*std::min(s_raf,s_sf) ); //width*Dx = area
  q_sz_out = width*l_sz;

}

double hillslope_hru::flz(double& x){ // compute saturated zone outflow
  //double l = l_sz_max*std::exp(-x*cosbeta_m);
  double l=-99.0;
  switch (opt) {
  case 1: // exponential transmissivity
    l = l_sz_max*std::exp(-x*cosbeta_m);
    break;
  case 2: // constant celerity
    l = c_sz*(D-x);
    break;
  case 3: // bounded exponential
    l = l_sz_max * ( (std::exp(-x*cosbeta_m) - std::exp( -D*cosbeta_m )) /
		     (1 - std::exp( -D*cosbeta_m ) ) );
    break;
  case 4: // double exponential
    l = l_sz_max * ( omega*std::exp(-x*cosbeta_m) + (1-omega)* std::exp( -x*cosbeta_m_2 ) );
    break;
  }
  return l;
}

double hillslope_hru::fsz(double& x, double& Dt){ // compute for saturated zone function
  //double l = l_sz_max*std::exp(-x*cosbeta_m);
  double l = flz(x);
  //double r = std::min( r_uz_sz_max , (s_uz + Dt*r_rz_uz)/(t_d*x + Dt) );
  double r = std::min( r_uz_sz_max , (s_uz + v_rz_uz)/(t_d*x + Dt) );
  
  return x - s_sz - Dt_Dx*(l - l_sz_in) + Dt*r;
}

double hillslope_hru::hsz(double& x, double& l, double& Dt){ // compute for saturated zone function
  //double l = l_sz_max*std::exp(-x*cosbeta_m);
  double r = std::min( r_uz_sz_max , (s_uz + v_rz_uz)/(t_d*x + Dt) );
  
  return x - s_sz - Dt_Dx*(l - l_sz_in) + Dt*r;
}

#include "hillslope_hru.h"
#include "Rcpp.h"
//#include <boost/math/tools/roots.hpp>
//#include <boost/bind/bind.hpp>
//#include <boost/bind/placeholders.hpp>

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

  //Dts = Dt;

  //const boost::uintmax_t opt_maxit = max_it; // maximum number of iterations of optimiser
  //boost::uintmax_t opt_it = opt_maxit; // used in call to optimiser
  //int opt_maxit = max_it;
  //int opt_it = max_it;
  
  std::pair<double, double> opt_res; // solution for output of optimiser
  
  
  // standardise inflow to the saturated zone by width
  l_sz_in = q_sz_in / width;

  // compute first downward flux estimate from surface ans root zone
  // these are given as \check{r} in documentation	
  r_sf_rz = std::min( r_sf_max , (s_sf + Dt*q_sf_in/area)/Dt );
  r_rz_uz = std::max(0.0 ,
		     (s_rz + Dt*(precip + r_sf_rz - pet) - s_rz_max)/Dt);

  // setup template function -TODO alter
  //fsz_exp<double> fnc(s_uz[ii], s_sz[ii],
  //		      l_sz_max[ii], cosbeta_m[ii],
  // t_d[ii], Dt,Dx[ii],  l_sz_in, r_rz_uz);

  double lwr(0), upr(D);

  
  // struct hillslope_hru::TerminationCondition tc;
  // struct hillslope_hru::Fsz fsz(s_uz, s_sz,
  // 				l_sz_max, cosbeta_m,
  // 				t_d, D, m,
  // 				Dt, Dx,
  // 				l_sz_in, r_rz_uz,
  // 				type_sz);
  
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
    double e,dlu;
    while( (it < max_it) and ((upr-lwr)>1e-16) ){

      x = (upr+lwr)/2.0;
      e = fsz(x,Dt);
      dlu = upr-lwr;
      //Rcpp::Rcout << it << " " << lwr << " " << upr << " " << dlu << " " << e << " " << std::endl;
      if(e<=0){ lwr=x; }
      if(e>=0){ upr=x; }
      it +=1;
    }
    if(it >max_it){
      Rcpp::Rcout << it << std::endl;
    }
    
    // solve saturated zone depending upon the type
    //opt_res = boost::math::tools::bisect(fsz,lwr, upr, tc,opt_it);
    
    // if(opt_it >= opt_maxit){
    //   std::cout << "Unable to locate solution in chosen iterations:" <<
    // 	" Current best guess is between " << opt_res.first << " and " <<
    // 	opt_res.second << std::endl;
    // }
    s_sz = (upr+lwr)/2.0;
    //s_sz = (opt_res.first+opt_res.second)/2.0; // solution as avaerage of bounds
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

double hillslope_hru::fsz(double& x, double& Dt){ // compute for saturated zone function
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
  double r = std::min( 1.0/t_d , (s_uz + Dt*r_rz_uz)/(t_d*x + Dt) );
  return x - s_sz - Dt*(l/Dx - l_sz_in/Dx - r);
}


// struct hillslope_hru::TerminationCondition  {
//   bool operator() (double min, double max)  {
//     return abs(min - max) <= 1e-16;
//   };
// };
bool hillslope_hru::TerminationCondition::operator() (double min, double max)  {
    return abs(min - max) <= 1e-16;
}

// annoyingly we appear to have to put all the values into the structure so they are fixed..
// using all lower case and dropping _ in names to differentiate
hillslope_hru::Fsz::Fsz(double const& suz_, double const& ssz_,
			double const& lszmax_, double const& cosbetam_,
			double const& td_, double const d_, double const m_,
			double const& dt_, double const& dx_,
			double const& lszin_, double const& rrzuz_,
			int const& typesz_):
  suz(suz_), ssz(ssz_),
  lszmax(lszmax_), cosbetam(cosbetam_),
  td(td_), d(d_), m(m_),
  dt(dt_),dx(dx_),
  lszin(lszin_), rrzuz(rrzuz_),
  typesz(typesz_){}

double hillslope_hru::Fsz::operator()(double const& x){
  double l = lszmax*std::exp(-x*cosbetam);
  // double l;
  // switch (typesz) {
  // case 1: // exponential transmissivity
  //   l = lszmax*std::exp(-x*cosbetam);
  //   break;
  // case 2: // constant celerity
  //   l = std::max( m*(d-x) , 0.0 );
  //   break;
  // case 3: // bounded exponential
  //   l = lszmax * ( std::exp(-x*cosbetam) - std::exp( -d *cosbetam ) );
  //   break;
  // }
  double r = std::min( 1.0/td , (suz + dt*rrzuz)/(td*x + dt) );
  return x - ssz - dt*(l/dx - lszin/dx - r);
}

// a test for passing data to the hillslope class
#include "dsu.h"

hsu::hsu(double& l_sf_, double& s_rz_, double& s_uz_, double& l_sz_,
	 double& q_sf_in_, double& q_sz_in_,double& p_, double& ep_,
	 double& w_, double& Dx_, double& beta_,
	 double& t_sf_, double& k_sf_,
	 double& s_rzmax_,
	 double& t_d_,
	 double& m_, double& ln_t0_,
	 double& timestep_):
  l_sf(l_sf_), s_rz(s_rz_), s_uz(s_uz_), l_sz(l_sz_),
  q_sf_in(q_sf_in_), q_sz_in(q_sz_in_),p(p_), ep(ep_),
  w(w_), Dx(Dx_), beta(beta_), // properties of HSU [m m rad]
  t_sf(t_sf_), k_sf(k_sf_), // properties of surface store 
  s_rzmax(s_rzmax_),  // properties of root zone
  t_d(t_d_),  // properties of unsat zone
  m(m_), ln_t0(ln_t0_),  // properties of sat zone
  timestep(timestep_)
{
  l_szmax = std::exp(ln_t0)*std::sin(beta);
  log_l_szmax = ln_t0 + std::log( std::sin(beta) );
  cosbeta_m = std::cos(beta) /m;
  lambda_szmax = hsu::flambda_sz(l_szmax);
  lambda_sf = timestep/(t_sf*Dx);
  //std::cout << lambda_sf << std::endl;
}

void hsu::step(){
  
  // standardise inflows by width
  l_sf_in = q_sf_in / w;
  l_sz_in = q_sz_in / w;
  
  // max downward flux times
  r_sf_rz = std::min( k_sf , (l_sf + lambda_sf*l_sf_in)/(Dx*lambda_sf) );
  r_rz_uz = std::max( 0.0 , ((s_rz-s_rzmax)/timestep) + p + r_sf_rz - ep );
  
  
  // max volume in unsat from drain down
  max_uz = s_uz + timestep*r_rz_uz;
  
  // maximum possible flow rate from uz to sz
  double r_uz_pos = std::min( 1/t_d , max_uz/timestep );
  // max rate before sz=0
  double r_uz_0 = ( (1+lambda_szmax)*l_szmax - l_sz - lambda_szmax*l_sz_in ) /(lambda_szmax*Dx); //  (lambda_szmax + lambda_szmax*(l_szmax - l_sz_in) - l_sz)/Dx;
  
  //std::cout << "Dx " << Dx << std::endl;
  //std::cout << "lambda_szmax " << lambda_szmax << std::endl;
  //						  std::cout << "l_szmax " << l_szmax << std::endl;
  //						  std::cout << "l_sz " << l_sz << std::endl;
	
  double s_sz = 0.0;
  if( r_uz_0 <= r_uz_pos ){
    // then sz=0
    r_uz_sz = r_uz_0;
    l_sz = l_szmax ;
    s_sz = 0.0;
  }else{
    // need to solve for a solution
    double fcr = fopt(l_sz);
    //std::cout << fcr << std::endl;
    if(fcr!=0.0){
      // initial range estimates
      double upr = l_sz, lwr = 0.0;
      if(fcr < 0){
	lwr = l_sz;
	upr = l_szmax;
      }
      
      // solve root problem
      // if( hsu::fopt(upr) < 0.0 ){
	
      // 	bool tmp_eq = lwr==upr;
	
      // 	std::cout << " lwr = " << lwr << " upr = " << upr << " diff " << tmp_eq << std::endl;
      // 	std::cout << " fsz(l_sz) = " << fsz(l_sz) << " s_uz = " << s_uz << std::endl;
      // 	std::cout << " fsz(l_szmax) = " << fsz(l_szmax) << std::endl;
      // 	std::cout << " pos = " << r_uz_pos << " zero = " << r_uz_0 << std::endl;
      // 	std::cout << "f(cr) " << fcr << " f(lwr) = " << hsu::fopt(lwr) << " f(upr) = " << hsu::fopt(upr) << std::endl;
      // }
      boost::math::tools::eps_tolerance<double> tol(get_digits);
      boost::uintmax_t it = maxit; // Initally our chosen max iterations, but updated with actual.
      std::pair<double, double> r = boost::math::tools::toms748_solve(boost::bind(&hsu::fopt,this,_1), lwr, upr, tol,it);
      if(it >= maxit)
	{ // 
	  std::cout << "Unable to locate solution in chosen iterations:"
	    " Current best guess is between " << r.first << " and " << r.second << std::endl;
	}
      l_sz = std::min(l_szmax,r.first + (r.second - r.first)/2);
    }
    
    s_sz = hsu::fsz(l_sz);
    // std::cout << "got here" << max_uz/(t_d*s_sz + timestep) << " " << 1/t_d << std::endl;
    r_uz_sz = std::min( max_uz/(t_d*s_sz + timestep), 1/t_d);
  }
  //std::cout << max_uz << " " << s_sz << std::endl;
      
  // solve unsaturated zone
  double tuz = r_uz_sz*t_d*s_sz;
  r_rz_uz = (tuz - s_uz)/timestep + r_uz_sz;
  s_uz = tuz;
  
  // solve root zone and surface
  r_sf_rz = std::min( r_sf_rz, ((s_rzmax - s_rz)/timestep) + ep - p + r_rz_uz);
  //std::cout << l_sf_in << " " << r_sf_rz << " " << l_sf << std::endl;
  l_sf = (l_sf + lambda_sf*(l_sf_in - Dx*r_sf_rz)) / (1 + lambda_sf);
  s_rz = (s_rz + timestep*(p+r_sf_rz-r_rz_uz)) / (1 + (ep*timestep/s_rzmax)) ;
  et = ep*(s_rz/s_rzmax);

  //std::cout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
}


std::vector<double> hsu::get_flux(){
  std::vector<double> out({r_sf_rz, r_rz_uz, r_uz_sz, et});
  return out;
}

std::vector<double> hsu::get_q(){
  std::vector<double> out({l_sf*w, l_sz*w});
  return out;
}

// compute lambda_sz from l
double hsu::flambda_sz(double& l){
  return cosbeta_m*l*timestep/Dx;
}

// compute s_sz from l_sz
double hsu::fsz(double& l){
  return (log_l_szmax - std::log(l))/cosbeta_m;
}

// function to optimise to find l_sz
double hsu::fopt(double l){
  double sz = hsu::fsz(l);
  double ruz = std::min( max_uz/(t_d*sz + timestep), 1/t_d);
  double lambda = hsu::flambda_sz(l);
  return l - (l_sz + lambda*(l_sz_in + Dx*ruz))/(1+lambda);
}

void hsu::init(double& s_rz_0, double& r_uz_sz_0){
  
  // standardise inflows by width
  l_sf_in = q_sf_in / w;
  l_sz_in = q_sz_in / w;

  // apply the initial value to the surface and root zones
  l_sf = 0;
  s_rz = s_rzmax*s_rz_0;

  // out flux under steady state
  l_sz = l_sz_in + Dx*r_uz_sz_0;
  if(l_sz > l_szmax ){
    l_sz = l_szmax;
    s_uz = 0.0;
  }else{
    double s_sz = hsu::fsz(l_sz);
    s_uz = std::min(s_sz,r_uz_sz_0*t_d*s_sz);
  }
}

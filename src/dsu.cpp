// a test for passing data to the hillslope class
#include "dsu.h"

hsu::hsu(double& l_sf_, double& s_rz_, double& s_uz_, double& l_sz_,
	 double& q_sf_in_, double& q_sz_in_,double& p_, double& ep_,
	 double& w_, double& Dx_, double& beta_,
	 double& c_sf_, double& k_sf_,
	 double& s_rzmax_,
	 double& t_d_,
	 double& m_, double& ln_t0_,
	 double& Dt_, std::vector<flink>& links_,
	 double& q_sf_out_,double& q_sz_out_ ):
  l_sf(l_sf_), s_rz(s_rz_), s_uz(s_uz_), l_sz(l_sz_),
  q_sf_in(q_sf_in_), q_sz_in(q_sz_in_),p(p_), ep(ep_),
  w(w_), Dx(Dx_), beta(beta_), // properties of HSU [m m rad]
  c_sf(c_sf_), k_sf(k_sf_), // properties of surface store 
  s_rzmax(s_rzmax_),  // properties of root zone
  t_d(t_d_),  // properties of unsat zone
  m(m_), ln_t0(ln_t0_),  // properties of sat zone
  Dt(Dt_), links(links_),
  q_sf_out(q_sf_out_), q_sz_out(q_sz_out_)
{
  l_szmax = std::exp(ln_t0)*std::sin(beta);
  log_l_szmax = ln_t0 + std::log( std::sin(beta) );
  cosbeta_m = std::cos(beta) /m;
  lambda_szmax = hsu::flambda_sz(l_szmax);
  lambda_sf = (c_sf*Dt)/Dx;
  n_link = links.size();
  s_sz = hsu::fopt(l_sz);
  //std::cout << lambda_sf << std::endl;
}

// fully implicit step
void hsu::step(){
  
  // standardise inflows by width
  l_sf_in = q_sf_in / w;
  l_sz_in = q_sz_in / w;
  
  // max downward flux times (superscript pos description)
  r_sf_rz = std::min( k_sf , (l_sf + lambda_sf*l_sf_in)/(Dx*lambda_sf) );
  r_rz_uz = std::max( 0.0 , ((s_rz-s_rzmax)/Dt) + p + r_sf_rz - ep );
   
  // max volume in unsat from drain down
  max_uz = s_uz + Dt*r_rz_uz;
  
  // maximum possible flow rate from uz to sz
  double r_uz_pos = std::min( 1/t_d , max_uz/Dt );
  // max rate before sz=0
  double r_uz_0 = ( (1+lambda_szmax)*l_szmax - l_sz - lambda_szmax*l_sz_in ) /(lambda_szmax*Dx);
	
  //double s_sz = 0.0;
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
      boost::math::tools::eps_tolerance<double> tol(get_digits);
      boost::uintmax_t it = maxit; // Initally our chosen max iterations, but updated with actual.
      std::pair<double, double> r = boost::math::tools::toms748_solve(boost::bind(&hsu::fopt,this,_1), lwr, upr, tol,it);
      if(it >= maxit)
	{ // commented out to pass CRAN tests 
	  //std::cout << "Unable to locate solution in chosen iterations:"
	  //  " Current best guess is between " << r.first << " and " << r.second << std::endl;
	}
      l_sz = std::min(l_szmax,r.first + (r.second - r.first)/2);
    }
    
    s_sz = hsu::fsz(l_sz);
    r_uz_sz = std::min( max_uz/(t_d*s_sz + Dt), 1/t_d);
  }

  // solve unsaturated zone
  double tuz = r_uz_sz*t_d*s_sz;
  r_rz_uz = (tuz - s_uz)/Dt + r_uz_sz;
  s_uz = tuz;
  
  // solve root zone and surface
  r_sf_rz = std::min( r_sf_rz, ((s_rzmax - s_rz)/Dt) + ep - p + r_rz_uz);
  l_sf = (l_sf + lambda_sf*(l_sf_in - Dx*r_sf_rz)) / (1 + lambda_sf);
  // this estimate of l_sf seems correct but returns small negative values due
  // to rounding differences between this expression and computation of r_sf_rz
  l_sf = std::max(l_sf, 0.0);

  s_rz = (s_rz + Dt*(p+r_sf_rz-r_rz_uz)) / (1 + (ep*Dt/s_rzmax)) ;
  et = ep*(s_rz/s_rzmax);

  // tranfer on the outflow
  q_sf_out = w*l_sf;
  q_sz_out = w*l_sz;
  for(uint i =0; i<n_link; ++i){
    links[i].eval(); //q_sf_out, q_sz_out );
  }
  //std::cout << r_sf_rz << " " << r_rz_uz << " " << r_uz_sz << std::endl;
}


// The semi implicit step
void hsu::astep(){
  
  // standardise inflows by width
  l_sf_in = q_sf_in / w;
  l_sz_in = q_sz_in / w;
  
  // max downward flux times (superscript pos description)
  r_sf_rz = std::min( k_sf , (l_sf + lambda_sf*l_sf_in)/(Dx*lambda_sf) );
  r_rz_uz = std::max( 0.0 , ((s_rz-s_rzmax)/Dt) + p + r_sf_rz - ep );
  
  // max volume in unsat from drain down
  max_uz = s_uz + Dt*r_rz_uz;

  // could be set during initialisation
  double alpha = cosbeta_m*Dt/Dx;
 
  // maximum possible flow rate from uz to sz
  double r_uz_pos = std::min( 1/t_d , max_uz/Dt );
  // max rate before sz=0
  double r_uz_0 = (1/(alpha*Dx))*(1+(alpha*l_szmax) - (l_sz/l_szmax) - (alpha*l_sz_in));

  //double s_sz = 0.0;
  if( r_uz_0 <= r_uz_pos ){
    // then sz=0
    r_uz_sz = r_uz_0;
    l_sz = l_szmax ;
    s_sz = 0.0;
  }else{
    // use approximate solution  based on previous storage
    r_uz_sz = std::min(max_uz / (t_d*s_sz + Dt),1/t_d); // flow based on current storage
    // solve for l_sz
    double gamma = alpha*(l_sz_in + Dx*r_uz_sz) - 1;
    l_sz = ( gamma + std::sqrt( std::pow(gamma,2) + 4*alpha*l_sz) )/(2*alpha);
    s_sz = hsu::fsz(l_sz); // storage based on current flow
  }

  // solve unsaturated zone
  double tuz = std::min(s_sz, s_uz + Dt*(r_rz_uz - r_uz_sz));
  r_rz_uz = (tuz - s_uz)/Dt + r_uz_sz;
  s_uz = tuz;
  
  // solve root zone and surface
  r_sf_rz = std::min( r_sf_rz, ((s_rzmax - s_rz)/Dt) + ep - p + r_rz_uz);
  l_sf = (l_sf + lambda_sf*(l_sf_in - Dx*r_sf_rz)) / (1 + lambda_sf);
  // this estimate of l_sf seems correct but returns small negative values due
  // to rounding differences between this expression and computation of r_sf_rz
  l_sf = std::max(l_sf, 0.0);

  s_rz = (s_rz + Dt*(p+r_sf_rz-r_rz_uz)) / (1 + (ep*Dt/s_rzmax)) ;
  et = ep*(s_rz/s_rzmax);

  // tranfer on the outflow
  q_sf_out = w*l_sf;
  q_sz_out = w*l_sz;
  for(uint i =0; i<n_link; ++i){
    links[i].eval();//q_sf, q_sz );
  }
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
  return cosbeta_m*l*Dt/Dx;
}

// compute s_sz from l_sz
double hsu::fsz(double& l){
  return (log_l_szmax - std::log(l))/cosbeta_m;
}

// function to optimise to find l_sz
double hsu::fopt(double l){
  double sz = hsu::fsz(l);
  double ruz = std::min( max_uz/(t_d*sz + Dt), 1/t_d);
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
    s_sz = 0.0;
    s_uz = 0.0;
  }else{
    s_sz = hsu::fsz(l_sz);
    s_uz = std::min(s_sz,r_uz_sz_0*t_d*s_sz);
  }

  // tranfer on the outflow
  q_sf_out = w*l_sf;
  q_sz_out = w*l_sz;
  //for(uint i =0; i<n_link; ++i){
  //  links[i].eval(q_sf, q_sz );
  //}
}

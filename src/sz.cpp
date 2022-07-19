#include "sz.h"

szc::szc(){}; //std::vector<double> const &param, double const &grd, double const &A, double const&w){};
double szc::fD(){ return(-9999.9); };
double szc::fq_szmax(){ return(-9999.9); };
double szc::fc(double& s){ return(-9999.9); };
void szc::init(double &s, double &q, double &qin, double &r, double &Dx){ };

// bounded exponential
szc_bexp::szc_bexp(std::vector<double> const &param, double const &grd, double const &A, double const &w){
  szc();
  const double &t_0(param[0]), &m(param[1]);
  D = param[2];
  // Rcpp::Rcout << "D is " << D << std::endl;
  double beta = std::atan(grd);
  kappa = std::cos(beta) / m ;
  omega = t_0*std::sin(beta);
  qmax = ( w * omega * ( 1 -  std::exp(-kappa*D) ) ) / A ;
};

double szc_bexp::fD(){ return(D); };
double szc_bexp::fq_szmax(){ return(qmax); };

double szc_bexp::fc(double& s){
  if(s > D){ return( 0.0 ); };
  return( kappa * omega *  std::exp(-s*kappa)  );
};

void szc_bexp::init(double &s, double &q, double &qin, double &r, double &Dx){
  q = qin + r;
  // bounds for search of of state
  std::pair<double,double> bnd(0.0,D), fbnd(0.0,0.0);
  double k = Dx / fc(bnd.first);
  fbnd.first = bnd.first + k*(qin+q)/2.0 - D;
  k = Dx / fc(bnd.second);
  fbnd.second = bnd.first + k*(qin+q)/2.0 - D;
  
  if( fbnd.first >= 0.0 ){
    // then saturated
    s = 0.0;
    q = std::max( (2*D)/k - qin, 0.0 );
  }else{
    if(fbnd.second <= 0.0){
      // then fully drained
      s = D;
      
    }else{
      // bisection
      double z = (bnd.first + bnd.second) / 2.0;
      k = Dx / fc(z);
      double e = z + k*(qin+q)/2.0 - D;
      int it(0.0);
      while( (it <= 1000) and ( (bnd.second - bnd.first)>1e-10 ) and (std::abs(e)>1e-10) ){
	if(e<=0.0){ bnd.first = z; fbnd.first = e;}
	if(e>=0.0){ bnd.second = z; fbnd.second = e;}
	//e = fbnd.first / (fbnd.second - fbnd.first);
	//z = (bnd.first * (1 + e)) - (e * bnd.second);
	z = (bnd.second + bnd.first)/ 2.0 ;
	k = Dx / fc(z);
	e = z + k*(qin+q)/2.0 - D;
	it +=1;
      }
      if(it > 1000){
	Rcpp::warning("No solution found within %i iterations. Difference between bounds is %d. Value of f(z) is %d",
		      1000, bnd.second - bnd.first,e);
      }
      s=z;
    }
    // update values cased on uz formulationa and mass balance ;
    // Rcpp::Rcout << "after bisection..." << std::endl;
    // Rcpp::Rcout << bnd.first << " " << bnd.second << " " << z << std::endl;
    // Rcpp::Rcout << fbnd.first << " " << fbnd.second << " " << std::endl;
  }
  
  if( (s < 0) | (s > D) ){
    Rcpp::Rcout << "HRU " << ": s_sz value " << s << " out of bounds" << std::endl;
  }
  if( (q < 0) ){
    Rcpp::Rcout << "HRU " << ": q_sz value " << q << " is less then 0" << std::endl;
  }
  
};

// // bounded double exponential
// // TODO
// qsz_dexp::qsz_dexp(std::vector<double> const &param, double const &grd, double const &A, double const &w){
//   qsz();
//   const double &t_0(param[0]), &m1(param[1]), &m2(param[2]);
//   double beta;
//   omega = param[3]; // weight
//   beta = std::atan(grd);
//   qmax = w*t_0*std::sin(beta);
//   kappa1 = std::cos(beta) / (m1*A) ;
//   kappa2 = std::cos(beta) / (m2*A) ;
// };
// double qsz_dexp::fq(double& s){
//   return( qmax * ( omega*std::exp(-kappa1*s) + (1-omega)*std::exp(-kappa2*s) ) );
// }
// void qsz_dexp::init(double &s, double &q, double &qin, double &r){
//   q = qin + r;
//   if( q >= qmax ){
//     q = qmax;
//     s = 0;
//     r = qmax - qin;
//   }else{
//     std::pair<double,double> bnd(0,0.1), qbnd(qmax,0.0);
//     qbnd.second = fq(bnd.second);
//     while(qbnd.second > q){
//       bnd.second += bnd.second;
//       qbnd.second = fq(bnd.second);
//     }
//     double z, e;
//     while( (qbnd.second - qbnd.first) > 1e-6 ){
//       z = (bnd.second + bnd.first) / 2.0;
//       e = fq(z);
//       if( e>=q) {
// 	bnd.first = z;
// 	qbnd.first = e;
//       }
//       if( e<=q ) {
// 	bnd.second = z;
// 	qbnd.second = e;
//       }
//     }
//     s = (bnd.second + bnd.first) / 2.0;
//   }
// };

// constant celerity/velocity
szc_cnst::szc_cnst(std::vector<double> const &param, double const &grd, double const &A, double const &w){
  szc();
  //const double &vsz(param[0]), &maxH(param[1]);
  celerity = param[0];
  D = param[1];
  qmax = (w*D*celerity) / A;
};
double szc_cnst::fD(){ return( D ); };
double szc_cnst::fq_szmax(){ return( qmax ); };
double szc_cnst::fc(double &s){ return( celerity ); };

void szc_cnst::init(double &s, double &q, double &qin, double &r){
  q = qin + r;
  if( q >= qmax ){
    s = 0;
    r = qmax - q;
    q = qmax;
  }else{
    s = D * ( 1- (q / qmax) );
  }
};

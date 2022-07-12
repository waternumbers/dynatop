#include "q_sz.h"

qsz::qsz(){}; //std::vector<double> const &param, double const &grd, double const &A, double const&w){};
double qsz::fq(double& s){ return(-9999.9); };
double qsz::fc(double& s){ return(-9999.9); };
void qsz::init(double &s, double &q, double &qin, double &r){ };

// bounded exponential
qsz_bexp::qsz_bexp(std::vector<double> const &param, double const &grd, double const &A, double const &w){
  qsz();
  const double &t_0(param[0]), &m(param[1]);
  double beta;
  D = A * param[2];
  // Rcpp::Rcout << "D is " << D << std::endl;
  beta = std::atan(grd);
  kappa = std::cos(beta) / (m*A) ;
  omega = exp(-kappa*D);
  qmax = w*t_0*std::sin(beta) * ( 1 -  omega );
};
double qsz_bexp::fq(double& s){
  if(s > D){ return( 0.0 ); };
  return( qmax * ( ( std::exp(-s*kappa) - omega) / (1-omega) ) );
}
void qsz_bexp::init(double &s, double &q, double &qin, double &r){
  q = qin + r;
  if( q >= qmax ){
    q = qmax;
    s = 0;
    r = qmax - qin;
  }else{
    s = -std::log( (1-omega)*(q/qmax) + omega ) / kappa;
  }
};

// bounded double exponential
// TODO
qsz_dexp::qsz_dexp(std::vector<double> const &param, double const &grd, double const &A, double const &w){
  qsz();
  const double &t_0(param[0]), &m1(param[1]), &m2(param[2]);
  double beta;
  omega = param[3]; // weight
  beta = std::atan(grd);
  qmax = w*t_0*std::sin(beta);
  kappa1 = std::cos(beta) / (m1*A) ;
  kappa2 = std::cos(beta) / (m2*A) ;
};
double qsz_dexp::fq(double& s){
  return( qmax * ( omega*std::exp(-kappa1*s) + (1-omega)*std::exp(-kappa2*s) ) );
}
void qsz_dexp::init(double &s, double &q, double &qin, double &r){
  q = qin + r;
  if( q >= qmax ){
    q = qmax;
    s = 0;
    r = qmax - qin;
  }else{
    std::pair<double,double> bnd(0,0.1), qbnd(qmax,0.0);
    qbnd.second = fq(bnd.second);
    while(qbnd.second > q){
      bnd.second += bnd.second;
      qbnd.second = fq(bnd.second);
    }
    double z, e;
    while( (qbnd.second - qbnd.first) > 1e-6 ){
      z = (bnd.second + bnd.first) / 2.0;
      e = fq(z);
      if( e>=q) {
	bnd.first = z;
	qbnd.first = e;
      }
      if( e<=q ) {
	bnd.second = z;
	qbnd.second = e;
      }
    }
    s = (bnd.second + bnd.first) / 2.0;
  }
};

// constant celerity/velocity
qsz_cnst::qsz_cnst(std::vector<double> const &param, double const &grd, double const &A, double const &w){
  const double &vsz(param[0]), &maxH(param[1]);
  qsz();
  qmax = w*vsz*maxH;
  kappa = w*vsz /A;
};
double qsz_cnst::fq(double& s){
  return( std::max( 0.0, qmax - kappa*s ) );
}
void qsz_cnst::init(double &s, double &q, double &qin, double &r){
  q = qin + r;
  if( q >= qmax ){
    q = qmax;
    s = 0;
    r = qmax - qin;
  }else{
    s = (qmax - q) / kappa;
  }
};

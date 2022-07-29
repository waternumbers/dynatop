#include "sf.h"

// solve 
sfc::sfc(){ }
double sfc::max_vert(double &s, double &qin, double const &Dt){ return( s/Dt + qin ); }
void sfc::solve(double &s, double &q, double &qin, double &r, double const &Dt){ }
void sfc::init(double &s, double &q, double &qin, double &r){ }

// constant celerity
sfc_cnstC::sfc_cnstC(std::vector<double> const &param, double const &Dx_):
  Dx(Dx_), celerity(param[0])
{}

double sfc_cnstC::max_vert(double &s, double &qin, double const &Dt){ return( s/Dt + qin ); }

void sfc_cnstC::solve(double &s, double &q, double &qin, double &r, double const &Dt){
  //Rcpp::Rcout << s << " " << q << " " << qin  << " " << r << " " << Dt << std::endl;
  double z = s + Dt*(qin -r);
  double ain = qin / celerity;
  //Rcpp::Rcout << s << " " << q << " " << qin  << " " << r << " " << Dt << " " << z << std::endl;
  if(z <= Dx*eta*ain){
    // trapezoid has broken down so mass balance with no outflow
    s = z;
    q = 0.0;
  }else{
    // solve assuming outflow
    s = ( (Dx*(1-eta)) / ( Dx*(1-eta) + celerity*Dt ) ) * z ;
    q = (z-s)/Dt;
  }
}

void sfc_cnstC::init(double &s, double &q, double &qin, double &r){
  q = qin - r;
  s = (Dx/celerity) * ( (eta * q)  + ( (1.0 - eta) * qin) );
}


// constant and diffusivity
sfc_cnstCD::sfc_cnstCD(std::vector<double> const &param, double const &Dx_):
  Dx(Dx_), celerity(param[0])
{
  double D = param[1];
  eta = std::max(0.0, 0.5 - D/(celerity*Dx));
}
double sfc_cnstCD::max_vert(double &s, double &qin, double const &Dt){ return( s/Dt + qin ); }
void sfc_cnstCD::solve(double &s, double &q, double &qin, double &r, double const &Dt){
  double z = s + Dt*(qin -r);
  double ain = qin / celerity;
  if(z <= Dx*eta*ain){
    // trapezoid has broken down so mass balance with no outflow
    s = z;
    q = 0.0;
  }else{
    // solve assuming outflow
    s = ( (Dx*(1-eta)) / ( Dx*(1-eta) + celerity*Dt ) ) * z ;
    q = (z-s)/Dt;
  }
}
void sfc_cnstCD::init(double &s, double &q, double &qin, double &r){
  q = qin - r;
  s = (Dx/celerity) * ( (eta * q)  + ( (1.0 - eta) * qin) );
}

// constant celerity with raf
sfc_cnstC_raf::sfc_cnstC_raf(std::vector<double> const &param, double const &Dx_):
  Dx(Dx_), celerity(param[0]),
  sraf(param[1]), traf(param[2])
{}

double sfc_cnstC_raf::max_vert(double &s, double &qin, double const &Dt){ return( s/Dt + qin ); }

void sfc_cnstC_raf::solve(double &s, double &q, double &qin, double &r, double const &Dt){
  double z = s + Dt*(qin - r);  // storage with inflow added
  Rcpp::Rcout << "raf evalutaion " << z << " " << s << " " << qin << " " << r <<std::endl;
  if( z - Dt*(sraf/traf) > sraf ){ // raf fills and some constant celerity flow
    double ain = std::max( 0.0, qin - ((sraf-s)/Dt)-r-(sraf/traf) ) / celerity;
    Rcpp::Rcout << "trapezoid " << ain <<std::endl;
    // then trapezoid
    if( (z-Dt*(sraf/traf)) <= Dx*eta*ain ){
      // trapezoid has broken down so mass balance with raf outflow
      s = z - Dt*(sraf/traf);
      q = (z-s)/Dt;
      Rcpp::Rcout << "broken trapezoid " << s << " " << q << " " << std::endl;
    }else{
      s = sraf + ( (Dx*(1-eta)) / ( Dx*(1-eta) + celerity*Dt ) ) * (z-Dt*(sraf/traf)-sraf) ;
      q = (z-s)/Dt;
      Rcpp::Rcout << "working trapezoid " << s << " " << q << " " << std::endl;
    }
  }else{
    // all water in raf
    s = z / (1+(Dt/traf));
    q = (z-s)/Dt;
    Rcpp::Rcout << "tank " << s << " " << q << " " << std::endl;
  }
}

void sfc_cnstC_raf::init(double &s, double &q, double &qin, double &r){
  q = qin - r;
  Rcpp::Rcout << "raf initialisation " << q << " " << qin << " " << r <<std::endl;
  Rcpp::Rcout << "raf initialisation " << sraf << " " << traf << " " << eta <<std::endl;
  if( q > (sraf/traf) ){
    // then some sort of trapezoid
    s = sraf + (Dx/celerity) * ( (eta * (q-(sraf/traf)))  + ( (1.0 - eta) * std::max(0.0,qin-(sraf/traf))) );
    Rcpp::Rcout << "trapezoid " << s <<std::endl;
  }else{
    // just in raf
    s = q*traf;
    Rcpp::Rcout << "tank " << s <<std::endl;
  }
}

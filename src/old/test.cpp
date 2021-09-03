#include <iostream> 
#include "hillslope_hru.h" 

using namespace std; 

int main() {   
  double const double_const(1.0);
  double double_var(0.0);
  int const int_const(0);
  
  hillslope_hru hs = hillslope_hru(double_var,double_var,double_var,double_var,
				   double_const,double_const,double_const,
				   double_var,double_var,
				   double_var,double_var,
				   double_var,
				   double_const,double_const,
				   double_const,
				   double_const,
				   double_const,double_const,double_const,
				   int_const);
  
  cout << "made hillslope" << endl;   
  return 0;
} 

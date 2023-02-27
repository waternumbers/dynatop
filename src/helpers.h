#ifndef HELPERS
#define HELPERS

#include <vector>
#include <utility>
#include <cmath>
#include "hru.h"
#include "Rcpp.h" // this is included since data types are used

std::vector<hru> makeHRUs(Rcpp::List);

Rcpp::List makeStateList(std::vector<hru>&);

class outFlux {
 protected: 
  std::vector<int> out_idx, idx, flux_type;
 public:
  outFlux(std::vector<int>, std::vector<int>, std::vector<int>);
  void apply(std::vector<hru>&, std::vector<double>&, double);
};

#endif

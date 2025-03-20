#include <iostream>
#include <vector>
#include <algorithm>
#include <execution>
#include <thread>
#include "Rcpp.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
//#include <tbb/global_control.h>

void Log(int &Number){
  using namespace std::chrono_literals;
  std::this_thread::sleep_for(1s);
  Rcpp::Rcout << Number << std::endl;
  //Number = std::this_thread::get_id();
  //  return( std::this_thread::get_id() );
}
// [[Rcpp::export]]
std::vector<int> test_function(unsigned int nt) {

  // void Log(int &Number){
  //   using namespace std::chrono_literals;
  //   std::this_thread::sleep_for(1s);
  //   Number = std::this_thread::get_id();
  //   //  return( std::this_thread::get_id() );
  // }
  
  // Determine the number of hardware threads available
  unsigned int num_threads = std::thread::hardware_concurrency();
  Rcpp::Rcout << "Number of hardware threads available: " << num_threads << std::endl;

  nt = std::min(nt,num_threads);
  Rcpp::Rcout << "Number of hardware threads used: " << nt << std::endl;
  // Set the number of threads for parallel execution
  //unsigned int desired_threads = 4; // Example: setting to 4 threads
  //std::execution::parallel_policy par = std::execution::par.with(std::execution::thread_pool(desired_threads));
  // somewhere

#if RCPP_PARALLEL_USE_TBB
  #include <tbb/global_control.h>
  Rcpp::Rcout << "using parallel" << std::endl;
  auto policy = std::execution::par;
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, nt);
#else
  Rcpp::Rcout << "using sequential" << std::endl;
  auto policy = std::execution::seq;
#endif

  // Example vector to sort
  std::vector<int> Numbers = {5, 3, 8, 1, 9, 2, 7, 4, 6};

  //int addMe = 23;
  
  std::for_each(policy,
		Numbers.rbegin(), Numbers.rend(),
		Log);
  return( Numbers );
}


// void Log(int &Number){
//   using namespace std::chrono_literals;
//   std::this_thread::sleep_for(1s);
//   Number = std::this_thread::get_id();
//   //  return( std::this_thread::get_id() );
// }


#include <iostream>
#include <vector>
#include <algorithm>
#include <execution>
#include <thread>
#include "Rcpp.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <tbb/global_control.h>

void Log(int &Number){
  using namespace std::chrono_literals;
  std::this_thread::sleep_for(1s);
  //Number = std::this_thread::get_id();
  //  return( std::this_thread::get_id() );
}
// [[Rcpp::export]]
std::vector<int> test_function() {

  // void Log(int &Number){
  //   using namespace std::chrono_literals;
  //   std::this_thread::sleep_for(1s);
  //   Number = std::this_thread::get_id();
  //   //  return( std::this_thread::get_id() );
  // }
  
  // Determine the number of hardware threads available
  unsigned int num_threads = std::thread::hardware_concurrency();
  Rcpp::Rcout << "Number of hardware threads available: " << num_threads << std::endl;
  
  // Set the number of threads for parallel execution
  //unsigned int desired_threads = 4; // Example: setting to 4 threads
  //std::execution::parallel_policy par = std::execution::par.with(std::execution::thread_pool(desired_threads));
  // somewhere
  tbb::global_control c(tbb::global_control::max_allowed_parallelism, 4);


  // Example vector to sort
  std::vector<int> Numbers = {5, 3, 8, 1, 9, 2, 7, 4, 6};
  
  std::for_each(std::execution::par, //par,
		Numbers.begin(), Numbers.end(),
		Log);
  
  return( Numbers );
}


// void Log(int &Number){
//   using namespace std::chrono_literals;
//   std::this_thread::sleep_for(1s);
//   Number = std::this_thread::get_id();
//   //  return( std::this_thread::get_id() );
// }


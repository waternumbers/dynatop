## R test function
rm(list=ls())
library(Rcpp)
setwd("./build_scripts/threads/")
sourceCpp("example.cpp")
test_function(1)

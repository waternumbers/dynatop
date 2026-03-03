## test cases for the hru
rm(list=ls())
library(Rcpp)
library(tinytest)

sourceCpp("hru.cpp")

## ###################################
## zero states to zero states
h <- new(hru,
         0,0,0,10, # s_sf, s_rz, s_uz, s_sz
         1,1, ## area, W_beta
         c(0,0,1),c(0,0,5/3),c(0,0), ## sf: a,b,s
         0.1, # rz: s_rz_max
         2e-5,0.2,10) # sz: K_0,m,D
h$update_inputs(0,0,0,0)
h$iterate(900,1e-4)
h$get_states()
expect_equal(h$get_states(),c(0,0,0,10))

## ##################################
## stationary sat s_sz saturation - sz inflow


## ##################################
## stationary and s_sz saturation - sf inflow
## required careful choice of K_0 and uz

## ##################################
##


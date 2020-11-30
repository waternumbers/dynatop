## some simple R test code for the single HSU
rm(list=ls())
graphics.off()
library(Rcpp)
sourceCpp("multi_hsu.cpp")
n <- 600
timestep <- 90.0
nstep <- 1

ext <- matrix(as.numeric(0),n,4)
states <- matrix(rep(as.numeric(c(0,0,0,6.555984e-06)),2),2,4,byrow=TRUE)
prop <- matrix(c(1.0,10.0,atan(6.554197e-05),
                 1.310684e+02,Inf,6.683900e-02,
                 1.391436e+05,1.205600e-02,1.845000e+00),
               2,9,byrow=TRUE)
ext_idx <- matrix(c(0,1,2,2),2,2) ## C++ indexes
ext[1:n,1] <- 0.0004147561/4

## all must be non-null but can be zero - or create a sink river channel
flow_dir <- list(list(idx = 2,frc=1.0),
                 list(idx = 3,frc=1.0))

channel_inflow <- matrix(0.0,n,1)
channel_id <- 3

keep_states <- rep(TRUE,n)
state_rec <- rep(list(NULL),n)
##ext[400:n,2] <- 0.0004147561/2
#ext[,3] <- as.numeric(1:500)


multi_hsu_cpp(1:nrow(states),states,prop,flow_dir,
              ext_idx,ext,
              channel_id,channel_inflow,
              keep_states,state_rec,
              timestep,nstep)



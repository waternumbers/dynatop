## some simple R test code for the single HSU
rm(list=ls())
graphics.off()
library(Rcpp)
sourceCpp("single_hsu.cpp")
n <- 600
timestep <- 90.0
nstep <- 1
sr <- matrix(as.numeric(NA),n,4)
ext <- matrix(as.numeric(0),n,4)

## ext is in order q_sf_in q_sz_in p ep
##sim <- 1
#ext[1:n,2] <- 0.0004147561
## simulation 2
#sim <- 2
#ext[1:n,1] <- 0.0004147561
## Simulation 3
##sim <- 3
ext[1:n,1] <- 0.0004147561/4
ext[400:n,2] <- 0.0004147561/2


#ext[,3] <- as.numeric(1:500)

flux <- matrix(as.numeric(NA),n,4)

sr[1,] <- as.numeric(c(1e-06,0,0,6.555984e-06))
## properties in order
## width, delta_x, beta, t_sf, k_sf,s_rzmax,t_d, m,ln_t0,timestep
prop <- c(1.0,10.0,atan(6.554197e-05),
          1.310684e+02,Inf,6.683900e-02,
          1.391436e+05,1.205600e-02,1.845000e+00)

single_hsu_cpp(sr,ext,flux,prop,timestep,nstep)

idx <- 390:400; sr[idx,]

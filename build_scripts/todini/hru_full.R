rm(list=ls())
library(R6)

chn <- R6Class(
    "channel",
    public = list(
        S0 = NA,
        B0 = NA,
        ca = NA, ## cot a
        sa = NA, ## sin angle
        Dx = NA,
        n = NA,
        kappa = NA,
        eta = NA,
        initialize = function(S0,B0,grd,Dx,n){
            self$S0 <- S0
            self$B0 <- B0
            self$ca <- 1/grd
            self$sa <- sin(atan(grd))
            self$Dx <- Dx
            self$n <- n
        },
        update = function(Q){
            ## solve for depth
            Qy <- function(y){ (sqrt(self$S0)/n) * (Ay(y)^(5/3)) / (Py(y)^(2/3)) }
            Ay <- function(y){ (self$B0 + y*self$ca)*y }
            By <- function(y){ self$B0 + 2*y*self$ca }
            Py <- function(y){ self$B0 + 2*(y/self$sa) }
            Qy <- function(y){ (sqrt(self$S0)/self$n) * (Ay(y)^(5/3)) / (Py(y)^(2/3)) }
            cy <- function(y){ (5/3)* (sqrt(self$S0)/self$n) * (Ay(y)^(2/3)) / (Py(y)^(2/3)) *
                                   ( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*self$sa) ) )
            }
            vy <- function(y){ (sqrt(self$S0)/self$n) * (Ay(y)/Py(y))^(2/3) }
            betay <- function(y){ (5/3)*( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*self$sa) ) ) }
            fa <- function(y,a){a - Ay(y)}
            fq <- function(y,q){q - Qy(y)}

            y <- uniroot(fq,c(0,100),q=Q)$root
            if(y<0.00001){
                self$kappa <- self$eta <- 0
            }else{
                self$kappa <- self$Dx / vy(y)
                D <- Q/ (2*By(y)*self$S0)
                Ds <- 2*D / (cy(y)*Dx)
                Ds <- Ds * (vy(y)/cy(y))
                self$eta <- 0.5*(1 - Ds)
##                self$eta <- (vy(y)/cy(y))*(0.5-Ds) ## this is wrong
            }
        }
    )
)

hru <- R6::R6Class(
               "hru",
               public = list(
                   ## states
                   s_sf = NA,
                   q_sf = NA,
                   s_rz = NA,
                   s_uz = NA,
                   s_sz = NA,
                   q_sz = NA,
                   sf = NA,
                   sz = NA,
                   e_sf = NA,
                   e_sz = NA,
                   ## initialisation
                   initialize = function(sf,sz){
                       self$sf <- sf
                       self$sz <- sz
                   },
                   init = function(q_sf_in,q_sz_in,

                       self$chn$update( q_in )
                       self$s_sf <- self$chn$kappa*(self$chn$eta*q_in + (1-self$chn$eta)*q_out)
                       self$q_sf <- q_out
                   },
                   ## evolve
                   evolve = function(q_in,Dt){
                       #browser()
                       q_out <- q_in

                       for(it in 1:10){
                           q_ref <- 0.5*(q_in + q_out)
                           self$chn$update(q_ref)
                           q_out <- max(0, (self$s_sf + (Dt - self$chn$kappa*self$chn$eta)*q_in) / (Dt + self$chn$kappa*(1-self$chn$eta)))
                       }
                       shat <- self$chn$kappa*(self$chn$eta*q_in + (1-self$chn$eta)*q_out)
                       stmp <- self$s_sf + Dt*(q_in - q_out)
                       self$e_sf <- shat - stmp
                       self$q_sf <- q_out
                       self$s_sf <- stmp #shat
                   }
               )
           )


## ##############################################################################
sim_time <- 96*60*60
sim_length <- 100*1000
## function to generate forcing
Qinflow <- function(tt){
    Qbase <- 100
    Qpeak <- 900
    beta <- 16
    Tp <- 24*60*60
    Qbase + (Qpeak-Qbase)*( (tt/Tp)*exp(1-(tt/Tp)) )^beta
}

## model steps
Dt <- 900
Dx <- 2000

## generate time steps
ts <- seq(0,sim_time,by=Dt)
nx <- ceiling(sim_length / Dx)

## initial outflow storage
Qrec <- rep(NA,length(ts)) ## outflow at end of simulation reach

## make HRUs
hrus <- list()
tt <- 1
q0 <- Qinflow(ts[tt])
for(ii in 1:nx){
    hrus[[ii]] <- hru$new(q0, chn$new(0.00025,50,Inf,Dx,0.035))
}

## simulate
for(tt in 2:length(ts)){
    qq <- Qinflow(ts[tt])
    ##if(tt==3){browser()}
    for(ii in 1:nx){
        hrus[[ii]]$evolve(qq,Dt)
         if( hrus[[ii]]$s_sf < 0 ){
            print(paste(c("negative",tt,ii,hrus[[ii]]$s_sf)))
        }
        if( hrus[[ii]]$e_sf > 1e-3 ){
            print(paste(c("error",tt,ii,hrus[[ii]]$e_sf,hrus[[ii]]$s_sf,hrus[[ii]]$q_sf)))
        }
        qq <- hrus[[ii]]$q_sf
    }
    Qrec[tt] <- qq
}

#x11()
#plot(ts/3600,Qinflow(ts),type="l")
#lines(ts/3600,Qrec,col="red")
lines(ts/3600,Qrec,col="blue")

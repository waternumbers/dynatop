rm(list=ls())
library(R6)

chn <- R6Class(
    "channel",
    public = list(
        S0 = NA,
        B0 = NA,
        ca = NA, ## cot a
        sa = NA, ## sin angle
        Cs = NA,
        Ds = NA,
        Dx = NA,
        n = NA,
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
            vy <- function(y){ (sqrt(self$S0)/n) * (Ay(y)/Py(y))^(2/3) }
            betay <- function(y){ (5/3)*( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*self$sa) ) ) }
            fy <- function(y,q){q - Qy(y)}
            y <- uniroot(fy,c(0,100),q=Q)$root
            #print(y)
            beta <- betay(y)
            cel <- cy(y)
            self$Cs <- cel/(beta*self$Dx) ## removed Dt compared to paper
            self$Ds <- ifelse(Q==0,0,Q/(beta*By(y)*self$S0*cel*Dx))
        }
    )
)

hru <- R6::R6Class(
               "hru",
               public = list(
                   ## states
                   q_sf_in = NA,
                   q_sf = NA,
                   s_sf = NA,
                   Cs_sf = NA,
                   Ds_sf = NA,
                   chn = NA,
                   e_sf = NA,
                   ## initialisation
                   initialize = function(qin,chn){
                       self$q_sf <- qin
                       self$q_sf_in <- qin
                       self$chn <- chn
                       self$chn$update( qin )
                       self$Cs_sf <- chn$Cs
                       self$Ds_sf <- chn$Ds
                       if(self$Cs_sf==0){
                           self$s_sf <- 0
                       }else{
                           self$s_sf <- (1/(2*self$Cs_sf))*( (1-self$Ds_sf)*self$q_sf_in + (1+self$Ds_sf)*self$q_sf )
                       }
                   },
                   ## evolve
                   evolve = function(q_in,Dt){
                       q_sf_hat <- q_in #self$q_sf + (q_in - self$q_sf_in)
                       shat <- self$s_sf + Dt*q_in
                       for(it in 1:10){
                           #browser()
                           Qref <- (q_sf_hat + q_in)/2
                           self$chn$update(Qref)
                           Cs <- self$chn$Cs
                           Ds <- self$chn$Ds

                           ## (1/(2*Cs)) * ((1-Ds)*q_in + (1+Ds)*qhat) = shat - Dt*qhat
                           ## ((1-Ds)*q_in + (1+Ds)*qhat) = 2*Cs*shat - 2*Cs*Dt*qhat
                           ## (1+Ds)*qhat = 2*Cs*shat - 2*Cs*Dt*qhat - (1-Ds)*q_in
                           ## (1+Ds+2*Cs*Dt)*qhat = 2*Cs*shat - (1-Ds)*q_in
                           ## qhat = (2*Cs*shat - (1-Ds)*q_in) / (1+Ds+2*Cs*Dt)

                           q_sf_hat <- (2*Cs*shat - (1-Ds)*q_in) / (1+Ds+2*Cs*Dt)

                           ## shat <-
                           ## K <- c(
                           ##     -1+(Dt*Cs)+Ds ,
                           ##     (1+(Dt*self$Cs_sf)-self$Ds_sf)*(Cs/self$Cs_sf),
                           ##     (1-(Dt*self$Cs_sf)+self$Ds_sf)*(Cs/self$Cs_sf)
                           ## ) / (1 + (Dt*Cs) + Ds)
                           ## q_sf_hat <- K[1]*q_in + K[2]*self$q_sf_in + K[3]*self$q_sf
                       }
                       stmp <- (1/(2*Cs))*( (1-Ds)*q_in + (1+Ds)*q_sf_hat )
                       self$e_sf <- self$s_sf + Dt*q_in - Dt*q_sf_hat - stmp
                       self$q_sf <- q_sf_hat
                       self$q_sf_in <- q_in
                       self$s_sf <- stmp
                       self$Cs_sf <- Cs
                       self$Ds_sf <- Ds
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
         if( hrus[[ii]]$s_sf <= 0 ){
            print(paste(c(tt,ii,hrus[[ii]]$s_sf)))
        }
        if( hrus[[ii]]$e_sf > 1e-3 ){
            print(paste(c(tt,ii,hrus[[ii]]$e_sf)))
        }
        qq <- hrus[[ii]]$q_sf
    }
    Qrec[tt] <- qq
}

x11()
plot(ts/3600,Qinflow(ts),type="l")
lines(ts/3600,Qrec,col="red")
##lines(ts/3600,Qrec,col="green")

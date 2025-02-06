## example from Ezio's paper
rm(list=ls())

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

### model steps
Dt <- 900
Dx <- 100

## generate time steps
ts <- seq(0,sim_time,by=Dt)
nx <- ceiling(sim_length / Dx)

## initial outflow storage
Qrec <- rep(NA,length(ts)) ## outflow at end of simulation reach

Q <- rep(NA,nx+1) ## flow at previous time step
Qcur <- rep(NA,nx+1) ## flow at current time step

## need to keep these for next timestep
Cstar <- rep(NA,nx) ## need to keep these for next
Dstar <- rep(NA,nx)


## channel definition
S0 <- rep(0.00025,nx)
n <- rep(0.035,nx)
B0 <- rep(50,nx)
ca <- rep(0,nx)
sa <- rep(1,nx)


Ay <- function(y){ (B0[ii] + y*ca[ii])*y }
By <- function(y){ B0[ii] + 2*y*ca[ii] }
Py <- function(y){ B0[ii] + 2*(y/sa[ii]) }
Qy <- function(y){ (sqrt(S0[ii])/n[ii]) * (Ay(y)^(5/3)) / (Py(y)^(2/3)) }
cy <- function(y){ (5/3)* (sqrt(S0[ii])/n[ii]) * (Ay(y)^(2/3)) / (Py(y)^(2/3)) *
                       ( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*sa[ii]) ) )
}
vy <- function(y){ (sqrt(S0[ii])/n[ii]) * (Ay(y)/Py(y))^(2/3) }
betay <- function(y){ (5/3)*( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*sa[ii]) ) ) }
fy <- function(y,q){q - Qy(y)}

## initialise as steady state
tt <- 1
Q[] <- Qinflow(ts[tt])
for(ii in 1:nx){
    Qref <- (Q[ii+1]+ Q[ii])/2
    y <- uniroot(fy,c(0,100),q=Qref)$root
    beta <- betay(y)
    cel <- cy(y)
    Cstar[ii] <- (cel*Dt)/(beta[tt]*Dx)
    Dstar[ii] <- Qref/(beta*By(y)*S0[ii]*cel*Dx)
}

Qrec[1] <- Q[nx+1]

for(tt in 2:length(ts)){
    Qcur[1] <- Qinflow(ts[tt])
    for(ii in 1:nx){
        Qcur[ii+1] <- Q[ii+1] + (Qcur[ii]-Q[ii])
        for(it in 1:10){
            Qref <- (Qcur[ii+1] + Qcur[ii])/2
            y <- uniroot(fy,c(0,100),q=Qref)$root
            beta <- betay(y)
            cel <- cy(y)
            Cs <- (cel*Dt)/(beta*Dx)
            Ds <- Qref/(beta*By(y)*S0[ii]*cel*Dx)

            K <- c(
                -1+Cs+Ds ,
                (1+Cstar[ii]-Dstar[ii])*(Cs/Cstar[ii]),
                (1-Cstar[ii]+Dstar[ii])*(Cs/Cstar[ii])
            ) / (1 + Cs + Ds)
            Qcur[ii+1] <- K[1]*Qcur[ii] + K[2]*Q[ii] + K[3]*Q[ii+1]
        }
        Cstar[ii] <- Cs
        Dstar[ii] <- Ds
    }

    Qrec[tt] <- Qcur[nx+1]
    Q <- Qcur
}

plot(ts/3600,Qinflow(ts),type="l")
lines(ts/3600,Qrec,col="red")
##lines(ts/3600,Qrec,col="green")

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
Dt <- 1800
Dx <- 2000

## generate time steps
ts <- seq(0,sim_time,by=Dt)
nx <- ceiling(sim_length / Dx)

## initial outflow storage
Qrec <- rep(NA,length(ts)) ## outflow at end of simulation reach

Q <- rep(NA,nx+1) ## flow at previous time step
V <- rep(NA,nx) ## current storage

## channel definition
S0 <- rep(0.00025,nx)
n <- rep(0.035,nx)
B0 <- rep(50,nx)
ca <- rep(0,nx)
sa <- rep(1,nx)


Ay <- function(y){ (B0[ii] + y*ca[ii])*y }
## By <- function(y){ B0[ii] + 2*y*ca[ii] }
Py <- function(y){ B0[ii] + 2*(y/sa[ii]) }
Qy <- function(y){ (sqrt(S0[ii])/n[ii]) * (Ay(y)^(5/3)) / (Py(y)^(2/3)) }
##cy <- function(y){ (5/3)* (sqrt(S0[ii])/n[ii]) * (Ay(y)^(2/3)) / (Py(y)^(2/3)) *
##                       ( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*sa[ii]) ) )
##}
vy <- function(y){ (sqrt(S0[ii])/n[ii]) * (Ay(y)/Py(y))^(2/3) }
##betay <- function(y){ (5/3)*( 1 - ( (4*Ay(y))/(5*By(y)*Py(y)*sa[ii]) ) ) }
fa <- function(y,a){a - Ay(y)}
fq <- function(y,q){q - Qy(y)}

## initialise as steady state
tt <- 1
Q[] <- Qinflow(ts[tt])
for(ii in 1:nx){
    y <- uniroot(fq,c(0,100),q=Q[ii+1])$root
    V[ii] <- Dx*Ay(y)
}

Qrec[1] <- Q[nx+1]


for(tt in 2:length(ts)){
    Q[1] <- Qinflow(ts[tt])
    for(ii in 1:nx){
        shat <- V[ii] + Dt*Q[ii]
        for(it in 1){
            #browser()
            a <- V[ii]/Dx
            y <- uniroot(fa,c(0,100000),a=a)$root
            v <- vy(y)
            V[ii] <- shat / (1 + (Dt*v/Dx))
        }
        Q[ii+1] <- V[ii]*v/Dx
    }
    Qrec[tt] <- Q[nx+1]
}

x11()
plot(ts/3600,Qinflow(ts),type="l")
lines(ts/3600,Qrec,col="red")
##lines(ts/3600,Qrec,col="green")

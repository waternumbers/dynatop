## example from Ezio's paper
rm(list=ls())

## generate inflows and time steps
Qbase <- 100
Qpeak <- 100 #900
beta <- 16
Tp <- 24*60*60
Dt <- 1800
ts <- seq(0,96*3600,by=Dt)
Qin <- Qbase + (Qpeak-Qbase)*((ts/Tp)*exp(1-(ts/Tp)))^beta
plot(ts/3600,Qin,type="l")

## channel
S0 <- 0.00025
n <- 0.035
Dx <- 2000

width <- 50

Qout <- rep(NA,length(Qin))
Qref <- cel <- beta <- Cstar <- Dstar <- y <- rep(NA,length(Qin))

## initialise
Qout[1] <- Qref[1] <- Qin[1] ## steady state

fy <- function(y,q,n,s,w){
    A <- y*w
    wp <- 2*y + w ## rectangular
    qhat <- (sqrt(S0)/n) * (A^(5/3)) / (wp^(2/3))
    q-qhat
}
tt <- 1
y[tt] <- uniroot(fy,c(0,100),q=Qref[tt],n=n,s=S0,w=width)$root
beta[tt] <- (5/3)*( 1- (4/5)*( y[tt]/(width + 2*y[tt])) )
cel[tt] <- beta[tt]*Qref[tt]/(width*y[tt])
Cstar[tt] <- (cel[tt]*Dt)/(beta[tt]*Dx)
Dstar[tt] <- Qref[tt]/(beta[tt]*width*S0*cel[tt]*Dx)

for(tt in 2:length(Qin)){
    Qout[tt] <- Qout[tt-1] + (Qin[tt]-Qin[tt-1])
    for(ii in 1:10){
        Qref[tt] <- (Qin[tt]+Qout[tt])/2
        y[tt] <- uniroot(fy,c(0,100),q=Qref[tt],n=n,s=S0,w=width)$root
        beta[tt] <- (5/3)*( 1- (4/5)*( y[tt]/(width + 2*y[tt])) )
        cel[tt] <- beta[tt]*Qref[tt]/(width*y[tt])
        Cstar[tt] <- (cel[tt]*Dt)/(beta[tt]*Dx)
        Dstar[tt] <- Qref[tt]/(beta[tt]*width*S0*cel[tt]*Dx)

        K <- c(
            -1+Cstar[tt-1]+Dstar[tt-1] ,
            (1+Cstar[tt-1]+Dstar[tt-1])*(Cstar[tt]/Cstar[tt-1]),
            (1-Cstar[tt-1]+Dstar[tt-1])*(Cstar[tt]/Cstar[tt-1])
        ) / (1 + Cstar[tt] + Dstar[tt])
        Qout[tt] <- K[1]*Qin[tt] + K[2]*Qin[tt-1] + K[3]*Qout[tt-1]
    }
}

plot(ts/3600,Qin)
lines(ts/3600,Qout)

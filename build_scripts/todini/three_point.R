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
Rinflow <- function(tt){ rep(0,nx) }

### model steps
Dt <- 1800
Dx <- 1000

## generate time steps
ts <- seq(0,sim_time,by=Dt)
nx <- ceiling(sim_length / Dx)

## initial outflow storage
Qrec <- rep(NA,length(ts)) ## outflow at end of simulation reach

Q <- rep(NA,nx+1) ## flow at previous time step

## need to keep these for next timestep
vel <- rep(NA,nx) ## need to keep these for next
D <- rep(NA,nx) ## cell diffusion
C <- rep(NA,nx) ## cell Courant
S <- rep(NA,nx) ## storage

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
Cy <- function(y){ cy(y)*Dt/Dx }
Dy <- function(y){ Qy(y) / (By(y)*S0[ii]*cy(y)*Dx) }

## initialise as steady state
tt <- 1
Q[] <- Qinflow(ts[tt])
for(ii in 1:nx){
    Qref <- (Q[ii+1]+ Q[ii])/2
    y <- uniroot(fy,c(0,100),q=Qref)$root
    S[ii] <- Ay(y)*Dx
}

Qrec[1] <- Q[nx+1]

for(tt in 2:length(ts)){
    Q[1] <- Qinflow(ts[tt])
    r <- Rinflow(ts[tt])

    for(ii in 1:nx){
        for(it in 1:10){
            Qref <- (Q[ii+1] + Q[ii])/2
            y <- uniroot(fy,c(0,100),q=Qref)$root
            C[ii] <- Cy(y)
            D[ii] <- Dy(y)
            vel[ii] <- vy(y)
            if( vel[ii]==0 | C[ii]==0 ){
                K <- rep(0,3)
            }else{
                K <- c(1,
                       Dt - 0.5*( (Dx/vel[ii])-(Dx*D[ii]/C[ii]) ),
                       Dt
                       ) / ( Dt + 0.5*( (Dx/vel[ii])+(Dx*D[ii]/C[ii]) ) )
            }
            Q[ii+1] <- K[1]*S[ii] + K[2]*Q[ii] + K[3]*r[ii]
            Q[ii+1] <- max(0,Q[ii+1],na.rm=T)
        }
        S[ii] <- S[ii] + Dt*Q[ii] - Dt* Q[ii+1] + Dt*r[ii]
        ##S[ii] <- max(S[ii],0)
    }
    Qrec[tt] <- Q[nx+1]
}

x11()
plot(ts/3600,Qinflow(ts),type="l")
lines(ts/3600,Qrec,col="green")

print(c( max(Qrec), which.max(Qrec)) )
##lines(ts/3600,Qrec,col="green")

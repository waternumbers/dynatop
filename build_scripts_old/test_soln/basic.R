## TVDish solution
rm(list=ls())
#graphics.off()


## get flow from storage and length
fq <- function(s,Dx){
    return( 10.3*s/Dx )
}

n <- 1100
s <- vin <- vr <- vout <- qout <- rep(NA,n)
s1 <- vin1 <- vr1 <- vout1 <- qout1 <- rep(NA,n)
qin <- rep(0,n)
qin[200:500] <- 3
r <- rep(0,n)
r[700:800] <-  2
Dt <- 10
Dx <- 10
## assume kinematic
eta <- 0.5
vel <- 1

## initialise at steady state
qout[1] <- qin[1] + r[1]
s[1] <- (Dx/vel)*(eta*qin[1] + (1-eta)*qout[1])
s1[1] <- s[1]
qout1[1] <- qout[1]

## loop
for(ii in 2:n){
    ## 2nd order solution
    ehat <- eta #min(eta,(Dt*vel)/(2*Dx))
    
    ## solve input volumes
    vin[ii] <- (Dt/2)*(qin[ii] + qin[ii-1])
    vr[ii] <- (Dt/2)*(r[ii] + r[ii-1])

    ## solve for qout
    qout[ii] <- s[ii-1] + (Dt/2)*(qin[ii-1] - qout[ii-1]) +
        qin[ii]*( (Dt/2)- (Dx*ehat/vel) ) + vr[ii]
    
    qout[ii] <- qout[ii] * ( (2*vel) / ( Dt*vel +2*Dx*(1-ehat)) )

    s[ii] <- (Dx/vel)*(ehat*qin[ii] + (1-ehat)*qout[ii])

    vout[ii] <- (Dt/2)*(qout[ii] + qout[ii-1])

    ## first order
    ehat <- eta #min(eta,(Dt*vel)/Dx)
    
    ## solve input volumes
    vin1[ii] <- Dt*qin[ii]
    vr1[ii] <- Dt*r[ii]

    ## solve for qout1
    qout1[ii] <- s1[ii-1] + 
        qin[ii]* ( Dt - (Dx*ehat/vel) ) + Dt*r[ii]
    qout1[ii] <- qout1[ii] * ( vel / ( Dt*vel + Dx*(1-ehat)) )

    s1[ii] <- (Dx/vel)*(ehat*qin[ii] + (1-ehat)*qout1[ii])

    vout1[ii] <- Dt*qout1[ii]


}
graphics.off()
x11()
layout(matrix(1:2,1,2))
matplot(cbind(qin,qout,qout1,r),type="p",xlim=c(180,220))
matplot(cbind(qin,qout1,r),type="l")#,xlim=c(180,220))

x11()
layout(matrix(1:2,1,2))
plot(s[-n]+vin[-1]+vr[-1]-vout[-1]-s[-1])
plot(s1[-n]+vin1[-1]+vr1[-1]-vout1[-1]-s1[-1])                                 

x11()
matplot(cbind(s,s1),type="l")







## ## get flow from storage and length
## fq <- function(s,Dx){
##     return( 10.3*s/Dx )
## }

## n <- 1100
## s <- vin <- vr <- vout <- qout <- qout2 <- rep(NA,n)
## qin <- rep(0,n)
## qin[200:500] <- 3
## r <- rep(0,n)
## r[700:800] <-  2
## Dt <- 4
## Dx <- 100
## ## assume kinematic
## eta <- 0.5

## ## initialise at steady state
## qout[1] <- qin[1] + r[1]
## qout2[1] <- qout[1]
## vin[1] <- qin[1]*Dt
## vout[1] <- qout[1]*Dt

## s[1] <- uniroot(function(x,trgt){fq(x,Dx)-trgt},c(0,100),
##                 trgt = (1-eta)*qout[1] + eta*qin[1] )$root




## ## loop
## for(ii in 2:n){
##     vin[ii] <- Dt*qin[ii]
##     vr[ii] <- Dt*r[ii]

##     s_max <- s[ii-1] + Dt*(qin[ii]-r[ii])

##     beta <- Dt/(1-eta)
##     alpha <- s[ii-1] + Dt*r[ii] + beta*qin[ii]

##     fe <- function(s){ s-alpha+beta*fq(s,Dx) }
##     s[ii] <- uniroot(fe,c(0,s[ii-1]+10),
##                      extendInt = "upX")$root
##     vout[ii] <- s[ii-1] + Dt*r[ii] + Dt*qin[ii] - s[ii]
##     qout[ii] <- vout[ii]/Dt
##     qout2[ii] <- (fq(s[ii],Dx) - eta*qin[ii])/(1-eta)
    
## }


## print( s[1] + sum(vin[-1]) + sum(vr[-1]) - sum(vout[-1]) - tail(s,1) )
## e <- c(0, s[-length(s)] + vin[-1] + vr[-1] - vout[-1] - s[-1])
## idx <- 1:length(s) #190:220 #1:length(s)

## x11()
## par(mfrow=c(4,1))
## plot(s[idx])
## plot(qout[idx]); lines(qin[idx])
## plot(vout[idx]); lines(vin[idx])
## plot(e[idx])

## ## Test the celerity based solution with additional dispersion
## rm(list=ls())
## graphics.off()
## n <- 1100
## s <- rep(NA,n)
## cel <- 1
## disp <- 0*100/2
## qin <- rep(1,n)
## qin[200:500] <- 0
## r <- rep(0.1,n)
## qout <- rep(NA,n)
## qav <- rep(NA,n)
## s[1] <- 10
## qout[1] <- qin[1]
## Dt <- 1
## Dx <- 100

## K <- cel*Dt/Dx
## eta <- min(K,0.5 - (disp/(cel*Dx)))

## qav[1] <- eta*qin[1] + (1-eta)*qout[1]

## qqout <- qout
## qqav <- qav
## ss <- s

## for(ii in 2:n){
    
##     ss[ii] <- max(0, ss[ii-1] + ( Dt / (1+K-eta) )*(qin[ii] + (1-eta)*r[ii] - qqav[ii-1]) )
##     qqout[ii] <- (ss[ii-1] + Dt*(qin[ii]+r[ii]) - ss[ii])/Dt
##     qqav[ii] <- eta*qin[ii] + (1-eta)*qqout[ii]
    
##     qout[ii] <- ( qav[ii-1] + (K-eta)*qin[ii] + K*r[ii] ) / ( 1-eta+K )
##     s[ii] <- s[ii-1] + Dt*(qin[ii] - qout[ii] + r[ii])
##     qav[ii] <- eta*qin[ii] + (1-eta)*qout[ii]
## }

## par(mfrow=c(2,1))
## plot(s); lines(ss)
## plot(qout); lines(qqout)


## ## Classic formulation
## rm(list=ls())
## n <- 1100
## s <- rep(NA,n)
## cel <- 100
## disp <- 100*100/2
## qin <- rep(1,n)
## qin[200:500] <- 0
## r <- rep(0,n)
## qout <- rep(NA,n)
## s[1] <- 1
## qout[1] <- qin[1]
## Dt <- 1
## Dx <- 100

## K <- cel*Dt/Dx
## R <- (2*disp)/(cel*Dx)

## theta <- c( 1+R-K, 1-R+K, K-1+R, 2*K)/(1+R+K)
    
## for(ii in 2:n){
##     qout[ii] <- theta[1]*qout[ii-1] + theta[2]*qin[ii-1] + theta[3]*qin[ii] + theta[4]*r[ii]
##     s[ii] <- s[ii-1] + (Dt/2)*(qin[ii] + qin[ii-1] - qout[ii] - qout[ii-1])
## }

## par(mfrow=c(2,1))
## plot(s)
## plot(qout)

## ## Area formulation
## n <- 1100
## s <- rep(NA,n)
## v <- 1.4
## qin <- rep(1,n)
## r <- rep(0,n)
## qout <- ain <- aout <- eta <- rep(NA,n)
## s0 <- 1
## Dt <- 1
## Dx <- 100

## for(ii in 1:n){
##     ain[ii] <- qin[ii] / v
##     eta[ii] <- min(0.5, ( s0+Dt*(qin[ii]-r[ii]) ) / (Dx*ain[ii]) )
##     lwr <- 0; upr <- ain
##     fopt <- function(x){
##         Dx*( (eta[ii]*ain[ii]) + (1-eta[ii])*x ) - (s0 + Dt*(qin[ii] - r[ii] - v*x))
##     }
##     tmp <- uniroot(fopt, c(0,ain[ii]),
##                    extendInt = "upX", check.conv = FALSE,
##                    tol = 2*.Machine$double.eps, maxiter = 1000, trace = 0)
##     aout[ii] <- tmp$root
##     qout[ii] <- v*aout[ii]
##     s[ii] <- s0 + Dt*(qin[ii]-r[ii]-qout[ii])
##     s0 <- s[ii]
## }


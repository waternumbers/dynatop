## TVDish solution
rm(list=ls())
#graphics.off()

## constant velocity
fqo <- function(s,qi){
    Dx <- 100
    v <- 10.3
    eta <- 0.5
    return( ( (s*v/Dx) - (eta*qi) )/(1-eta) )
}

fE <- function(s,s0,vr,qi,ql,qu,Dt){
    ##determine qo
    qo <- fqo(s,qi)
    qo <- min(qu,max(ql,qo))
    vo <- Dt*qo
    vi <- Dt*qi
    return( s+vo-s0-vi-vr )
}

fs <- function(s,qi,qo){
    fqo(s,qi) - qo
}


n <- 1100
s <- vin <- vr <- vout <- qout <- rep(NA,n)
qin <- rep(1,n)
qin[200:500] <- 3
r <- rep(0,n)
r[700:800] <-  2
Dt <- 4

## initialise at steady state
qout[1] <- qin[1] + r[1]
vin[1] <- qin[1]*Dt
vout[1] <- qout[1]*Dt
s[1] <- uniroot(fs,c(-1,100),
                qi=qin[1], qo=qout[1])$root



## loop
for(ii in 2:n){
    vin[ii] <- Dt*qin[ii]
    vr[ii] <- Dt*r[ii]

    vr[ii] <- max(vr[ii], -(s[ii-1]+vin[ii])) ##??
    ## mock up for no limits
    ql <- 0
    qu <- 1000
        ql <- max(0, min(qin[ii],qout[ii-1]) + min(vr[ii]/Dt,0))
    qu <- max(qin[ii],qout[ii-1]) + max(vr[ii]/Dt,0)


    if( fE(s[ii-1],s[ii-1],vr[ii],qin[ii],ql=0,qu=1000,Dt=Dt)==0 ){
        s[ii] <- s[ii-1]
        qout[ii] <- qout[ii-1]
    }else{
        s[ii] <- uniroot(fE,c(0,s[ii-1]+10),
                         s0=s[ii-1],
                         vr=vr[ii],qi=qin[ii],
                         ql=ql,qu=qu,Dt=Dt,
                         extendInt = "upX")$root
        qout[ii] <- min(qu,max(ql,fqo(s[ii],qin[ii])))
    }
    
    #if(ii == 200){browser()}
    vout[ii] <- s[ii-1] + vin[ii] + vr[ii] - s[ii]
##        Dt*qout[ii]
}

print( s[1] + sum(vin[-1]) + sum(vr[-1]) - sum(vout[-1]) - tail(s,1) )
e <- c(0, s[-length(s)] + vin[-1] + vr[-1] - vout[-1] - s[-1])
idx <- 190:220 #1:length(s)

x11()
par(mfrow=c(4,1))
plot(s[idx])
plot(qout[idx]); lines(qin[idx])
plot(vout[idx]); lines(vin[idx])
plot(e[idx])

## ## test the celerity based solution with additional dispersion
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


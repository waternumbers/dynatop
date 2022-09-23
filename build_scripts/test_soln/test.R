## TVDish solution
rm(list=ls())
graphics.off()

## constant velocity
fqo <- function(s,qi){
    Dx <- 100
    v <- 1.3
    eta <- 0.5
    return( ( (s*v/Dx) - (eta*qi) )/(1-eta) )
}

fE <- function(s,s0,vi,vr,qi,q0,Dt,omega){
    ##determine qo
    qo <- fqo(s,qi)
    vo <- Dt*(omega*qo + (1-omega)*q0)
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
r[700:800] <- -2
Dt <- 4

## initialise at steady state
qout[1] <- qin[1] + r[1]
vin[1] <- qin[1]*Dt
vout[1] <- qout[1]*Dt
s[1] <- uniroot(fs,c(-1,100),
                qi=qin[1], qo=qout[1])$root



## loop
for(ii in 2:n){
    vin[ii] <- Dt*0.5*(qin[ii]+qin[ii-1])
    vr[ii] <- Dt*0.5*(r[ii]+r[ii-1])

    vr[ii] <- max(vr[ii], -(s[ii-1]+vin[ii]))
    ## mock up for no limits
    ql <- 0
    qu <- 1000
    ql <- max(0, min(qin[ii],qin[ii-1],qout[ii-1]) + min(vr[ii]/Dt,0))
    qu <- max(qin[ii],qin[ii-1],qout[ii-1]) + max(vr[ii]/Dt,0)
    
    ## find sl
    sl <- uniroot(fs,c(-1,100),
                  qi=qin[ii], qo=ql,
                  extendInt="upX"
                  )$root
    
    ## find su
    su <- uniroot(fs,c(-1,100),
                  qi=qin[ii], qo=qu,
                  extendInt="upX"
                  )$root

    ## check can evaluate
    if(ql==qu){ ## steady state...
        qout[ii] <- qu
        s[ii] <- s[ii-1]
        omega <- 0.5
    }else{
        if( fE(su,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,1) <=0 ){
            #print(paste(ii, "Upper Q limited"))
            qout[ii] <- qu
            s[ii] <- s[ii-1] + vin[ii] + vr[ii] - Dt*qu
            omega <- 1
        }
        if( fE(su,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,1) > 0 &
            fE(su,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,0.5) <=0 ){   
            qout[ii] <- qu
            s[ii] <- su
            omega <- ( su + (Dt*qout[ii-1]) - s[ii-1] - vin[ii] -vr[ii] ) / (Dt*(qout[ii-1]-qu))
            #print(paste(ii, "Upper Q limited on omega", omega))
        }
        if( fE(su,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,0.5) > 0 &
            fE(sl,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,0.5) < 0 ){
            if(ii %in% c(400,750)){ print(paste(ii, "optim")) }
            ## numeric search
            s[ii] <- uniroot(fE,c(sl,su),
                             s0=s[ii-1],vi=vin[ii],
                             vr=vr[ii],qi=qin[ii],
                             q0=qout[ii-1],Dt=Dt,omega=0.5)$root
            qout[ii] <- fqo(s[ii],qin[ii])
            omega <- 0.5
        }
        if( fE(sl,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,0.5) >= 0 &
            fE(sl,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,1) < 0 ){
            qout[ii] <- ql
            s[ii] <- sl
            omega <- ( sl + (Dt*qout[ii-1]) - s[ii-1] - vin[ii] -vr[ii] ) / (Dt*(qout[ii-1]-ql))
            #print(paste(ii, "Lower Q limited on omega", omega))
        }
        if( fE(sl,s[ii-1],vin[ii],vr[ii],qin[ii],qout[ii-1],Dt,1) >= 0 ){
            qout[ii] <- ql
            s[ii] <- s[ii-1] + vin[ii] + vr[ii] - Dt*ql
            omega <- 1
            #print(paste(ii, "Lower Q limited"))
        }
    }
    #if(ii == 200){browser()}
    vout[ii] <- s[ii-1] + vin[ii] + vr[ii] - s[ii]
##        Dt*(omega*qout[ii] + (1-omega)*qout[ii-1])
}

print( s[1] + sum(vin[-1]) + sum(vr[-1]) - sum(vout[-1]) - tail(s,1) )
e <- c(0, s[-length(s)] + vin[-1] + vr[-1] - vout[-1] - s[-1])
idx <- 1:length(s)

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


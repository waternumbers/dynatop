## TVDish solution
rm(list=ls())
#graphics.off()

## first order with muskingham

## get flow from storage and length
fq <- function(s){
    return( 10.3*s/100 )
}
finit <- function(s,q){ fq(s)-q }

n <- 1100
s <- vin <- vr <- vout <- qout <- q <- rep(NA,n)
qin <- rep(0,n)
qin[200:500] <- 3
r <- rep(0,n)
r[700:800] <-  2
Dt <- 40

## assume kinematic
eta <- 0.5

## initialise at steady state
qout[1] <- qin[1] + r[1]
q[1] <- eta*qin[1] + (1-eta)*qout[1]
s[1] <- uniroot(finit,c(0,100),q=q[1])$root

## loop
for(tt in 2:n){
    ##if(tt==212){ browser() }
    vr[tt] <- Dt*r[tt]
    vin[tt] <- qin[tt]*Dt

    kappa <- s[tt-1] + vr[tt] + vin[tt]/(1-eta)
    omega <- Dt/(1-eta)
    
    qmin <- eta*qin[tt]
    smax <- s[tt-1] + vr[tt] + vin[tt]

    if(smax < kappa - omega*fq(smax)){
        print(tt)
        s[tt] <- smax
        q[tt] <- fq(smax)
        next
    }else{
        
        rng <- c(0,smax)
        shat <- mean(rng)
        e <- Inf
        while(abs(e) > 1e-10){
            qq <- max(qmin,fq(shat))
            ## if(q<qmin){
            ##     idx <- 1
            ##     e <- Inf
            ## }else{
            e <- shat - kappa + omega*qq #fq(shat)
            idx <- ifelse(e<0,1,2)
            ##}
            rng[idx] <- shat
            shat <- mean(rng)
        }
        s[tt] <- shat
        q[tt] <- max(qmin,fq(shat))
    }
    
    qout[tt] <- (s[tt-1] + Dt*(r[tt]+qin[tt]) - s[tt])/Dt ##
    ##qout[tt] <- (q[tt] - eta*qin[tt])/(1-eta)##(kappa - s[tt])/omega
    
    vout[tt] <- Dt*qout[tt]
}

    
graphics.off()
x11()
layout(matrix(1:2,1,2))
matplot(cbind(qin,qout,r,q),type="l",xlim=c(300,400))
plot(s[-n]+vin[-1]+vr[-1]-vout[-1]-s[-1])


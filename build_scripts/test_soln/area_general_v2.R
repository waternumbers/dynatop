## Solution based on Area
rm(list=ls())
#graphics.off()

## first order with muskingham

## get flow from storage and length
fq <- function(a){ 1.1 * (a^(5/3)) } #1.3*a }
fa <- function(q){ (q/1.1)^(3/5) } ##q/1.3 }

eta <- 0.3
Dx <- 3000 #2340 #10000#10
n <- 1100
s <- vin <- vr <- vout <- qout <- rep(NA,n)
qin <- rep(1,n)
qin[200:500] <- 3
r <- rep(0,n)
r[700:800] <-  2
nstep <- 10
Dt <- 900

## initialise at steady state
qout[1] <- qin[1] + r[1]
qq <- eta*qin[1] + (1-eta)*qout[1]
s[1] <- fa(qq)*Dx

## loop
for(tt in 2:n){
    ##if(tt==212){ browser() }
    vr[tt] <- Dt*r[tt]
    vin[tt] <- qin[tt]*Dt

    ## omega <- (vel*(1-eta))/ ( Dx*(1-eta) + vel*Dt )
    ## kappa <- s[tt-1] + vr[tt] + vin[tt]/(1-eta)
    ## q[tt] <- omega*kappa
    ## s[tt] <- q[tt]*Dx/vel
    ## qout[tt] <- (q[tt] - eta*qin[tt]) / (1-eta)
    ## vout[tt] <- Dt*qout[tt]

    vinStep <- vin[tt]/nstep
    vrStep <- vr[tt]/nstep
    DtStep <- Dt/nstep
    ain <- fa(qin[tt])

    ss <- s[tt-1]
    vout[tt] <- 0
    for(ii in 1:nstep){
        ##browser()
        smax <- ss + vinStep + vrStep
        rng <- c(0,smax)

        while( diff(rng) > 1e-10 ){
            shat <- mean(rng)
            ahat <- shat/Dx
            qq <- max( 0, (fq(ahat) - eta*qin[tt])/(1-eta) )
            e <- shat - smax + DtStep*qq
            if( e <=0 ){ rng[1] <- shat } else { rng[2] <- shat }
        }

        ss <- mean(rng)
        ahat <- shat/Dx
        qq <- max( 0, (fq(ahat) - eta*qin[tt])/(1-eta) )
        vout[tt] <- vout[tt] + DtStep*qq
    }
    s[tt] <- ss
    qout[tt] <- qq
}
    
#graphics.off()
x11()
layout(matrix(1:3,1,3))
matplot(cbind(qin,qout,r),type="l",main=paste("flow",nstep))
matplot(cbind(qin,qout,r),type="l",main=paste("flow",nstep),xlim=c(195,203))
plot(s[-n]+vin[-1]+vr[-1]-vout[-1]-s[-1])


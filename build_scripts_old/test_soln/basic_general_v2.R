## TVDish solution
rm(list=ls())
#graphics.off()

## first order with muskingham

## get flow from storage and length
etaOpt <- 0.5
vel <- 1.3
Dx <- 3000 #2340 #10000#10

n <- 1100
s <- vin <- vr <- vout <- qout <- q <- rep(NA,n)
qin <- rep(1,n)
qin[200:500] <- 3
r <- rep(0,n)
r[700:800] <-  2
nstep <- 3
Dt <- 900

Dt*vel/(nstep*Dx)

## initialise at steady state
qout[1] <- qin[1] + r[1]
q[1] <- etaOpt*qin[1] + (1-etaOpt)*qout[1]
s[1] <- Dx*q[1]/vel ## uniroot(finit,c(0,100),q=q[1])$root

## no osscilation if (vel*(1-eta))/ ( Dx*(1-eta) + vel*Dt ) < eta*vel/Dx I think...
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

    ss <- s[tt-1]
    vout[tt] <- 0
    for(ii in 1:nstep){
        ##browser()
        eta <- etaOpt##min(etaOpt, DtStep*vel/Dx)
        omega <- (vel*(1-eta))/ ( Dx*(1-eta) + vel*DtStep )
        kappa <- ss + vrStep + vinStep/(1-eta)
        qq <- omega*kappa
        while(qq < eta*qin[tt]){
            eta <- 0.9*eta #DtStep*vel/Dx
            print(eta)
            omega <- (vel*(1-eta))/ ( Dx*(1-eta) + vel*DtStep )
            kappa <- ss + vrStep + vinStep/(1-eta)
            qq <- omega*kappa
        }
        ss <- qq*Dx/vel
        
        #if(tt==200){browser()}
        qqout <- (qq - eta*qin[tt]) / (1-eta)
        if(qqout<0){ print(qqout) }
        vout[tt] <- vout[tt] + DtStep*qqout
    }
    q[tt] <- qq
    s[tt] <- ss
    qout[tt] <- (q[tt] - eta*qin[tt]) / (1-eta)

    
}
    
#graphics.off()
x11()
layout(matrix(1:3,1,3))
matplot(cbind(qin,qout,r,q),type="l",main=nstep)
matplot(cbind(qin,qout,r,q),type="l",main=nstep,xlim=c(195,203))
plot(s[-n]+vin[-1]+vr[-1]-vout[-1]-s[-1])


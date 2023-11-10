## test compound tank solution for stability... and optimisiom
rm(list=ls())

qin <- c(seq(10,450,length=1000),seq(450,400,length=1000))
r <- rep(100,length(qin))##*sign(rnorm(length(qin)))
rout <- r
qout <- rep(NA,length(qin))
s <- rep(NA,length(qin))
s[1] <- 1000

flt <- function(s0,qi,r,Dt,kappa,eta){
    smax <- s0 + Dt*(qi-r)
    sc <- eta*qi/kappa
    if( smax <= sc ){ return( max(0,smax)) }
    else{ return( (s0 + Dt*( (qi/(1-eta)) - r )) / (1+(Dt*kappa/(1-eta))) ) }
}

kappa1 <- 1/1000; eta1 <- 0.5; s1 <- 4e4 ##3000
kappa2 <- 0.1/1000; eta2 <- 0
Dt <- 900

for(ii in 2:length(qin)){
    if(ii==17){ browser() }
    s0 <- s[ii-1]
    shat <- flt(s0,qin[ii],r[ii],Dt,kappa1,eta1)
    if( shat > s1 ){
        ## then reached top part of channel
        qhat <- ((s1-s0)/Dt) + r[ii]
        if( kappa1*s1/eta1 > qhat ){
            qhat <- ( (1+(Dt*kappa1/(1-eta1)))*s1 - s0 + Dt*r[ii] ) /(Dt/(1-eta1))
        }
        if( qhat > 0 ){
            s2 <- 0
        }else{
            qhat <- 0
#            s2 <- s1 + s0 - (Dt*r[ii]) - (Dt/(1-eta1))*kappa1*s1
            s2 <- s0 - (Dt*r[ii]) - (Dt/(1-eta1))*kappa1*s1 - s1
            
        }
        s2 <- flt(s2,qin[ii]-qhat,0,Dt,kappa2,eta2)
        shat <- s1 + s2
    }
    s[ii] <- shat
    qout[ii] <- (s0 + Dt*(qin[ii]-r[ii]) - shat)/Dt
    if( qout[ii] < 0 ){
        rout[ii] <- rout[ii] + qout[ii]
        qout[ii] <- 0
    }
    
}


layout(matrix(1:2))
plot(s)
matplot(cbind(qin,qout),type="l")

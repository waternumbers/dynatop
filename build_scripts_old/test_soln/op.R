## test compound tank solution with operator splitting
rm(list=ls())

q <- 100
r <- -100

eta_1 <- 0.5
kappa_1 <- 0.4
eta_2 <- 0.2
kappa_2 <- 0.4

s_c <- eta_1*q/kappa_1
s_1 <- 1.1*s_c #max(1.1*s_c,10*q)

s <- seq(0,2*max(s_1,s_c),length=10000)
dsdt <- rep(NA,length(s))
ds1dt <- rep(NA,length(s))
q_ss <- max(0, kappa_1*s_1 + (1-eta_1)*r) ## inflow needed at steady state in lower part
q <- q_ss
#q_ss <- kappa_1*s_1
for(ii in 1:length(s)){
    ss <- s[ii]

#    if(ss > s_1){
        q_1 <- min(q,q_ss)
        q_2 <- q - q_1
#    }else{
#        q_1 <- q
#        q_2 <- 0
#    }
    dsdt[ii] <- q - r - max(0,kappa_1*min(ss,s_1) - eta_1*q_1)/(1-eta_1) -
        max(0,kappa_2*max(ss-s_1,0) - eta_2*q_2)/(1-eta_2)
    ds1dt[ii] <- q_1 - r - max(0,kappa_1*min(ss,s_1) - eta_1*q_1)/(1-eta_1)
}

Dt <- 900
s0 <- s - 900*dsdt

##     }else{
        
    
##     if( ss > s_1 ){
        
##         ## then in upper part of channel
##         q_ss <- kappa_1*s_1 + (1-eta_1)*r ## indlow needed for steady state
##         if( q_ss < 0){ ## then should be zero with lateral inflow to upper store
##             r_2 <- q_ss / (1-eta_1)
##             q_ss <- 0
##         }else{
##             r_2 <- 0
##         }
##         q_1 <- min(q_ss,q)
##         q_2 <- q-q_1
##         ss <- ss - s_1
##         if (ss > eta_2*q_2/kappa_2) {
##             dsdt[ii] <- q_2 - r_2 - (kappa_2*ss - eta_2*q_2)/(1-eta_2)
##         }else{
##             dsdt[ii] <- q_2 - r_2
##         }
##     }else{
##         ## then in lower part
##         if (ss > eta_1*q/kappa_1) {
##             dsdt[ii] <- q - r - (kappa_1*ss - eta_1*q)/(1-eta_1)
##         }else{
##             dsdt[ii] <- q - r
##         }
##     }
## }

plot(s,dsdt)
abline(v=s_1)
abline(v=s_c,col="red")



s[which.min(abs(dsdt))]
s[which.min(abs(ds1dt))]

plot(s0,s,xlim=c(0,max(s0)))
## qin <- c(seq(10,450,length=1000),seq(450,400,length=1000))
## r <- rep(100,length(qin))##*sign(rnorm(length(qin)))
## rout <- r
## qout <- rep(NA,length(qin))
## s <- rep(NA,length(qin))
## s[1] <- 1000

## flt <- function(s0,qi,r,Dt,kappa,eta){
##     smax <- s0 + Dt*(qi-r)
##     sc <- eta*qi/kappa
##     if( smax <= sc ){ return( max(0,smax)) }
##     else{ return( (s0 + Dt*( (qi/(1-eta)) - r )) / (1+(Dt*kappa/(1-eta))) ) }
## }

## kappa1 <- 1/1000; eta1 <- 0.5; s1 <- 4e4 ##3000
## kappa2 <- 0.1/1000; eta2 <- 0
## Dt <- 900

## for(ii in 2:length(qin)){
##     if(ii==17){ browser() }
##     s0 <- s[ii-1]
##     shat <- flt(s0,qin[ii],r[ii],Dt,kappa1,eta1)
##     if( shat > s1 ){
##         ## then reached top part of channel
##         qhat <- ((s1-s0)/Dt) + r[ii]
##         if( kappa1*s1/eta1 > qhat ){
##             qhat <- ( (1+(Dt*kappa1/(1-eta1)))*s1 - s0 + Dt*r[ii] ) /(Dt/(1-eta1))
##         }
##         if( qhat > 0 ){
##             s2 <- 0
##         }else{
##             qhat <- 0
## #            s2 <- s1 + s0 - (Dt*r[ii]) - (Dt/(1-eta1))*kappa1*s1
##             s2 <- s0 - (Dt*r[ii]) - (Dt/(1-eta1))*kappa1*s1 - s1
            
##         }
##         s2 <- flt(s2,qin[ii]-qhat,0,Dt,kappa2,eta2)
##         shat <- s1 + s2
##     }
##     s[ii] <- shat
##     qout[ii] <- (s0 + Dt*(qin[ii]-r[ii]) - shat)/Dt
##     if( qout[ii] < 0 ){
##         rout[ii] <- rout[ii] + qout[ii]
##         qout[ii] <- 0
##     }
    
## }


## layout(matrix(1:2))
## plot(s)
## matplot(cbind(qin,qout),type="l")



## ## operator splitting

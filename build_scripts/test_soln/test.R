## test the celerity based solution with additional dispersion
rm(list=ls())
n <- 1100
s <- rep(NA,n)
cel <- 10
disp <- 000*100/2
qin <- rep(1,n)
qin[200:500] <- 0
r <- rep(0.1,n)
qout <- rep(NA,n)
qav <- rep(NA,n)
s[1] <- 10
qout[1] <- qin[1]
Dt <- 1
Dx <- 100

K <- cel*Dt/Dx
eta <- min(K,0.5 - (disp/(cel*Dx)))

qav[1] <- eta*qin[1] + (1-eta)*qout[1]

qqout <- qout
qqav <- qav
ss <- s

for(ii in 2:n){
    
    ss[ii] <- max(0, ss[ii-1] + ( Dt / (1+K-eta) )*(qin[ii] + (1-eta)*r[ii] - qqav[ii-1]) )
    qqout[ii] <- (ss[ii-1] + Dt*(qin[ii]+r[ii]) - ss[ii])/Dt
    qqav[ii] <- eta*qin[ii] + (1-eta)*qqout[ii]
    
    qout[ii] <- ( qav[ii-1] + (K-eta)*qin[ii] + K*r[ii] ) / ( 1-eta+K )
    s[ii] <- s[ii-1] + Dt*(qin[ii] - qout[ii] + r[ii])
    qav[ii] <- eta*qin[ii] + (1-eta)*qout[ii]
}

par(mfrow=c(2,1))
plot(s); lines(ss)
plot(qout); lines(qqout)


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


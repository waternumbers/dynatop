## netwon iteration for powerlaw tank
rm(list=ls())
Dt <- 900 ## s

## inputs
q_sf_in <- 0 #1e-4 ## m3/s
q_sz_in <- 0 #1e-4
p <- 0 #0.01
e_p <- 0 #0.001
A <- 1
tol <- 1e-4

## states
s_sf <- 0
s_rz <- 0
s_uz <- 0
s_sz <- 0.5

## sf param
a <- c(0,0,1)
b <- c(0,0,5/3)
s <- c(0,0)

## rz param
s_rz_max <- 0.1

## us, sz param
K_0 <- 2e-5
m <- 0.002
W_beta <- 1

## single evolution
v_sf_rz <- s_sf + Dt*q_sf_in
v_rz_uz <- min(A*K_0*Dt, max(0,s_rz + Dt*(p-e_p) + v_sf_rz - s_rz_max))


## This bit isn't working
phi <- A*K_0*exp(-s_sz/m)
Hz <- function(z){
    omega <- min(1, (s_uz + v_rz_uz)/(z + Dt*phi))
    q <- W_beta*m*K_0*exp(-z/m)
    out <- c(
        s_sz + Dt*(q - q_sz_in - phi*omega) - z,
        -(q/m) - 1)
    if(omega < 1){
        out[2] <- out[2] + omega*((Dt*phi)/(z+Dt*phi))
    }
    return(out)
}

z <- max(0,s_uz + v_rz_uz - phi*Dt)
h <- Hz(z)
rng <- c(0,Inf)
if(h[1]>=0){rng[1] <- z}
if(h[2]<=0){rng[2] <- z}
it <- 0
while( (-h[1]>tol) & (it<100) ){
    zz <- z
    z <- z - (h[1]/h[2])
    if( z < rng[1] ){ shat <- (rng[1] + zz)/2 }
    if( z > rng[2] ){ shat <- (rng[2] + zz)/2 }
    h <- Hz(z)
    if(h[1]>=0){rng[1] <- z}
    if(h[2]<=0){rng[2] <- z}
    it <- it+1
}


####
v_uz_sz <- Dt*phi*min(1, (s_uz + v_rz_uz)/(z + Dt*phi))
q_sz <- (z - s_sz + v_uz_sz + q_sz_in) / Dt
s_sz <- z

z <- min(s_sz,s_uz + v_rz_uz - v_uz_sz)
v_rz_uz <- z - s_uz + v_uz_sz
s_uz <- z

v_sf_rz <- min(v_sf_rz,s_rz_max - s_rz - Dt*(p-e_p) + v_rz_uz)
z <- ( s_rz_max / (s_rz_max + e_p*Dt) ) * (s_rz + Dt*p + v_sf_rz - v_rz_uz)
v_ep <- s_rz + v_sf_rz - v_rz_uz + Dt*p - z
s_rz <- z

Sw <- function(w){
    if(w <= s[1]){
        q <- a[1]*(w^b[1])
        dq <- a[1]*b[1]*(w^(b[1]-1))
    }else if(w <= s[1]){
        q <- a[1]*(s[1]^b[1]) + a[2]*((w-s[1])^b[2])
        dq <- a[2]*b[2]*((w-s[1])^(b[2]-1))
    }else{
        q <- a[1]*(s[1]^b[1]) + a[2]*((s[2]-s[1])^b[2]) + a[3]*(w-s[2])*b[3]
        dq <- a[3]*b[3]*((w-s[3])^(b[3]-1))
    }
    out <- c(
        s_sf + Dt*(s - q_sf_in - q) - v_sf_rz - w,
        - Dt*dq - 1
    )
    return(out)
}

z <- s_sf
h <- Sw(z)
rng <- c(0,Inf)
if(h[1]>=0){rng[1] <- z}
if(h[2]<=0){rng[2] <- z}
it <- 0
while( (h[1]>tol) & (it<100) ){
    zz <- z
    z <- z - (h[1]/h[2])
    if( z < rng[1] ){ shat <- (rng[1] + zz)/2 }
    if( z > rng[2] ){ shat <- (rng[2] + zz)/2 }
    h <- Hz(z)
    if(h[1]>=0){rng[1] <- z}
    if(h[2]<=0){rng[2] <- z}
    it <- it+1
}

q_sf <- q_sf_in + (s_sf - v_sf_rz - z)/Dt
s_sf <- z




## ## ## iteration
## ## tol <- 1e-4
## ## rng <- c(0,smax)
## ## shat <- smax / (1 + Dt*a)
## ## eta <- 2*tol
## ## it <- 0
## ## while( (abs(eta)> tol) & (it<100) ){
## ##     eta <- shat
## ##     ## this newtom works
## ##     q <- a*(shat^b)
## ##     f <- smax - Dt*q - shat
## ##     df <- -Dt*b*q/shat - 1 ##-Dt*b*a*(shat^(b-1)) - 1
## ##     shat <- shat - (f/df)
## ##     ## shat <- smax / (1 + Dt*a*(shat^(b-1))) ## simpler but slower
## ##     if( shat < 0 ){ shat <- eta/2 }
## ##     if( shat > smax ){ shat <- (smax+eta)/2 }
## ##     print(shat)
## ##     eta <- eta - shat
## ##     it <- it+1
## ## }
## ## print(shat)

## ## iteration
## tol <- 1e-4
## rng <- c(0,smax)
## shat <- smax / (1 + Dt*a)
## eta <- 2*tol
## it <- 0
## while( (abs(eta)>tol) & (it<100) ){
##     ## this newtom works
##     eta <- shat
##     q <- a*(shat^b)
##     f <- smax - Dt*q - shat
##     df <- -Dt*b*q/shat - 1 ##-Dt*b*a*(shat^(b-1)) - 1
##     if(f>=0){rng[1] <- shat}
##     if(f<=0){rng[2] <- shat}
##     shat <- shat - (f/df)
##     ## shat <- smax / (1 + Dt*a*(shat^(b-1))) ## simpler but slower
##     if( shat < rng[1] ){ shat <- (rng[1] + eta)/2 }
##     if( shat > rng[2] ){ shat <- (rng[2] + eta)/2 }
##     eta <- eta - shat
##     print(shat)
##     it <- it+1
## }

## ## 0.007369111
## ## 0.003522388
## ## [1] 0.002730145
## ## [1] 0.002678513
## ## [1] 0.002678271
## ## [1] 0.002678271
## ## [1] 0.002678271
## ## [1] 0.002678271
## ## 0.002678271
## print(shat)
## ## plot for a range of s
## s <- seq(0.9*shat,1.1*shat,length=10000)
## e <- smax - Dt*a*(s^b) - s
## plot(s,e,type="l"); abline(h=0,col="red"); abline(v=shat,lty=2)


## print(s[which.min(abs(e))])


## ## newton iteration for saturated zone
## rm(list=ls())

## s0 <- 3
## umax <- 2
## qr <- 0.0001
## Dt <- 900
## qin <- 0
## m <- 0.02

## u <- function(s){
##     ##pmin(s, (umax*s) / (s + Dt*qr) )
##     pmin(s, umax / (1 + qr*exp(-s/m)/s) )
## }

## ## look at gradient of u(s) for different s

## s <- seq(0,2*umax,length=1000)
## us <- u(s)
## plot(s,us,type="l")
## #abline(v=u0)

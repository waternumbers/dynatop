## netwon iteration for powerlaw tank
rm(list=ls())
library(Rcpp)
Dt <- 900 ## s

## inputs
q_sf_in <- 0 #1e-4 ## m3/s
q_sz_in <- 0 #3.2834e-07 #1e-4
p <- 0 #0.01
e_p <- 0 #0.001
A <- 1
tol <- 1e-5

## states
s_sf <- 1
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
m <- 0.2
W_beta <- 1
D <- 10

## single evolution
v_sf_rz <- s_sf + Dt*q_sf_in
v_rz_uz <- min(A*K_0*Dt, max(0,s_rz + Dt*(p-e_p) + v_sf_rz - s_rz_max))


## solve saturated zone
phi <- A*K_0*exp(-s_sz/m)
Hz <- function(z){
    omega <- min(1, (s_uz + v_rz_uz)/(z + Dt*phi))
    q <- W_beta*m*K_0*exp(-z/m)
    h <- s_sz + Dt*(q - q_sz_in - phi*omega) - z
    dh <- -(q/m) - 1
    if(omega < 1){
        dh <- dh + omega*((Dt*phi)/(z+Dt*phi))
    }
    c(h,dh)
}
## test 0
z <- 0
h <- Hz(z)
if( h[1] <= 0 ){
    q_sz <- W_beta*m*K_0
    ## some other stuff
    v_uz_sz <- s_sz + Dt*(q_sz - q_sz_in)
    s_sz <- z
}else{
    ## not saturated
    rng <- c(0,1000*m)
    z <- s_sz
    h <- Hz(z)
    if(h[1]>=0){rng[1] <- z}
    if(h[1]<=0){rng[2] <- z}
    it <- 0
    while( ((h[1]>0)|(-h[1]>tol)) & (it<100) ){
        zz <- z
        z <- z - (h[1]/h[2])
        if( z < rng[1] ){ z <- (rng[1] + zz)/2 }
        if( z > rng[2] ){ z <- (rng[2] + zz)/2 }
        h <- Hz(z)
        if(h[1]>=0){rng[1] <- z}
        if(h[1]<=0){rng[2] <- z}
        it <- it+1
    }
    v_uz_sz <- Dt*phi*min(1, (s_uz + v_rz_uz)/(z + Dt*phi))
    q_sz <- (z - s_sz + v_uz_sz + q_sz_in) / Dt
    s_sz <- z
}

print(paste("sz iterations =",it))

#### upward pass
z <- min(s_sz,s_uz + v_rz_uz - v_uz_sz)
v_rz_uz <- z - s_uz + v_uz_sz
s_uz <- z

v_sf_rz <- min(v_sf_rz,s_rz_max - s_rz - Dt*(p-e_p) + v_rz_uz)
z <- ( s_rz_max / (s_rz_max + e_p*Dt) ) * (s_rz + Dt*p + v_sf_rz - v_rz_uz)
v_ep <- s_rz + v_sf_rz - v_rz_uz + Dt*p - z
s_rz <- z

print(paste("v_sf_rz =",v_sf_rz))
Sw <- function(w){
    if(w <= s[1]){
        q <- a[1]*(w^b[1])
        dq <- a[1]*b[1]*(w^(b[1]-1))
    }else{
        if(w <= s[1]){
            q <- a[1]*(s[1]^b[1]) + a[2]*((w-s[1])^b[2])
            dq <- a[2]*b[2]*((w-s[1])^(b[2]-1))
        }else{
            q <- a[1]*(s[1]^b[1]) + a[2]*((s[2]-s[1])^b[2]) + a[3]*(w-s[2])^b[3]
            dq <- a[3]*b[3]*((w-s[2])^(b[3]-1))
        }
    }
    print(paste("q is ",q))
    out <- c(
        s_sf + Dt*(q_sf_in - q) - v_sf_rz - w,
        - Dt*dq - 1
    )
    return(out)
}

z <- s_sf
h <- Sw(z)
print(paste(c(z,h)))
rng <- c(0,Inf)
if(h[1]>=0){rng[1] <- z}
if(h[2]<=0){rng[2] <- z}
it <- 0
while( ((h[1]>0)|(-h[1]>tol)) & (it<100) ){
    zz <- z
    z <- z - (h[1]/h[2])
    if( z < rng[1] ){ z <- (rng[1] + zz)/2 }
    if( z > rng[2] ){ z <- (rng[2] + zz)/2 }
    h <- Sw(z)
    if(h[1]>=0){rng[1] <- z}
    if(h[1]<=0){rng[2] <- z}
    it <- it+1
    print(paste(c(z,h)))
}
print(paste("R sf iteration =",it))
q_sf <- q_sf_in + (s_sf - v_sf_rz - z)/Dt
s_sf <- z



library(Rcpp)

sourceCpp("newton.cpp")
nout <- newton(c(1,0,0,0.5),p,e_p,q_sf_in,q_sz_in,900)



sourceCpp("hru.cpp")

h <- new(hru,
         1,0,0,0.5,
         A,W_beta,
         a,b,s,
         s_rz_max,
         K_0,m,D)
h$iterate(900,1e-4)
mout <- h$get_states()

c(s_sf,s_rz,s_uz,s_sz)
nout
mout

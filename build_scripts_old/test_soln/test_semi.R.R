## Semi implicit solution based on M-C with implicit outflow but implict v,eta
rm(list=ls())
graphics.off()

## inputs and time step
n <- 1100

q_sf_in <- rep(1,n)
q_sf_in[200:500] <- 3

r_sf_rz <- rep(0,n)

Dt <- 30
Dx <- 100
s_sf <- s_sf_2 <- rep(NA,n)
q_sf_out <- rep(NA,n)

q_sf_out[1] <- max(0, q_sf_in[1] - r_sf_rz[1])
v <- eta <- 0.5
s_sf[1] <- s_sf_2[1] <- (Dx/v)*( eta*q_sf_in[1] + (1-eta)*q_sf_out[1] ) 

for(tt in 2:n){
    ## compute explicit velocity estiamtes
    v_sf <- 0.5 # *(s_sf[tt-1]/Dx)^0.1
    eta_sf <- 0.5
    kappa_sf <- ( (1-eta_sf)*Dx ) / ( ((1-eta_sf)*Dx) + Dt*v_sf )
    rho_sf <- eta_sf/(1-eta_sf)
    tmp <- s_sf[tt-1] + Dt*(q_sf_in[tt] - r_sf_rz[tt])

    s_sf[tt] <- kappa_sf * (tmp + Dt*rho_sf*q_sf_in[tt])
    
    q_sf_out[tt] <- (tmp - s_sf[tt])/Dt
    if( q_sf_out[tt] < 0 ){
        s_sf[tt] <- tmp
        q_sf_out[tt] <- 0
    }
    s_sf_2[tt] <- (Dx/v_sf)*( (eta_sf*q_sf_in[tt]) + ((1-eta_sf)*q_sf_out[tt]) )
}

##s_sf_2 <- (Dx/v)*( (eta*q_sf_in) + ((1-eta)*q_sf_out) )
plot( s_sf_2, type="l",ylim=c(0,1.05*max(s_sf_2))); points(s_sf)
matplot(cbind(q_sf_in,q_sf_out),type="l")


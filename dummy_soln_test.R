## simple test of stability
rm(list=ls())
graphics.off()

f_rz <- 0.001
f_in <- 0.001

Dt <- 90
n <- ceiling(5000 / Dt)
a <- 0.01; b <- 0.1
T <- 900

u <- s <- f_out <- rep(NA,n)
f_rz_in <- rep(f_rz,n)

u[1] <- 0
s[1] <- 0.5

ui <- u
si <- s
f_rz_in_i <- rep(f_rz,n)


##system.time({
for(tt in 2:n){
   rho = Dt*a*(s[tt-1]^b)
   Ts = T*s[tt-1]
   s[tt] = min(1,(s[tt-1] + Dt*f_in + (Dt/(Ts+Dt))*(u[tt-1] + Dt*f_rz_in[tt]))/(1+rho))
   v_uz_sz = s[tt] - s[tt-1] - Dt*f_in + rho*s[tt]
   u[tt] = min( 1 - s[tt],  u[tt-1] + Dt*f_rz_in[tt] - v_uz_sz )
   f_rz_in[tt] = (u[tt] - u[tt-1] + v_uz_sz)/Dt

   fopt = function(x){ ( (si[tt-1] + Dt*f_in + (Dt/((T*x)+Dt))*(ui[tt-1] + Dt*f_rz_in[tt]))/(1+Dt*a*(x^b)) ) - x }
   if ( fopt(1) >= 0 ){ si[tt] <- 1 }
   else { si[tt] <- uniroot(fopt,c(0,1))$root }
   v_uz_sz = si[tt] - si[tt-1] - Dt*f_in + Dt*a*(si[tt]^b)*si[tt]
   ui[tt] = min( 1 - si[tt],  ui[tt-1] + Dt*f_rz_in[tt] - v_uz_sz )
   f_rz_in_i[tt] = (ui[tt] - ui[tt-1] + v_uz_sz)/Dt
}
##})
#plot((0:(n-1))*Dt,f_rz_in)

plot((0:(n-1))*Dt,s,ylim=c(0,1),type="l",col="red",main="s")
lines((0:(n-1))*Dt,si)

x11();
plot((0:(n-1))*Dt,u,ylim=c(0,1),type="l",col="red",main="u")
lines((0:(n-1))*Dt,ui)

print(n/1e5)

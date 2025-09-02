rm(list=ls())

Dt <- 900
uz <- 0.1

t_d <- 13

q_sz_in <- 0
q_szmax <- 123
psi <- 233
sz <- 0.1
A <- 1
v_rz_uz <- 0.8*Dt*A/t_d

z <- seq(0.03,0.06,length=1000)
r_u <- A*(uz+v_rz_uz)/((t_d*z)+(A*Dt))
idx <- r_u > (A/t_d)
r_u[idx] <- A/t_d
drdz <- ( (-t_d)/(t_d*z + A*Dt) ) * r_u * (!idx)

q_sz <- q_szmax*exp(-psi*z)
dqdz <- -psi*q_sz

f <- sz - (Dt*q_sz_in) - Dt*r_u + Dt*q_sz - z
dfdz <- -Dt*drdz + Dt*dqdz - 1

par(mfrow=c(3,1)); plot(z,f,type="l"); plot(z,dfdz,type="l");lines(head(z,-1),diff(f)/diff(z),col="red",lty=2); plot(z,r_u,type="l"); abline(h =A/t_d,col="red",lty=2)

range(drdz)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
## L is length of channel
## vch is channel velocity
## Dt is the time step
## tau0 is the time delay between the gauge and foot of the reach
fp <- function(L,vch,Dt,tau0){
  tauL <- L/vch
  rL <- floor((tau0+tauL)/Dt)
  irL <- rL+1 ## index in vector since R starts at 1 not zero
  b <- rep(0,irL+1)
  b[irL] <- ((rL+1)*Dt - tau0-tauL)/Dt
  b[irL+1] <- (tau0+tauL - rL*Dt)/Dt
  return(b)
}

## -----------------------------------------------------------------------------
fp(L=100,vch=0.5,Dt=200,tau0=0)

## -----------------------------------------------------------------------------
fp(L=100,vch=1,Dt=200,tau0=600)

## -----------------------------------------------------------------------------
fp(L=100,vch=1.5,Dt=200,tau0=600)

## -----------------------------------------------------------------------------
fd <- function(L,vch,Dt,tau0){
  tauL <- L/vch
  r0 <- floor(tau0/Dt)
  ir0 <- r0+1 ## index in vector since R starts at 1 not zero
  rL <- floor((tau0+tauL)/Dt)
  irL <- rL+1 ## index in vector since R starts at 1 not zero
  b <- rep(0,irL+1)
  
  if(rL>r0){
    b[ir0:(irL-1)] <- Dt # inital values valid unless over written
    b[ir0] <- ( ((r0+1)*Dt - tau0)^2 ) / (2*Dt)
    b[ir0+1] <- b[ir0] + ( ((r0+1)*Dt - tau0)*(tau0 - (r0*Dt)) ) / Dt
    b[irL+1] <- ( (tau0+tauL - (rL*Dt))^2 ) / (2*Dt)
    b[irL] <- b[irL] + b[irL+1] + ( ((tau0+tauL - (rL*Dt))*((rL+1)*Dt - tau0-tauL)) / Dt ) ## added to self since rL could equal r0+1
    if( rL > (r0+1) ){
      b[ir0+1] <- b[ir0+1] + Dt/2
      b[irL] <- b[irL] + Dt/2
    }
  }else{
    b[ir0] <- (tauL/Dt)*( (r0+1)*Dt - tau0 - (tauL/2) )
    b[ir0+1] <- (tauL/Dt)*(tau0 + (tauL/2) - r0*Dt )
  }
  
  return(b*vch/L)
}

## -----------------------------------------------------------------------------
b <- fd(L=200,vch=1,Dt=200,tau0=0)
sum(b)
barplot(b)

## -----------------------------------------------------------------------------
b <- fd(L=200,vch=1.5,Dt=200,tau0=0)
sum(b)
barplot(b)

## -----------------------------------------------------------------------------
b <- fd(L=200,vch=0.5,Dt=200,tau0=0)
sum(b)
barplot(b)

## -----------------------------------------------------------------------------
b <- fd(L=600,vch=1,Dt=200,tau0=800)
sum(b)
barplot(b)

## -----------------------------------------------------------------------------
b <- fd(L=600,vch=1,Dt=200,tau0=850)
sum(b)
barplot(b)


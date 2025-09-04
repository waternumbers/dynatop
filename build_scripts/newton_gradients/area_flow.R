rm(list=ls())

n <- 0.03
B0 <- 5
tan_alpha <- 0.01
q_crit <- 15
grd <- 0.001
Dx <- 1
 ##         area        Dx  gradient
## 373.0799 1638.0000    0.0010

ca <- 1/tan_alpha
sa <- sin(atan(tan_alpha))
beta <- sqrt(grd)/n

Ay <- function(y){ return( (B0*y) + pmax(0.0,y-y_crit)*ca*pmax(0.0,y-y_crit) ); }
#Py <- function(y){ return( B0 + 2.0*pmin(y,y_crit) + 2.0*(pmax(0.0,y-y_crit)/sa) ); }

Qy <- function(y){
    Atri <- 0.5*pmax(0.0,y-y_crit)*ca*pmax(0.0,y-y_crit)
    Ptri <- pmax(1e-6,y-y_crit)/sa ## since is y<y_crit Atri=0
    Achn <- B0*y
    Pchn <- B0 + 2*pmin(y,y_crit)
    Qchn <- beta * (Achn^(5/3)) / Pchn^(2/3)
    Qtri <- beta * (Atri^(5/3)) / Ptri^(2/3)
    return(Qchn + 2*Qtri)
}

y_crit <- Inf
y_crit <- uniroot(function(x){q_crit - Qy(x)},c(0,100))$root

y <- seq(0,3,length=10000)
A <- Ay(y)
Q <- Qy(y)
plot(y,A)
abline(a=0,b=B0)

plot(y,Q)
plot(A,Q)

Qs <- function(s){
    A <- s/Dx
    Ac <- s_crit/Dx
    yc <- Ac/B0
    dA <- pmax(0,A-Ac)
    dy <- ( -B0 + sqrt( B0^2 + 4*dA*ca ) ) / (2*ca)
    y <- pmin(yc+dy, A/B0)
    qrect <- beta*((B0*y)^(5/3))/( (B0+2*pmin(yc,y))^(2/3) )
    qtri <- beta*((ca/2)^(5/3))*(sa^(2/3))*(dy^(8/3))
    q <- qrect + 2*qtri
    dqds <- rep(NA,length(qrect))
    idx <- y<yc
    dqds[idx] <- ((5/3)*(q[idx]/y[idx]) - (4/3)*(q[idx]/(B0+2*y[idx])))/ (B0*Dx)
    dqds[!idx] <- ((5/3)*(qrect[!idx]/y[!idx]) + (16/3)*(qtri[!idx]/dy[!idx]))*(1/sqrt(B0^2 + 4*dA[!idx]*ca))/ Dx
    return(list(y=y,
                q=q,
                dqds=dqds))
}

s_crit <- y_crit*Dx*B0

S <- A*Dx
plot(S,Q,type="l",col="red")
tmp <- Qs(S)
lines(S,tmp$q,lty=2)


plot(head(S,-1),diff(Q)/diff(S),type="l",col="red",xlim=c(12,13))
lines(A*Dx,tmp$dqds,lty=2)

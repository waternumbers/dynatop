
## x' = Ax
rm(list=ls())

## a real case
A <- matrix( c(2,0,0,-2,3,-1,4,-2,2),3 )
x0 <- c(2,0,2)

tmp <- eigen(A)
C <- solve(tmp$vectors,x0)

dt <- 1
## my solution
x <- tmp$vectors %*% (C*exp(-tmp$values*dt))
## solution should be
B <- matrix(c(1, 0, 0, -2,1,1,4,-2,1),3)
xhat <- B %*% (solve(B,x0)*exp(-c(2,1,4)*dt))

x-xhat



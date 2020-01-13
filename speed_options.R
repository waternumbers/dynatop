rm(list=ls())

n <- 20000
nxt <- sample(n,n,replace=TRUE)
W <- Matrix::sparseMatrix(1:n,nxt,x=1,dims=c(n,n))
ex <- runif(n)
bnd <- rep(1:2000,length.out=n)
tex <- runif(n)


## create data objects
d1 <- list(ex = ex,
           bnd = bnd,
           tex=tex,
           rd=list())

for(ii in unique(bnd)){
    idx <- which(d1$bnd==ii)
    d1$rd[[ii]] <- list(idx=idx,
                        X=W[idx,])
}


d2 <- list()
for(ii in 1:n){
    d2[[ii]] <- list(ex = ex[ii],
                     tex = tex[ii],
                     rd = list(idx = nxt[ii],
                               w=1))
}

f1 <- function(d,t){
    lex <- rep(0,length(d$ex))
    for(ii in unique(d$bnd)){
        idx <- d$rd[[ii]]$idx
        X <- d$rd[[ii]]$X
        lex_in <- X %*% lex
        ebt <- exp(-t/d$tex[idx])
        ex <- d$ex[idx]*ebt + lex_in*d$tex[idx]*(1-ebt)
        lex[idx] <- d$ex[idx] + lex_in - ex
        d$ex[idx] <- ex
    }
    return(d)
}

f2 <- function(lst,t){
    lex <- rep(0,length(lst))
    for(ii in 1:length(lst)){
        ebt <- exp(-t/lst[[ii]]$tex[1])
        ex <- lst[[ii]]$ex*ebt + lex[ii]*lst[[ii]]$tex*(1-ebt)
        lx <- lst[[ii]]$ex + lex[ii] - ex
        lex[lst[[ii]]$rd$idx] <- lex[lst[[ii]]$rd$idx] + lx*lst[[ii]]$rd$w
        lst[[ii]]$ex <- ex
    }
    return(lst)
}

f2_cpm <- compiler::cmpfun(f2)

#out1 <- f1(d1,1)
#out2 <- f2(d2,1)
#out2_cpm <- f2_cpm(d2,1)

system.time(f1(d1,1))
system.time(f2(d2,1))
system.time(f2_cpm(d2,1))

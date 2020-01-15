## test for a leaf and branch method
rm(list=ls())
library("Rcpp")
## leaf types
transfers <- function(){c("precip"=0,"pet"=0,"lex"=0,"lsz"=0,"qch"=0)}

external <- function(id){
    list(id=id,
         type = "ex",
         series = c("precip"="a","pet"="b"),
         output = transfers())
}

leaf <- function(id,us){
    list(id =id,
         type="leaf",
         uptree=us,
         output=transfers())
}

branch <- function(id,us,w){
    list(id=id,
         type="branch",
         uptree=us,
         weight=w,
         output=transfers())
}
mux <- function(id,us){
    list(id=id,
         type="mux",
         uptree=us,
         weight=rep(1/length(us),length(us)),
         output=transfers())
}

hillslope <- function(id,us,prop,param){
    list(id=id,
         type="hillslope",
         uptree=us,
         output=transfers(),
         state=c("ex"=0,
                 "rz"=0,
                 "uz"=0,
                 "sz"=0,
                 "lsz_in"=0,
                 "lsz"=0,
                 "lsz_max"=0),
         prop=prop,
         param=param)
}

channel <- function(id,us,prop){
    list(id=id,
         type="channel",
         uptree=us,
         output=transfers(),
         prop=prop)
}

test_prop <- c("mapid"=1,"area"=1,"atb_bar"=0.25,"s_bar"=0.003,"delta_x"=1)
test_param <- c("qex_max"=34,"srz_max"=1,"srz_0"=0.5,"ln_t0"=25,"m"=19,"td"=100,"tex"=92)

hru <- list(external("ex"),
            hillslope("h1","ex",test_prop,test_param),
            hillslope("h2","ex",test_prop,test_param),
            mux("m1",c("ex","h1","h2")),
            channel("c1","m1",c("area"=1)))

names(hru) <- vapply(hru,function(x){x$id},character(1))

sourceCpp("landb.cpp")
rcpp_dynatop(hru,c("a"=4,"b"=34),list(time_step=1,omega=07,theta=0.7))

hru2 <- rep(hru,25000)
system.time({rcpp_dynatop(hru2,c("a"=4,"b"=34))})
## system.time({
##     for(ii in names(hru)){
##         #print(ii)
##         if(hru[[ii]]$type=="ex"){
##             ##hru[[ii]] <- evolve_ex(hru[[ii]])
##             hru[[ii]]$output[] <- rnorm(5)
##         }
##         if(hru[[ii]]$type=="leaf"){
##             ## hru[[ii]] <- evolve_leaf(hru[[ii]],hru[[hru[[ii]]$uptree]]$output)
##             hru[[ii]]$output[] <- hru[[hru[[ii]]$uptree]]$output + 3
##             #for(jj in 1:20){
##             #    hru[[ii]]$output <- hru[[ii]]$output + 2
##             #}

##         }

##         if(hru[[ii]]$type=="branch"){
##             #hru[[ii]] <- evolve_branch(hru[[ii]],hru[[hru[[ii]]$uptree]]$output)
##             hru[[ii]]$output[] <- hru[[hru[[ii]]$uptree]]$output*hru[[ii]]$weight
##         }
##         if(hru[[ii]]$type=="mux"){
##             hru[[ii]]$output[] <- 0
##             for(jj in hru[[ii]]$uptree){
##                 hru[[ii]]$output <- hru[[ii]]$output + hru[[jj]]$output
##             }
##         }
##     }
## })

## system.time({tmp <- rcpp_dynatop(hru)})

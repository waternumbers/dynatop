## test for a leaf and branch method
rm(list=ls())
library("Rcpp")
## leaf types
transfers <- function(){c("precip"=0,"pet"=0,"lex"=0,"lch"=0,"qch"=0)}

external <- function(id){
    list(id=id,
         type = "ex",
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
                 "lsz"=0),
         prop=prop,
         param=param)
}

             
hru <- list(external("ex"),
            leaf("l1","ex"),
            mux("m2",c("ex","l1")),
            leaf("l2","m2"),
            branch("b1","l1",0.5),
            mux("m1",c("b1","l2")))

names(hru) <- vapply(hru,function(x){x$id},character(1))

sourceCpp("landb.cpp")
tmp <- rcpp_dynatop(hru)

hru <- rep(hru,100000)
#system.time({
    for(ii in names(hru)){
        #print(ii)
        if(hru[[ii]]$type=="ex"){
            ##hru[[ii]] <- evolve_ex(hru[[ii]])
            hru[[ii]]$output[] <- rnorm(5)
        }
        if(hru[[ii]]$type=="leaf"){
            ## hru[[ii]] <- evolve_leaf(hru[[ii]],hru[[hru[[ii]]$uptree]]$output)
            hru[[ii]]$output[] <- hru[[hru[[ii]]$uptree]]$output + 3
            #for(jj in 1:20){
            #    hru[[ii]]$output <- hru[[ii]]$output + 2
            #}
            
        }
        
        if(hru[[ii]]$type=="branch"){
            #hru[[ii]] <- evolve_branch(hru[[ii]],hru[[hru[[ii]]$uptree]]$output)
            hru[[ii]]$output[] <- hru[[hru[[ii]]$uptree]]$output*hru[[ii]]$weight
        }
        if(hru[[ii]]$type=="mux"){
            hru[[ii]]$output[] <- 0
            for(jj in hru[[ii]]$uptree){
                hru[[ii]]$output <- hru[[ii]]$output + hru[[jj]]$output
            }
        }
    }
#})
    

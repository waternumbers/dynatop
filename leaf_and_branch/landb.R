## test for a leaf and branch method

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

evolve_ex <- function(h){
    h$output[] <- rnorm(5)
    return(h)
}
evolve_leaf <- function(h,input){
    h$output <- input + 3
    return(h)
}
evolve_branch <- function(h,input){
    h$output <- input*h$weight
    return(h)
}
evolve_mux <- function(h,input){
    h$output[] <- 0
    for(jj in input){
        h$output <- h$output + jj
    }
    return(h)
}


hru <- list(external("ex"),
            leaf("l1","ex"),
            leaf("l2","ex"),
            branch("b1","l1",0.5),
            mux("m1",c("b1","l2")))
names(hru) <- vapply(hru,function(x){x$id},character(1))
hru <- rep(hru,10000)
system.time({
    for(ii in names(hru)){
        #print(ii)
        if(hru[[ii]]$type=="ex"){
            ##hru[[ii]] <- evolve_ex(hru[[ii]])
            hru[[ii]]$output[] <- rnorm(5)
        }
        if(hru[[ii]]$type=="leaf"){
                                        #hru[[ii]] <- evolve_leaf(hru[[ii]],hru[[hru[[ii]]$uptree]]$output)
            hru[[ii]]$output[] <- hru[[hru[[ii]]$uptree]]$output + 3
            for(jj in 1:20){
                hru[[ii]]$output <- hru[[ii]]$output + 2
            }
            
        }
        if(hru[[ii]]$type=="branch"){
            #hru[[ii]] <- evolve_branch(hru[[ii]],hru[[hru[[ii]]$uptree]]$output)
            hru[[ii]]$output[] <- hru[[hru[[ii]]$uptree]]$output*hru[[ii]]$weight
        }
        if(hru[[ii]]$type=="mux"){
            #input <- lapply(hru[ hru[[ii]]$uptree ],function(x){x$output})
            #hru[[ii]] <- evolve_mux(hru[[ii]],input)
            hru[[ii]]$output[] <- 0
            for(jj in hru[[ii]]$uptree){
                hru[[ii]]$output <- hru[[ii]]$output + hru[[jj]]$output
            }
        }
    }
})
    

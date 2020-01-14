## quick code to help debug the latest dynamic topmodel
rm(list=ls())
library(Matrix)
#devtools::load_all("./dynatop")

#mdl <- readRDS("./dynatop/data/Swindale_model.rds")
mdl <- readRDS("Swindale_ordered.rds")
load("./dynatop/data/test_catchment.rda")


for(ii in c("hillslope","channel")){
    names(mdl[[ii]]) <- gsub("_input","_series",names(mdl[[ii]]))
    mdl[[ii]]$precip_series <- "rain"
    mdl[[ii]]$pet_series <- "pet"
}

mdl$channel[['band']] <- 1e4

## mdl$Dsz <- rbind(mdl$Fsz,mdl$Wsz)
## mdl$Dsz <- cbind( matrix(0,nrow(mdl$Dsz),nrow(mdl$Dsz)-ncol(mdl$Dsz)),mdl$Dsz )
## mdl$Wsz <- mdl$Fsz <- mdl$Wch <- NULL

## mdl$hillslope[,'qex_max'] <- 'qex_max_default'
## mdl$param <- c(mdl$param,'qex_max_default'=Inf)

## mdl$hillslope[,'delta_x'] <- 1
## names(mdl$gauge) <- gsub("channel_id","id",names(mdl$gauge))
## names(mdl$point_inflow) <- gsub("channel_id","id",names(mdl$point_inflow))

## standardise column sums
tmp <- colSums(mdl$Dsz)
while(sum(tmp>1)>0){
    tmp[tmp<1] <- 1
    mdl$Dsz <- mdl$Dsz %*% Diagonal(length(tmp),1/tmp)
    tmp <- colSums(mdl$Dsz)
}
mdl$Dex <- mdl$Dsz

mdl$param <- c(mdl$param,'qex_max_default'=Inf)

devtools::load_all("./dynatop")
check_model(mdl)
##mdl$states <- initialise_dynatop(mdl,0.01)
profvis::profvis({mdl$states <- initialise_dynatop(mdl,0.01)})

#dynatop_cp <- compiler::cmpfun(dynatop)

profvis::profvis({tmp <- dynatop(mdl,test_catchment$obs[1:2,],use_states=TRUE)})
#profvis::profvis({tmp <- dynatop_cp(mdl,test_catchment$obs[1:10,],use_states=TRUE)})
tmp <- dynatop(mdl,test_catchment$obs[1:2,],0.1,use_states=TRUE)


## tmp <- lapply(mdl$state,function(x){x})

## it <- 1
## system.time({ tmp <- lapply(mdl$state,function(x,y){
##     x$input[] <- 0
##     x$input['precip'] = y[x$prop$precip_series]
##     x$input['pet'] = y[x$prop$pet_series]
##     return(x)},
##     y=as.matrix(test_catchment$obs)[it,]) })

## hillslope <- mdl$hillslope
## h0 <- create_hillslope()
## prop_names <- intersect(names(h0$prop),names(hillslope))
## prop <- rep(list(NULL),nrow(hillslope))
## for(ii in prop_names){   
##     tmp <- lapply(hillslope[[ii]],
##                   function(x,nm){setNames(list(x),nm)},
##                   nm=ii)
##     prop <- Map(c,prop,tmp)
## }
## prop <- lapply(prop,function(x){list(prop=x)})

## param_names <- intersect(names(h0$param),names(hillslope))
## param <- rep(list(NULL),nrow(hillslope))
## for(ii in param_names){
##     tmp <- lapply(hillslope[[ii]],
##                   function(x,nm){setNames(list(x),nm)},
##                   nm=ii)
##     param <- Map(c,param,tmp)
## }
## param <- lapply(param,function(x){list(param=x)})

## hru <- Map(c,param,prop)

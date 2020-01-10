## quick code to help debug the latest dynamic topmodel
rm(list=ls())

#devtools::load_all("./dynatop")

mdl <- readRDS("./dynatop/data/Swindale_model.rds")
load("./dynatop/data/test_catchment.rda")

for(ii in c("hillslope","channel")){
    mdl[[ii]]$precip_input <- "rain"
    mdl[[ii]]$pet_input <- "pet"
}

mdl$Dsz <- rbind(mdl$Fsz,mdl$Wsz)
mdl$Dsz <- cbind( matrix(0,nrow(mdl$Dsz),nrow(mdl$Dsz)-ncol(mdl$Dsz)),mdl$Dsz )
mdl$Wsz <- mdl$Fsz <- mdl$Wch <- NULL

mdl$hillslope[,'qex_max'] <- 'qex_max_default'
mdl$param <- c(mdl$param,'qex_max_default'=Inf)

mdl$hillslope[,'delta_x'] <- 1
names(mdl$gauge) <- gsub("channel_id","id",names(mdl$gauge))
names(mdl$point_inflow) <- gsub("channel_id","id",names(mdl$point_inflow))


for(ii in 1:ncol(mdl$Dsz)){
    sm <- sum(mdl$Dsz[,ii])
    if( sm > 1 ){
        mdl$Dsz[,ii] <- mdl$Dsz[,ii] / sum(mdl$Dsz[,ii])
    }
    sm <- sum(mdl$Dex[,ii])
    if( sm > 1 ){
        mdl$Dex[,ii] <- mdl$Dex[,ii] / sum(mdl$Dex[,ii])
    }
}



devtools::load_all("./dynatop")
check_model(mdl)
tmp <- dynatop_nwtn(mdl,test_catchment$obs,0.1)
#profvis::profvis({tmp <- dynatop_nwtn(mdl,test_catchment$obs,0.1)})

#Rprof()
#tmp <- dynatop(mdl,test_catchment$obs,0.001)
#Rprof(NULL)
#summaryRprof()

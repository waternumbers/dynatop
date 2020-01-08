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
tmp <- dynatop(mdl,test_catchment$obs)

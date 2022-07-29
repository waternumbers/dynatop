## ##############################
## simple run for debugging
rm(list=ls())
devtools::load_all("../")
data("Swindale");

mdl <- Swindale$model

## temp fix to get correct model form
odfn <- data.frame(name=c("outlet","under_outlet"),id=c(0,0),flux=c("q_sf","q_sz"))
h <- list()
for(ii in 1:nrow(mdl$hru)){
    tmp <- list()
    tmp$id <- as.integer( mdl$hru$id[ii]-1 )
    tmp$states <- setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz"))
    tmp$properties <- c(width =  mdl$hru$area[ii]/mdl$hru$length[ii], area = mdl$hru$area[ii], gradient = mdl$hru$s_bar[ii])

    ## mdl$hru$width[ii], area = mdl$hru$area[ii], gradient = mdl$hru$s_bar[ii])
    ##if( tmp$properties["width"] ==0 ){ tmp$properties["width"] <- 1 }
    tmp$sf <- list(type = "cnstC", #mdl$hru$sf[[ii]]$type,
                   parameters = c("c_sf" = 0.3))
    tmp$sz <- list(type = "exp",
                   parameters = c(t_0=0.0001,m=0.08)) #c(mdl$hru$sz[[ii]]$param, "D" = 0.05))
    if(mdl$hru$is_channel[ii]){
        tmp$sz$parameters["t_0"] <- 1e-60
        tmp$sf$parameters["c_sf"] <- 0.7
    }
    tmp$precip <- mdl$hru$precip[[ii]]
    names(tmp$precip) <- c("name","fraction")
    tmp$pet <- mdl$hru$pet[[ii]]
    names(tmp$pet) <- c("name","fraction")
    tmp$rz <- list(type="orig", parameters = c("s_rzmax" = 0.1))#mdl$hru$s_rzmax[ii]))
    tmp$uz <- list(type="orig", parameters = c("t_d" = 8*60*60)) #mdl$hru$t_d[ii]))
    tmp$sf_flow_direction <- list(id = as.integer(mdl$hru$sf_flow_direction[[ii]]$id - 1),
                                  fraction = mdl$hru$sf_flow_direction[[ii]]$frc)
    tmp$sz_flow_direction <- list(id = as.integer(mdl$hru$sz_flow_direction[[ii]]$id - 1),
                                  fraction = mdl$hru$sz_flow_direction[[ii]]$frc)
    tmp$initialisation <- c("s_rz_0" = 0.98, #mdl$hru$s_rz0[ii], "r_uz_sz_0" = mdl$hru$r_uz_sz0[ii])
                            "r_uz_sz_0" = 1.75e-7) #mdl$hru$r_uz_sz0[ii])
    h[[ii]] <- tmp
}

hh <- h
hh[[1]]$id <- hh[[1]]$id + 0
system.time({
    dt <- dynatop$new(h)
    obs <- Swindale$obs
    dt$add_data(obs)
    dt$initialise()

    tail(dt$get_states())
    head(dt$get_states())
    dt$sim(odfn)
})


x11(); plot(Swindale$obs$flow, ylim=c(0,60))
lines(dt$get_output())

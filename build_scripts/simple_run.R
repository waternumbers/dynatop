## ##############################
## simple run for debugging
rm(list=ls())
devtools::load_all("../")
data("Swindale");

mdl <- Swindale$model

## temp fix to get correct model form
odfn <- data.frame(name="outlet",id=0,flux="q_sf")
h <- list()
for(ii in 1:nrow(mdl$hru)){
    tmp <- list()
    tmp$id <- as.integer( mdl$hru$id[ii]-1 )
    tmp$states <- setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz"))
    tmp$properties <- c(width =  mdl$hru$area[ii]/mdl$hru$length[ii], area = mdl$hru$area[ii], gradient = mdl$hru$s_bar[ii])

    ## mdl$hru$width[ii], area = mdl$hru$area[ii], gradient = mdl$hru$s_bar[ii])
    ##if( tmp$properties["width"] ==0 ){ tmp$properties["width"] <- 1 }
    tmp$sf <- list(type = mdl$hru$sf[[ii]]$type,
                   parameters = c("c_sf" = 0.8))
    tmp$sz <- list(type = "bexp",
                   parameters = c(t_0=exp(7.46),m=0.0063,D=0.5)) #c(mdl$hru$sz[[ii]]$param, "D" = 0.05))
    if(mdl$hru$is_channel[ii]){ tmp$sz$parameters["t_0"] <- 0 }
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

    dt$add_data(Swindale$obs[1:10,])
    dt$initialise()

    tail(dt$get_states())
         dt$sim(odfn)
})

plot(dt$get_output())
plot(Swindale$obs$flow)

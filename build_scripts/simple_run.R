## ##############################
## simple run for debugging
rm(list=ls())
graphics.off()
devtools::load_all()
data("Swindale");

mdl <- Swindale$model

## temp fix to get correct model form
odfn <- data.frame(name=c("outlet","under_outlet"),id=c(0,0),flux=c("q_sf","q_sz"))
h <- list()
for(ii in 1:nrow(mdl$hru)){
    tmp <- list()
    tmp$id <- as.integer( mdl$hru$id[ii]-1 )
    tmp$states <- setNames(as.numeric(rep(NA,6)), c("s_sf","s_rz","s_uz","s_sz","q_sf","q_sz"))
    tmp$properties <- c(width =  mdl$hru$area[ii]/mdl$hru$length[ii], area = mdl$hru$area[ii], gradient = mdl$hru$s_bar[ii], Dx = mdl$hru$length[ii])

    ## mdl$hru$width[ii], area = mdl$hru$area[ii], gradient = mdl$hru$s_bar[ii])
    ##if( tmp$properties["width"] ==0 ){ tmp$properties["width"] <- 1 }
    tmp$sf <- list(type = "cnstCD", #mdl$hru$sf[[ii]]$type,
                   parameters = c("c_sf" = 0.3, "d_sf"=0.0))
    tmp$sz <- list(type = "exp",
                   parameters = c(t_0=0.08,m=0.009,D=5)) #c(mdl$hru$sz[[ii]]$param, "D" = 0.05))
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

##system.time({
dt <- dynatop$new(h)
obs <- Swindale$obs
dt$add_data(obs)
dt$initialise()
s0 <- dt$get_states()
#    print(head(dt$get_states()))
dt$sim(odfn)
sn <- dt$get_states()
##    print(head(dt$get_states()))
##})

x11(); plot(Swindale$obs$flow, ylim=c(0,60))
lines(dt$get_output())



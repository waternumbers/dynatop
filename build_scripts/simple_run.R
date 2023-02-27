## ##############################
## simple run for debugging
rm(list=ls())
graphics.off()
devtools::load_all()
data("Swindale");

hru <- Swindale$model$hru
for(ii in 1:length(hru)){
    if(is.na(hru[[ii]]$class$endNode)){
        ## then HRU is not a channel
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["m"] <- 0.0063
        hru[[ii]]$sz$parameters["t_0"] <- exp(7.46)
        ## unsaturated zone parameters
        hru[[ii]]$uz$parameters["t_d"] <- 8*60*60
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.1
        ## surface parameters
        hru[[ii]]$sf$parameters["c_sf"] <- 0.4
    }else{
        ## then HRU is a channel - set so no subsurface response
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["t_0"] <- 0.001
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.001
        ## surface parameters
        hru[[ii]]$sf$parameters["c_sf"] <- 0.8
    }
    ## initialisation parameters
    hru[[ii]]$initialisation["s_rz_0"] <- 0.98
    hru[[ii]]$initialisation["r_uz_sz_0"] <- 1.755582e-07 ## initial outflow divided by catchment area
}

dt <- dynatop$new(hru,map=system.file("extdata","Swindale.tif",package="dynatop",mustWork=TRUE))

dt$add_data(Swindale$obs)
dt$initialise()

s0 <- dt$get_states()

out_def <- Swindale$model$output_flux
dt$sim(out_def,keep_states=index(Swindale$obs))

sn <- dt$get_states()
colSums(s0)
colSums(sn)
all(s0$s_uz<=s0$s_sz)


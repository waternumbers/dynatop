## ----create_object------------------------------------------------------------
rm(list=ls())
devtools::load_all()
data("Swindale")
swindale_model = Swindale$model
#swindale_model <- readRDS(file.path(".","/build_scripts","vignette_debug","new_model.rds"))
hru <- swindale_model$hru
for(ii in 1:length(hru)){
    hru[[ii]]$sz$type <- "exp"
##    hru[[ii]]$sz$parameters <- hru[[ii]]$sz$parameters,"h_max"=2)
    if(!("endNode" %in% names(hru[[ii]]$class))){
        ## then HRU is not a channel
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["m"] <- 0.0063
        hru[[ii]]$sz$parameters["t_0"] <- 0.1 #exp(7.46) ##0.135 #exp(-1) #5) #)
        ## unsaturated zone parameters
        hru[[ii]]$uz$parameters["t_d"] <- 8*60*60
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.1
        ## surface parameters
        hru[[ii]]$sf$type <- "kin"
        hru[[ii]]$sf$parameters <- c("n"=0.08,"t_raf" = 999.9, "s_raf" = 0)
        ##hru[[ii]]$sf$parameters["v_sf"] <- 0.4
        ## test of raf
        ##hru[[ii]]$sf$parameters["s_raf"] <- 9000
    }else{
        ## then HRU is a channel - set so no subsurface response
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["t_0"] <- 0.000
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.1
        ## surface parameters
        hru[[ii]]$sf$type <- "kin"
        hru[[ii]]$sf$parameters <- c("n"=0.03,"s_raf" = 0,"t_raf" = 999.9)
        ##hru[[ii]]$sf$parameters["v_sf"] <- 0.8
    }
    ## initialisation parameters
    hru[[ii]]$initialisation["s_rz_0"] <- 0.98
    hru[[ii]]$initialisation["r_uz_sz_0"] <- 1.755582e-07 ## initial outflow divided by catchment area
}

## ctch_mdl <- dynatop$new(hru) #,map=swindale_model$map)
## ## ----add_data-----------------------------------------------------------------
## #data("Swindale")
swindale_obs <- Swindale$obs
ctch_mdl <- dynatop$new(hru)
ctch_mdl$add_data(swindale_obs)
ctch_mdl$initialise()
outDefn <- swindale_model$output_flux
outDefn <- rbind(outDefn,outDefn)
outDefn[2,"name"] <- "s_sf_0"
outDefn[2,"flux"] <- "s_sf"
system.time({ ctch_mdl$sim(outDefn) })
o1 <- ctch_mdl$get_states()
y1 <- ctch_mdl$get_output()

o1[1,]

plot(y1$s_sf_0)
plot(swindale_obs$flow);lines(y1$q_sf_0)
lines(yold$q_sf_0,col="red")

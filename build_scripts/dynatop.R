knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(dynatop)
data("Swindale")

names(Swindale)

swindale_model <- Swindale$model
swindale_obs <- Swindale$obs

names(swindale_model)

swindale_model$map <- system.file("extdata","Swindale.tif",package="dynatop",mustWork=TRUE)

head(swindale_obs)

head(swindale_model$hru[[1]]$precip)
head(swindale_model$hru[[1]]$pet)

hru <- swindale_model$hru
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

ctch_mdl <- dynatop$new(hru,map=swindale_model$map)

ctch_mdl$add_data(swindale_obs)



ctch_mdl$initialise()$plot_state("s_sz")

sim1 <- ctch_mdl$sim(swindale_model$output_flux)$get_output()

ctch_mdl$plot_state("s_sz")

sim2 <- ctch_mdl$sim(swindale_model$output_flux)$get_output()
out <- merge( merge(swindale_obs,sim1),sim2)
names(out) <- c(names(swindale_obs),'sim_1','sim_2')
plot(out[,c('flow','sim_1','sim_2')], main="Discharge",ylab="m3/s",legend.loc="topright")

mb <- ctch_mdl$get_mass_errors()
plot( mb[,6] , main="Mass Error", ylab="[m^3]")

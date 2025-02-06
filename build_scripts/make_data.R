###########################################################################
## This converts the data for Swindale into the data object used in the examples.
rm(list=ls())
devtools::load_all()
setwd("./build_scripts")
model <- readRDS("Swindale_exp.rds")

model$hru <- lapply(model$hru,function(h){h$sf$parameters <- c("v_sf_1" = 999, "s_1" = 0, "v_sf_2" = 0.3);h})

qr <- read.csv( "start=2009-11-18_end=2009-11_4_int=0.25-hours_units=mm.hr-1.tsv",sep="\t")
obs <- as.xts(qr[,c("Flow","Rainfall")],order.by= as.POSIXct(qr[,'Date'],format="%d/%m/%Y %H:%M",tz='GMT'))
## According to original notes and code Flow in cumecs and precip in mm/timestep
obs$Rainfall <- obs$Rainfall/1000 # convert to m/timestep
obs$PET <- evap_est(index(obs),0,5/1000) # in m
names(obs) <- c("flow","precip","pet")
Swindale <- list(model=model,obs=obs)
save("Swindale",file="../data/Swindale.rda")

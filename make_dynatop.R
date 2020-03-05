## Script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- './dynatop'
##devtools::load_all(pacPath)
## Rcpp::compileAttributes(pacPath)
devtools::document(pacPath)
devtools::check(pacPath)
tmp <- devtools::build(pacPath)
install.packages(tmp)
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)
#pkgdown::build_article("Time_Delay_Channel_Routing",pacPath)

## This converts the exdata for Swindale into the data object used in the examples.
devtools::load_all(pacPath)
model <- readRDS( system.file("extdata","Swindale.rds",package="dynatop") )

qr <- read.csv( system.file("extdata","start=2009-11-18_end=2009-11_4_int=0.25-hours_units=mm.hr-1.tsv",package="dynatop") ,sep="\t")
obs <- as.xts(qr[,c("Flow","Rainfall")],order.by= as.POSIXct(qr[,'Date'],format="%d/%m/%Y %H:%M",tz='GMT'))
## According to original notes and code Flow in cumecs and precip in mm/timestep
obs$Rainfall <- obs$Rainfall/1000 # convert to m/timestep
obs$PET <- evap_est(index(obs),0,5/1000) # in m

Swindale <- list(model=model,obs=obs)
save("Swindale",file="./dynatop/data/Swindale.rda")


#### Converting the old test_catchments to new ones
rm(list=ls())
load("./test_catchment.rda")

simple_hillslope <- test_catchment$model
simple_hillslope$hillslope$id <- 3:10
simple_hillslope$hillslope$sf_band <- simple_hillslope$hillslope$sz_band <- 1:8
simple_hillslope$hillslope$qsf_max <- "qsf_max_default"
simple_hillslope$hillslope$delta_x <- 1
tmp <- names(simple_hillslope$hillslope)
tmp <- gsub("tex","tsf",tmp)
tmp <- gsub("_input","",tmp)
names(simple_hillslope$hillslope) <- tmp

tmp <- names(simple_hillslope$channel)
tmp <- gsub("tex","tsf",tmp)
tmp <- gsub("_input","",tmp)
names(simple_hillslope$channel) <- tmp
simple_hillslope$channel$sf_band <- simple_hillslope$channel$sz_band <- 9:10

tmp <- names(simple_hillslope$point_inflow)
tmp <- gsub("channel_","",tmp)
names(simple_hillslope$point_inflow) <- tmp

tmp <- names(simple_hillslope$gauge)
tmp <- gsub("channel_","",tmp)
names(simple_hillslope$gauge) <- tmp

simple_hillslope$param['qsf_max_default'] <- Inf

simple_hillslope$Fex <- simple_hillslope$Fsat <- simple_hillslope$Wex <- simple_hillslope$Wsat <- NULL
simple_hillslope$Fsf <- simple_hillslope$Fsz <-
    Matrix::sparseMatrix(i=c(4:10,1,2),j=c(3:10,10),x=c(rep(1,8),0))

obs <- test_catchment$obs
obs$point_1 <- obs$point_2 <- 0


test_catchments <- list(simple_hillslope=simple_hillslope,
                        obs = obs)

save(test_catchments,file="./dynatop/data/test_catchments.rda")





###########################################
## trying to get a decent cimulation for Swindale
rm(list=ls())
pacPath <- './dynatop'
devtools::load_all(pacPath)

hillslope_volume <- function(mdl,ignore_sz=FALSE){
  tmp <- mdl$hillslope
  if(!ignore_sz){
    sum( tmp$area * (tmp$s_sf + tmp$s_rz + tmp$s_uz - tmp$s_sz) )
  }else{
    sum( tmp$area * (tmp$s_sf + tmp$s_rz + tmp$s_uz) )
  }
}

data("Swindale")
#model <- band_model(Swindale$model,5)
model <- Swindale$model
## to fix elsewhere
model$channel$sz_band <- model$channel$sf_band <- 1e4


model$hillslope[,'precip'] <- model$channel[,'precip'] <- "Rainfall"
model$hillslope[,'pet'] <- model$channel[,'pet'] <- "PET"

model$param <- c(q_sfmax_default=0,
                 m_default=0.0063,
                 ln_t0_default=7.46,
                 s_rz0_default=0.98,
                 s_rzmax_default=0.1,
                 v_ch_default=3000,
                 t_d_default=8*60*60,
                 t_sf_default=3.6e+05
                 )

mdlb <- band_model(model,5)
mdlb$hillslope[,'precip'] <- mdlb$channel[,'precip'] <- "Rainfall"
mdlb$hillslope[,'pet'] <- mdlb$channel[,'pet'] <- "PET"

obs <- Swindale$obs
obs$PET[] <- 0
#obs$Rainfall[] <- 0

model0 <- initialise(model,0.000001)
v0 <- hillslope_volume(model0)
qhat <- dynatop(model,obs,0.000001)
vend <- hillslope_volume(qhat$model)
v0 + (sum(obs$Rainfall)*sum(model$hillslope$area)) + (sum(obs$Rainfall)*sum(model$channel$area)) - (sum(qhat$channel_input)*15*60) - vend

mdlb0 <- initialise(mdlb,0.000001)
vb0 <- hillslope_volume(mdlb0)
qhatb <- dynatop(mdlb0,obs,use_states=TRUE)
vbend <- hillslope_volume(qhatb$model)

vb0 + (sum(obs$Rainfall)*sum(mdlb$hillslope$area)) + (sum(obs$Rainfall)*sum(mdlb$channel$area)) - (sum(qhatb$channel_input)*15*60) - vbend

x11();plot(cbind(Swindale$obs[,'Flow'],rowSums(qhatb$channel_input)))

                 rowSums(qhat$channel_input),rowSums(qhatb$channel_input)))



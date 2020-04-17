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

###########################################################################
## This converts the exdata for Swindale into the data object used in the examples.
devtools::load_all(pacPath)
model <- readRDS( system.file("extdata","Swindale_from_r6.rds",package="dynatop") )
for(ii in c("precip","pet")){
    for(jj in c("hillslope","channel")){
        model[[jj]][[ii]] <- switch(ii,precip="Rainfall",pet="PET")
    }
}

qr <- read.csv( system.file("extdata","start=2009-11-18_end=2009-11_4_int=0.25-hours_units=mm.hr-1.tsv",package="dynatop") ,sep="\t")
obs <- as.xts(qr[,c("Flow","Rainfall")],order.by= as.POSIXct(qr[,'Date'],format="%d/%m/%Y %H:%M",tz='GMT'))
## According to original notes and code Flow in cumecs and precip in mm/timestep
obs$Rainfall <- obs$Rainfall/1000 # convert to m/timestep
obs$PET <- evap_est(index(obs),0,5/1000) # in m


Swindale <- list(model=model,obs=obs)
save("Swindale",file="./dynatop/data/Swindale.rda")

devtools::load_all(pacPath); m1 <- dynatop$new(model)$add_data(obs)$initialise(1e-6)$sim_hillslope(mass_check=TRUE)$sim_channel()

devtools::load_all(pacPath)
devtools::load_all(pacPath); m1 <- dynatop$new(model)$add_data(obs)$initialise(1e-6)$sim_hillslope(mass_check=TRUE)
m1$sim_hillslope(mass_check=TRUE)
#profvis::profvis({
    m1 <- dynatop$new(model)
    m1$add_data(obs)
    m1$initialise(1e-6)
    m1$sim_hillslope(mass_check=TRUE)
#})


##################################
## Bodging back together the test catchments
rm(list=ls())
pacPath <- "./dynatop"
devtools::load_all(pacPath)

data("test_catchments")

## change hillslope names
tmp <- names(test_catchments$simple_hillslope$hillslope)
tmpr <- c("qsf_max"="q_sfmax",
          "srz_max"="s_rzmax",
          "srz_0"="s_rz0",
          "td"="t_d",
          "tsf"="t_sf")
for(ii in names(tmpr)){
    tmp <- gsub(ii,tmpr[ii],tmp)
}
names(test_catchments$simple_hillslope$hillslope) <- tmp

chng <- data.frame(t = "hillslope",
                   v="sz_band",
                   c="numeric",
                   stringsAsFactors=FALSE)
chng[1,] <- c("hillslope","sz_band","numeric")
chng[2,] <- c("hillslope","sf_band","numeric")
chng[3,] <- c("channel","sz_band","numeric")
chng[4,] <- c("channel","sf_band","numeric")
chng[5,] <- c("channel","id","integer")
chng[6,] <- c("channel","next_id","integer")
chng[7,] <- c("gauge","id","integer")
chng[8,] <- c("point_inflow","id","integer")

for(ii in 1:nrow(chng)){
    tmp <- test_catchments$simple_hillslope[[ chng[ii,'t'] ]][[ chng[ii,'v'] ]]
    tmp <- switch(chng[ii,'c'],
                  "numeric" = as.numeric(tmp),
                  "integer" = as.integer(tmp),
                  stop("missing type"))
    test_catchments$simple_hillslope[[ chng[ii,'t'] ]][[ chng[ii,'v'] ]] <- tmp
}



check_model(test_catchments$simple_hillslope)


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



model$hillslope[,'precip'] <- model$channel[,'precip'] <- "Rainfall"
model$hillslope[,'pet'] <- model$channel[,'pet'] <- "PET"

mdlb$param <- c(q_sfmax_default=Inf,
                 m_default=0.0063,
                 ln_t0_default=3.46,
                 s_rz0_default=0.98,
                 s_rzmax_default=0.1,
                 v_ch_default=0.8,
                 t_d_default=8,#1,#8*60*60,
                 t_sf_default=3#.6e+05
                 )

mdlb <- band_model(model,5)
mdlb$hillslope[,'precip'] <- mdlb$channel[,'precip'] <- "Rainfall"
mdlb$hillslope[,'pet'] <- mdlb$channel[,'pet'] <- "PET"

obs <- Swindale$obs
#obs$PET[] <- 0
#obs$Rainfall[] <- 0

model0 <- initialise(model,0.000001)
v0 <- hillslope_volume(model0)
qhat3 <- dynatop(model0,obs,sim_time_step=300,use_states=TRUE)
vend <- hillslope_volume(qhat3$model)
v0 + (sum(obs$Rainfall)*sum(model$hillslope$area)) + (sum(obs$Rainfall)*sum(model$channel$area)) - (sum(qhat3$channel_input)*15*60) - vend

mdlb0 <- initialise(mdlb,0.000001)
vb0 <- hillslope_volume(mdlb0)
qhatb <- dynatop(mdlb0,obs,sim_time_step=3600,use_states=TRUE,mass_check=TRUE) #0.000001)
qhatb2 <- dynatop(qhatb$model,obs,sim_time_step=3600,use_states=TRUE,mass_check=TRUE)
vbend <- hillslope_volume(qhatb$model)

vb0 + (sum(obs$Rainfall)*sum(mdlb$hillslope$area)) + (sum(obs$Rainfall)*sum(mdlb$channel$area)) - (sum(qhatb$channel_input)*15*60) - vbend

x11();plot(cbind(Swindale$obs[,'Flow'],rowSums(qhatb$channel_input),
                 rowSums(qhatb2$channel_input)))

chhat <- time_delay_routing(mdlb,qhatb$channel_input)
x11();plot(cbind(Swindale$obs[,'Flow'],rowSums(qhatb$channel_input),chhat))
plot(chhat)


######################

rm(list=ls())
pacPath <- "./dynatop"
devtools::load_all(pacPath)
fn <- list.files(file.path(pacPath,"vignettes"),pattern="*.Rmd$",full.names=TRUE)
rc <- rep(NA,length(fn))
for(ii in 1:length(fn)){
    print(fn[ii])
    rc[ii] <- system.time({ rmarkdown::render(fn[ii],quiet=TRUE) })['elapsed']
    print(rc)
}

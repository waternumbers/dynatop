## Script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- '../'
Rcpp::compileAttributes(pacPath)
devtools::document(pacPath)
devtools::check(pacPath)

## check documentation build
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)
pkgdown::clean_site(pacPath)

## build, populate drat
## linux
dratPath <- "~/Documents/Software/drat"
tmp <- devtools::build(pacPath)
install.packages(tmp)
drat::insertPackage(tmp,dratPath)#,action="prune")

## mac and windows
rhub::validate_email() # for first time that session
pkgName <- sub('\\.tar.gz$', '', basename(tmp)) 
## rhub::platforms()[,1] # lists platforms
mch <- rhub::check(path = tmp,
                   platform = c("macos-highsierra-release-cran","windows-x86_64-release"))

tmp <- paste0(pkgName,".tgz")
ftmp <- file.path("../..",tmp)
download.file(file.path(mch$urls()$artifacts[1],tmp),ftmp)
drat::insertPackage(ftmp,dratPath,action="prune")

tmp <- paste0(pkgName,".zip")
ftmp <- file.path("../..",tmp)
download.file(file.path(mch$urls()$artifacts[2],tmp),ftmp)
drat::insertPackage(ftmp,dratPath,action="prune")

## tidy up drat
drat::pruneRepo(dratPath,pkg="dynatop",remove="git")## this only does source files




###########################################################################
## This converts the data for Swindale into the data object used in the examples.
rm(list=ls())
pacPath <- '..'
devtools::load_all(pacPath)
model <- "Swindale.rds"
model$precip_input$name <- "Rainfall"
model$pet_input$name <- "PET"

qr <- read.csv( "start=2009-11-18_end=2009-11_4_int=0.25-hours_units=mm.hr-1.tsv",sep="\t")
obs <- as.xts(qr[,c("Flow","Rainfall")],order.by= as.POSIXct(qr[,'Date'],format="%d/%m/%Y %H:%M",tz='GMT'))
## According to original notes and code Flow in cumecs and precip in mm/timestep
obs$Rainfall <- obs$Rainfall/1000 # convert to m/timestep
obs$PET <- evap_est(index(obs),0,5/1000) # in m

Swindale <- list(model=model,obs=obs)
save("Swindale",file="../data/Swindale.rda")







## ########################################
## This code uses Swindale and calls all the different transmissivity profiles
rm(list=ls())
devtools::load_all(".."); data("Swindale");

dt <- dynatop$new(Swindale$model)$add_data(Swindale$obs)
dt$initialise(1e-6)$sim_hillslope()
range(dt$get_mass_errors())

mdl <- Swindale$model
mdl$hillslope$c_sz <- 0.5; mdl$hillslope$D <- 5
mdl$options["transmissivity_profile"] <- "constant"
dt <- dynatop$new(mdl)$add_data(Swindale$obs)
dt$initialise(1e-6)$sim_hillslope()
range(rowSums(dt$get_mass_errors()))

mdl <- Swindale$model
mdl$hillslope$D <- 5
mdl$options["transmissivity_profile"] <- "bounded_exponential"
dt <- dynatop$new(mdl)$add_data(Swindale$obs)
dt$initialise(1e-6)$sim_hillslope()
range(rowSums(dt$get_mass_errors()))


############################################################################

profvis::profvis({m1 <- dynatop$new(Swindale$model)$add_data(Swindale$obs)$initialise(1e-6)$sim_hillslope(mass_check=TRUE)$sim_channel(mass_check=TRUE)})


m1 <- dynatop$new(Swindale$model)$add_data(Swindale$obs)$initialise(1e-6)$sim_hillslope(mass_check=TRUE)$sim_channel(mass_check=TRUE)

head( m1$get_channel_inflow() )
head( m1$get_channel_inflow(TRUE) )
m1$plot_channel_inflow()
m1$plot_channel_inflow(TRUE)
head( m1$get_gauge_flow() )
head( m1$get_gauge_flow("channel_18") )
m1$plot_gauge_flow()
m1$plot_gauge_flow("channel_18")
head( m1$get_obs_data() )
tmp <- m1$get_model()
head( m1$get_mass_errors() )
head( m1$get_states() )
head( m1$get_states(TRUE) )
m1$plot_state("s_sf")

## ################################
## This code checks the cpp and r versions of the simulation
rm(list=ls())
devtools::load_all("./dynatop"); data("Swindale")

##profvis::profvis({

idx <- 1:273
m2 <-  dynatop$new(Swindale$model)$add_data(Swindale$obs[idx,])$initialise(1e-6)$sim(mass_check=FALSE,keep_states=zoo::index(Swindale$obs))
m1 <- dynatop$new(Swindale$model)$add_data(Swindale$obs[idx,])$initialise(1e-6)$sim(mass_check=FALSE,use_R=TRUE,keep_states=zoo::index(Swindale$obs))


range(m2$get_states()-m1$get_states())

m1$get_states()[37,]
m2$get_states()[37,]

tail(m2$get_mass_errors())
tail(m1$get_mass_errors())
##})
plot(merge( Swindale$obs[,'Flow'],m1$get_gauge_flow(),m2$get_gauge_flow()))

s2 <- m2$get_states(TRUE)
s1 <- m1$get_states(TRUE)

range(s2[[67]] - s1[[67]])

## ###########################################
## This hits more of the models compoents for testing
rm(list=ls())
devtools::load_all("./dynatop"); data("Swindale"); 
Swindale$model$param <- c(r_sfmax_default=Inf,
                          m_default=0.006, ## 0.05 - 0.6 m
                          ln_t0_default=0.746, ## 0.1 - 8 m2/h
                          s_rz0_default=0.98,
                          s_rzmax_default=0.1,
                          v_ch_default=0.4, ## 1000 -5000 m/h
                          t_d_default=80*60*60, ## 0.1-120 m/h - are these units correct - or is cobceptualisation other way round?
                          c_sf_default=0.4
                          )
## bits that I shouldn't have to do...
Swindale$model$options=c("transmisivity_profile"="exponential","channel_solver"="histogram")

m1 <- dynatop$new(Swindale$model)$add_data(Swindale$obs)$initialise(1e-6)

m1$sim_hillslope()

$sim_channel(mass_check=TRUE)
## $add_data(Swindale$obs)$initialise(1e-6)$sim_hillslope(mass_check=TRUE)$sim_channel(mass_check=TRUE)
m1$plot_gauge_flow()
m1$sim()$plot_gauge_flow()

head(m1$get_mass_errors())


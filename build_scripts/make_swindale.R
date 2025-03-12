## make the Swindale data a set
rm(list=ls())
devtools::load_all()

## ----tempory_dir--------------------------------------------------------------
demo_dir <- tempdir()

## ----initialisation-----------------------------------------------------------
ctch <- dynatopGIS$new(file.path(demo_dir,"example.tif"))

## ----data_files---------------------------------------------------------------
dem_file <- system.file("extdata", "gis","SwindaleDTM40m.tif", package="dynatop", mustWork = TRUE)
channel_file <- system.file("extdata", "gis","SwindaleRiverNetwork.shp", package="dynatop", mustWork = TRUE)


## ----add_catchment------------------------------------------------------------
dem <- terra::rast(dem_file)
dem <- terra::extend(dem,1) ## pad with NA values
catchment_outline <- terra::ifel(is.finite(dem),1,NA)
ctch$add_catchment(catchment_outline)


## ----add_dem------------------------------------------------------------------
ctch$add_dem(dem)


## ----channel_current----------------------------------------------------------
sp_lines <- terra::vect(channel_file)
head(sp_lines)


## ----channel_properties-------------------------------------------------------
property_names <- c(name="identifier",
                    endNode="endNode",
                    startNode="startNode",
                    length="length")
chn <- convert_channel(sp_lines,property_names)


## ----add_channel--------------------------------------------------------------
ctch$add_channel(chn)

## ----sink_fill----------------------------------------------------------------
ctch$sink_fill()

ctch$compute_properties()

## ----atb_split----------------------------------------------------------------
ctch$add_layer( ctch$classify("atb_20","atb",cuts=20) )
ctch$plot_layer("atb_20")

## ----model_atb_split----------------------------------------------------------
ctch$create_model(file.path(demo_dir,"new_model"),"atb_20")
##ctch$create_model(file.path(".","/build_scripts","vignette_debug","new_model"),"atb_20")


model <- readRDS(file.path(demo_dir,"new_model.rds"))

qr <- read.csv( "./build_scripts/start=2009-11-18_end=2009-11_4_int=0.25-hours_units=mm.hr-1.tsv",sep="\t")
obs <- as.xts(qr[,c("Flow","Rainfall")],order.by= as.POSIXct(qr[,'Date'],format="%d/%m/%Y %H:%M",tz='GMT'))
## According to original notes and code Flow in cumecs and precip in mm/timestep
obs$Rainfall <- obs$Rainfall/1000 # convert to m/timestep
obs$PET <- evap_est(index(obs),0,5/1000) # in m
names(obs) <- c("flow","precip","pet")
Swindale <- list(model=model,obs=obs)

save("Swindale",file="./data/Swindale.rda")
file.copy(file.path(demo_dir,"new_model.tif"), "./inst/extdata/mdl/Swindale.tif",overwrite=TRUE)
file.copy(file.path(demo_dir,"example.geojson"), "./inst/extdata/mdl/Swindale.geojson",overwrite=TRUE)

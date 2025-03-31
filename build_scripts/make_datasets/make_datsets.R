## make the data sets
rm(list=ls())
graphics.off()

devtools::load_all()
setwd("./build_scripts/make_datasets")

## make the model
ctch <- dynatopGIS$new("example.tif")
dem_file <- system.file("extdata", "gis","SwindaleDTM40m.tif", package="dynatop", mustWork = TRUE)
channel_file <- system.file("extdata", "gis","SwindaleRiverNetwork.shp", package="dynatop", mustWork = TRUE)

dem <- terra::rast(dem_file)
dem <- terra::extend(dem,1) ## pad with NA values
catchment_outline <- terra::ifel(is.finite(dem),1,NA)
ctch$add_catchment(catchment_outline)
ctch$add_dem(dem)
sp_lines <- terra::vect(channel_file)
property_names <- c(name="identifier",
                    endNode="endNode",
                    startNode="startNode",
                    length="length")
chn <- convert_channel(sp_lines,property_names)
ctch$add_channel(chn)
ctch$sink_fill()
ctch$compute_properties()
tmp <- ctch$get_layer("filled_dem")
tmp <- terra::ifel(tmp<=500,0,1)
ctch$add_layer(tmp, "greater_500")
ctch$add_layer( ctch$classify("atb_20","atb",cuts=20) )
ctch$add_layer( ctch$combine_classes("atb_20_band",c("atb_20","band")) )
ctch$create_model("new_model","atb_20")

## create swindale data object
model <- readRDS("new_model.rds")
qr <- read.csv( "start=2009-11-18_end=2009-11_4_int=0.25-hours_units=mm.hr-1.tsv",sep="\t")
obs <- as.xts(qr[,c("Flow","Rainfall")],order.by= as.POSIXct(qr[,'Date'],format="%d/%m/%Y %H:%M",tz='GMT'))
## According to original notes and code Flow in cumecs and precip in mm/timestep
obs$Rainfall <- obs$Rainfall/1000 # convert to m/timestep
obs$PET <- evap_est(index(obs),0,5/1000) # in m
names(obs) <- c("flow","precip","pet")
Swindale <- list(model=model,obs=obs)
save(Swindale, file="Swindale.rda")
tools::resaveRdaFiles("Swindale.rda")


## copy files across to correct locations
## copy the test outputs
ov = TRUE
file.copy("new_model.rds","../../inst/tinytest/test_output/",overwrite=ov)
file.copy("new_model.tif","../../inst/tinytest/test_output/",overwrite=ov)
file.copy("example.tif","../../inst/tinytest/test_output/",overwrite=ov)
file.copy("example.geojson","../../inst/tinytest/test_output/",overwrite=ov)
file.copy("Swindale.rda","../../data/",overwrite=ov)
file.copy("new_model.tif","../../inst/extdata/mdl/Swindale.tif",overwrite=ov)
file.copy("example.geojson","../../inst/extdata/mdl/Swindale.geojson",overwrite=ov)

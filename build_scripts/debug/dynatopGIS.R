rm(list=ls())
graphics.off()
devtools::load_all(".")

## ----tempory_dir-------------------------------------------------------------------------------------------------
demo_dir <- "./build_scripts/debug/demo"
unlink(demo_dir, recursive=TRUE)
dir.create(demo_dir)

## ----data_files--------------------------------------------------------------------------------------------------
dem_file <- system.file("extdata", "GIS", "SwindaleDTM.tif", package="dynatop", mustWork = TRUE)
channel_file <- system.file("extdata", "GIS", "SwindaleRiverNetwork.gpkg", package="dynatop", mustWork = TRUE)
catchment_outline <- system.file("extdata", "GIS", "SwindaleBoundary.gpkg", package="dynatop", mustWork = TRUE)

## ----initialisation----------------------------------------------------------------------------------------------
ctch <- dynatopGIS$new(demo_dir)

## ----add_catchment-----------------------------------------------------------------------------------------------
ctch$add_catchment(catchment_outline,dem_file)


## ----channel_current---------------------------------------------------------------------------------------------
channel_lines <- terra::vect(channel_file)
head(channel_lines)


## ----channel_properties------------------------------------------------------------------------------------------
property_names <- c(uid = "identifier",
                    endNode = "endNode",
                    startNode = "startNode",
                    name="name1")
chn <- convert_channel(channel_lines,property_names)

## ----check channel------------

print(check_channel(chn,outlets = "F6D9CBCC-436F-46E0-A631-1F7A5F007FBF"))



## ----add_channel-------------------------------------------------------------------------------------------------
ctch$add_channel(chn)


## ----list_layers-------------------------------------------------------------------------------------------------
ctch$get_layer()


## ----plot--------------------------------------------------------------------------------------------------------
ctch$plot_layer("dem", add_channel=TRUE)


## ----get_layer---------------------------------------------------------------------------------------------------
ctch$get_layer("dem")


## ----sink_fill---------------------------------------------------------------------------------------------------
ctch$sink_fill()


terra::plot( ctch$get_layer('filled_dem') - ctch$get_layer('dem'),
            main="Changes to height")


## ----band--------------------------------------------------------------------------------------------------------
ctch$plot_layer("band")
ctch$plot_layer("hru")

#devtools::load_all(".")
#ctch <- dynatopGIS$new(demo_file)
ctch$create_model("./inst/extdata/mdl/SwindaleModel")


## ----extract_filled----------------------------------------------------------------------------------------------
tmp <- ctch$get_layer("filled_dem")


## ----height layer------------------------------------------------------------------------------------------------
## T
tmp <- terra::ifel(tmp<=500,0,1)


## ----add_height_layer--------------------------------------------------------------------------------------------
ctch$add_layer(tmp, "greater_500")
ctch$get_layer()


rm(list=ls())
devtools::load_all(".")


## ----tempory_dir-------------------------------------------------------------------------------------------------
demo_file <- "./build_scripts/debug/dygis.tif"
unlink(demo_file)


## ----initialisation----------------------------------------------------------------------------------------------
ctch <- dynatopGIS$new(demo_file)


## ----data_files--------------------------------------------------------------------------------------------------
dem_file <- system.file("extdata", "gis", "SwindaleDTM40m.tif", package="dynatop", mustWork = TRUE)
channel_file <- system.file("extdata", "gis", "SwindaleRiverNetwork.shp", package="dynatop", mustWork = TRUE)


## ----add_catchment-----------------------------------------------------------------------------------------------
dem <- terra::rast(dem_file)
dem <- terra::extend(dem,1) ## pad with NA values
catchment_outline <- terra::as.polygons(terra::ifel(is.finite(dem),1,NA),dissolve=TRUE)
ctch$add_catchment(catchment_outline,dem)


## ----channel_current---------------------------------------------------------------------------------------------
sp_lines <- terra::vect(channel_file)
head(sp_lines)


## ----channel_properties------------------------------------------------------------------------------------------
property_names <- c(endNode="endNode",
                    startNode="startNode")
chn <- convert_channel(sp_lines,property_names)


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


ctch$create_model("test_model")


## ----calc_atb----------------------------------------------------------------------------------------------------
ctch$compute_properties()


## ----plot_atb----------------------------------------------------------------------------------------------------
## plot of topographic index (log(a/tan b))
ctch$plot_layer('atb')


## ----flow_length-------------------------------------------------------------------------------------------------
ctch$compute_flow_lengths(flow_routing="shortest")


## ----flow_length_plot--------------------------------------------------------------------------------------------
ctch$get_layer()
ctch$plot_layer("shortest_flow_length")


## ----extract_filled----------------------------------------------------------------------------------------------
tmp <- ctch$get_layer("filled_dem")


## ----height layer------------------------------------------------------------------------------------------------
## T
tmp <- terra::ifel(tmp<=500,NA,-999)


## ----add_height_layer--------------------------------------------------------------------------------------------
ctch$add_layer(tmp, "greater_500")
ctch$get_layer()


## ----atb_split---------------------------------------------------------------------------------------------------
ctch$classify("atb_20","atb",cuts=20)
ctch$plot_layer("atb_20")


## ----atb_splt_get_class------------------------------------------------------------------------------------------
ctch$get_method("atb_20")


## ----atb_20_band-------------------------------------------------------------------------------------------------
ctch$combine_classes("atb_20_band",c("atb_20","band"))
ctch$plot_layer("atb_20_band")


## ----atb_20_band_burn--------------------------------------------------------------------------------------------
ctch$combine_classes("atb_20_band_500",pairs=c("atb_20","band"),burns="greater_500")
ctch$plot_layer("atb_20_band_500")


## ----see_class---------------------------------------------------------------------------------------------------
head( ctch$get_method("atb_20_band_500")$groups )


## ----model_atb_split---------------------------------------------------------------------------------------------
ctch$create_model(file.path(demo_dir,"new_model"),"atb_20")


## ----model files-------------------------------------------------------------------------------------------------
list.files(demo_dir,pattern="new_model*")


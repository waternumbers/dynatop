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

bnd <- ctch$get_layer("band")
chn <- ctch$get_layer("channel")
hru <- ctch$get_layer("hru")

range(bnd)
range(chn$band)
chn$band

#devtools::load_all(".")
#ctch <- dynatopGIS$new(demo_file)
ctch$create_model("./inst/extdata/mdl/SwindaleModel")

## ----compute properties
devtools::load_all(".")
ctch <- dynatopGIS$new(demo_dir)
ctch$compute_properties()



## ----extract_filled----------------------------------------------------------------------------------------------
tmp <- ctch$get_layer("filled_dem")
tmp[] <- 40*40
names(tmp) <- "area"
ctch$add_layer(tmp, "area")
ctch$accumulate_layer("area")
ua <- ctch$get_layer("upslope_area")
aa <- ctch$get_layer("acc_area")
all( aa - ua == 0 )


## ###############################################################
## see what happens if we use lines for computation
rm(list=ls())
graphics.off()
devtools::load_all(".")

## ----tempory_dir-------------------------------------------------------------------------------------------------
demo_dir <- "./build_scripts/debug/demo"
unlink(demo_dir, recursive=TRUE)
dir.create(demo_dir)

ctch <- dynatopGIS$new(demo_dir)
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
chn <- channel_lines[,property_names]
names(chn) <- names(property_names)
chn$width <- 2
chn$slope <- 0.001
chn$depth <- 1
## ----check channel------------
print(check_channel(chn,outlets = "F6D9CBCC-436F-46E0-A631-1F7A5F007FBF", chn_is_lines=TRUE))

## ----add_channel-------------------------------------------------------------------------------------------------
ctch$add_channel(chn,chn_is_lines=TRUE)


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

bnd <- ctch$get_layer("band")
chn <- ctch$get_layer("channel")
hru <- ctch$get_layer("hru")

range(bnd)
range(chn$band)
chn$band

#devtools::load_all(".")
#ctch <- dynatopGIS$new(demo_file)
ctch$create_model("./inst/extdata/mdl/SwindaleModel")

## ----compute properties
ctch$compute_properties()



## ----extract_filled----------------------------------------------------------------------------------------------
tmp <- ctch$get_layer("filled_dem")
tmp[] <- 40*40
names(tmp) <- "area"
ctch$add_layer(tmp, "area")
ctch$accumulate_layer("area")
ua <- ctch$get_layer("upslope_area")
aa <- ctch$get_layer("acc_area")
tmp <- ctch$get_layer("channel")
all( tmp$upslope_area - tmp$acc_area == 0 )
abs( max(tmp$upslope_area) / terra::global(is.finite(ctch$get_layer("catchment")),sum) - 1600 ) < 1e-10


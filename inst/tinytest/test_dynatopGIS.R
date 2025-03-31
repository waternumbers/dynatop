##setwd("./inst/tinytest") ## comment out
library(tinytest)
##devtools::load_all() ## comment out
library(dynatop)

demo_dir <- tempfile("dygis")
on.exit( unlink(demo_dir) )
dir.create(demo_dir)

## load the comparative output
expect_silent({
    brk <- terra::rast("./test_output/example.tif")
    chn <- terra::vect("./test_output/example.geojson")
})

## test creation
expect_silent({
    ctch <- dynatopGIS$new(file.path(demo_dir,"demo.tif"))
    dem_file <- system.file("extdata", "gis", "SwindaleDTM40m.tif", package="dynatop", mustWork = TRUE)
    channel_file <- system.file("extdata", "gis", "SwindaleRiverNetwork.shp", package="dynatop", mustWork = TRUE)
})

## test adding a catchment
expect_silent({
    dem <- terra::rast(dem_file)
    dem <- terra::extend(dem,1)
    catchment_outline <- terra::ifel(is.finite(dem),1,NA)
    ctch$add_catchment(catchment_outline)
})
## test adding dem
expect_silent({
    ctch$add_dem(dem)
})

expect_true( terra::identical(ctch$get_layer("dem"), brk[["dem"]]) )

## test adding channel
expect_silent({
    suppressWarnings({ sp_lines <- terra::vect(channel_file) })
    property_names <- c(name="identifier",
                        endNode="endNode",
                        startNode="startNode",
                        length="length")
    suppressWarnings({ chn <- convert_channel(sp_lines,property_names) })
    ctch$add_channel(chn)
})

expect_true( terra::identical(ctch$get_layer("channel"), brk[["channel"]]) )

## terra identical and compareGeom don't appear to work for SpatVector objects
## expect_silent({
##     tmp <- ctch$get_layer("channel_vect")
##     tmp$slope <- NULL
##     ttmp <- terra::vect("./test_output/demo/channel.shp")
##     ttmp$to_keep <- NULL
## })
## expect_true( terra::identical( tmp, ttmp ) )

## test dem filling
expect_silent({ ctch$sink_fill() })
expect_true( terra::identical(ctch$get_layer("filled_dem"), brk[["filled_dem"]]) )



## Check compute properties
expect_silent({ ctch$compute_properties() })
expect_true( terra::identical(ctch$get_layer("band"), brk[["band"]]) )
expect_true( terra::identical(ctch$get_layer("gradient"), brk[["gradient"]]) )
expect_true( terra::identical(ctch$get_layer("atb"), brk[["atb"]]) )
expect_true( terra::identical(ctch$get_layer("upslope_area"), brk[["upslope_area"]]) )


## test adding a layer
expect_silent({ 
    tmp <- ctch$get_layer("filled_dem")
    tmp <- terra::ifel(tmp<=500,0,1)
    names(tmp) <- "greater_500"
    ctch$add_layer(tmp, "greater_500")
})
expect_true( terra::identical(ctch$get_layer("greater_500"), tmp) )

## test a classification
expect_silent({ ctch$add_layer( ctch$classify("atb_20","atb",cuts=20) )})
expect_true( terra::identical(ctch$get_layer("atb_20"), brk[["atb_20"]]) )

## test combining classes (simple)
expect_silent({ ctch$add_layer( ctch$combine_classes("atb_20_band",c("atb_20","band")) ) })
expect_true( terra::identical(ctch$get_layer("atb_20_band"), brk[["atb_20_band"]]) )

## test combining classes (complex)
expect_silent({ ctch$add_layer( ctch$combine_classes("atb_20_band_500",pairs=c("atb_20","band"),burns="greater_500") ) })
## TODO - generate expect_true( terra::identical(ctch$get_layer("atb_20_band_500"), terra::rast("./test_output/demo/atb_20_band_500.tif")) )

expect_silent({ ctch$create_model(file.path(demo_dir,"new_model"),"atb_20") })
expect_silent({
    tmp <- readRDS( file.path(demo_dir,"new_model.rds") )
    ttmp <- readRDS( "./test_output/new_model.rds")
    tmp$map <- ttmp$map <- "no mapfor testing"
})
expect_true( terra::identical(terra::rast( file.path(demo_dir,"new_model.tif") ),
                              terra::rast("./test_output/new_model.tif")) )
expect_identical(tmp, ttmp)

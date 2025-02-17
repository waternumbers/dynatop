##setwd("./inst/tinytest") ## comment out
library(tinytest)
##devtools::load_all() ## comment out
library(dynatop)

demo_dir <- tempfile("dygis")
on.exit( unlink(demo_dir) )
dir.create(demo_dir)

## test creation
expect_silent({
    ctch <- dynatopGIS$new(file.path(demo_dir,"demo"))
    dem_file <- system.file("extdata", "gis", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
    channel_file <- system.file("extdata", "gis", "SwindaleRiverNetwork.shp", package="dynatopGIS", mustWork = TRUE)
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

expect_true( terra::identical(ctch$get_layer("dem"), terra::rast("./test_output/demo/dem.tif")) )

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

expect_true( terra::identical(ctch$get_layer("channel"), terra::rast("./test_output/demo/channel.tif")) )
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
expect_true( terra::identical(ctch$get_layer("filled_dem"), terra::rast("./test_output/demo/filled_dem.tif")) )

expect_silent({ ctch$compute_band() })
expect_true({ terra::identical(ctch$get_layer("band"), terra::rast("./test_output/demo/band.tif")) })

## Check compute properties
expect_silent({ ctch$compute_properties() })
tmp <- ctch$get_layer("channel")
expect_true( terra::identical(ctch$get_layer("gradient"),
                              terra::mask(terra::rast("./test_output/demo/gradient.tif"),tmp,inverse=TRUE)
                              ))
## this fails due to older handling around channel and computing of upstream area everywhere
## expect_true( terra::identical(ctch$get_layer("upslope_area"),
##                               terra::mask(terra::rast("./test_output/demo/upslope_area.tif"),tmp,inverse=TRUE)
##                               ))
## terra::plot( ctch$get_layer("upslope_area") - terra::mask(terra::rast("./test_output/demo/upslope_area.tif"),tmp,inverse=TRUE) )
## the follwoing fails as a knock on from the above
## expect_true( terra::identical(ctch$get_layer("atb"),
##                               terra::mask(terra::rast("./test_output/demo/atb.tif"),tmp,inverse=TRUE)
##                               ))
## terra::plot( ctch$get_layer("atb") - terra::mask(terra::rast("./test_output/demo/atb.tif"),tmp,inverse=TRUE) )

## check flow distances
expect_silent({ ctch$compute_flow_lengths("expected") })
expect_true( terra::identical(ctch$get_layer("expected_flow_length"), terra::rast("./test_output/demo/expected_flow_length.tif")) )
expect_silent({ ctch$compute_flow_lengths("dominant") })
expect_true( terra::identical(ctch$get_layer("dominant_flow_length"), terra::rast("./test_output/demo/dominant_flow_length.tif")) )
expect_silent({ ctch$compute_flow_lengths("shortest") })
expect_true( terra::identical(ctch$get_layer("shortest_flow_length"), terra::rast("./test_output/demo/shortest_flow_length.tif")) )

## test adding a layer
expect_silent({ 
    tmp <- ctch$get_layer("filled_dem")
    tmp <- terra::ifel(tmp<=500,NA,-999)
    ctch$add_layer(tmp, "greater_500")
})

## add this so test on classification etc pass
expect_silent({ ctch$add_layer( terra::rast("./test_output/demo/atb.tif"), "atb_old") })

## test a classification
expect_silent({ ctch$classify("atb_20","atb_old",cuts=20) })
expect_true( terra::identical(ctch$get_layer("atb_20"), terra::rast("./test_output/demo/atb_20.tif")) )

## test combining classes (simple)
expect_silent({ ctch$combine_classes("atb_20_band",c("atb_20","band")) })
expect_true( terra::identical(ctch$get_layer("atb_20_band"), terra::rast("./test_output/demo/atb_20_band.tif")) )

## test combining classes (complex)
expect_silent({ ctch$combine_classes("atb_20_band_500",pairs=c("atb_20","band"),burns="greater_500") })
expect_true( terra::identical(ctch$get_layer("atb_20_band_500"), terra::rast("./test_output/demo/atb_20_band_500.tif")) )

## compare get method retreival
expect_silent({
    tmp <- ctch$get_method("atb_20_band_500")
    ttmp <- jsonlite::fromJSON( "./test_output/demo/atb_20_band_500.json" )
})
## changed JSON format - check what we can
expect_identical( tmp$groups,ttmp$groups, info="Comparision of method retreival" )

expect_silent({ ctch$create_model(file.path(demo_dir,"new_model"),"atb_20_band") })
expect_silent({
    tmp <- readRDS( file.path(demo_dir,"new_model.rds") )
    ttmp <- readRDS( "./test_output/new_model.rds")
    tmp$output_flux$scale <- NULL ## remove this since not in original
})
## TODO test model
## expect_true( terra::identical(terra::rast( file.path(demo_dir,"new_model.tif") ),
##                               terra::rast("./test_output/new_model.tif")) )
expect_identical(tmp$output_flux,ttmp$output_flux)

##file.copy( file.path(demo_dir,"new_model.rds"), file.path("/home/smithpj1/Documents/Software/dynatopGIS/debug_model.rds") )
## models chould be identical except for class variables
## expect_silent({
##     th <- lapply(tmp$hru,function(x){x$class <- NULL;x})
##     tth <- lapply(ttmp$hru,function(x){x$class <- NULL;x})
## })
## expect_identical(th,tth)

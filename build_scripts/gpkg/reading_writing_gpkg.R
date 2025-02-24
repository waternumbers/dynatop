## reading a wring gpkg files
rm(list=ls())
library(terra)

rst <- rast("./inst/extdata/gis/SwindaleDTM40m.tif")
vec <- vect("./inst/extdata/gis/SwindaleRiverNetwork.shp")


writeRaster(rst,"test.gpkg",gdal = c("RASTER_TABLE=dem","FIELD_NAME=dem"),overwrite=TRUE)

names(rst) <- "egg"
writeRaster(
  rst,
  "test.gpkg",
  filetype = "GPKG",
  gdal = c("APPEND_SUBDATASET=YES", "RASTER_TABLE=dem2", "FIELD_NAME=dem2")
)

writeRaster(
  rst,
  "test.gpkg",
  filetype = "GPKG",
  gdal = c("APPEND_SUBDATASET=YES", "RASTER_TABLE=dem", "FIELD_NAME=egg")
)


writeVector(
  vec,
  "test.gpkg",
  layer="channel",
  filetype = "GPKG",
  insert=TRUE
)


names(rast("test.gpkg"))
names(vect("test.gpkg"))

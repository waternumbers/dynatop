rm(list=ls())
library(terra)

dem <- rast("./inst/extdata/gis/SwindaleDTM40m.tif")
chn <- vect("./inst/extdata/gis/SwindaleRiverNetwork.shp")

devtools::load_all()

prop_names=c(name = "identifier",
             length = "length",
             startNode = "startNode",
             endNode = "endNode")
tmp <- convert_channel(chn,prop_names)


chn_frac <- rasterize(tmp,dem,background=0,cover=TRUE)
chn_frac <- mask(chn_frac,dem)

#prc <- mask(prc,hru)
#names(prc) <- "prc_con"

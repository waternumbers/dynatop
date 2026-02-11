## see about snap to grid
rm(list=ls())
library(terra)
library(sf)
library(lwgeom)

chn <- st_read("./inst/extdata/GIS/SwindaleRiverNetwork.gpkg")
dem <- rast("./inst/extdata/GIS/SwindaleDTM.tif")

## get reslution from DEM
rs <- res(dem)
orgn <- as.vector(ext(dem))[c("xmin","ymin")]
orgn <- orgn - rs/2
rs <- c(rs,0,0)
orgn <- c(orgn,0,0)

chn <- st_segmentize(chn,rs[1])
snp <- st_snap_to_grid(chn,c(rs,0,0),c(orgn,0,0))
st_write(snp,"snapped_channel.gpkg",append=FALSE)


idx <- st_is_empty(snp)
for(ii in which(idx)){
    snp$endNode[ snp$endNode == snp$startNode[ii] ] <- snp$endNode[ii]
}
snp <- snp[!idx,]



##tmp <- extractAlong(dem,snp,cells=TRUE)

plot(st_geometry(chn))
plot(st_geometry(snp),col="red",add=TRUE)

plot(dem)
plot(st_geometry(snp),col="red",add=TRUE)



tmp <- as.data.frame(st_coordinates(snp))

tmp <- split(tmp[,c("X","Y")], tmp[,setdiff(names(tmp),c("X","Y"))])

f_simplify <- function(xy){
    idx <- cellFromXY(dem,xy)
    jdx <- idx
    dup <- jdx[duplicated(jdx)]
    while( length(dup) > 0 ){
        tmp <- sapply(dup,function(x){ range(which(idx==x)) })
        tmp <- tmp[,which.min(tmp[1,])]
        jdx <- jdx[ -(tmp[1]:(tmp[2]-1))]
        dup <- jdx[duplicated(jdx)]
    }
    return(jdx)
}

out <- lapply(tmp,f_simplify)

dout <- lapply(out,diff)
all( unique(abs(unlist(dout))) %in% c(1,ncol(dem),ncol(dem)-1,ncol(dem)+1) )

## check connectivity - assumes snp is correct
se <- sapply(out,function(x){c(x[1],tail(x,1))})
srt <- !(se[1,] %in% se[2,])
nd <- !(se[2,] %in% se[1,])
all(nd == !(snp$endNode %in% snp$startNode))
all(srt == !(snp$startNode %in% snp$endNode))

for(ii in which(nd)){ out[[ii]] <- c(out[[ii]],-999) }

chn_next <- dem
chn_next[] <- NA
for(ii in 1:length(out)){
    chn_next[head(out[[ii]],-1)] <- tail(out[[ii]],-1)
}

## process
## dem - give each cell a number
## burn in snapped channel - give neighbouring cell number - store width & depth
## add water body - renumber cells + alter neighbouring cells in channel - crate wb layer of 1
## sink fille - make sure of drainge
## compute band - renumber the cells etc to hru number
## add other classes (lu cannot overwrite lakes)
## write out model inc edges which match gradient / channel flow direction
##


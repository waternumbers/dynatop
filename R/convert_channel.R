## to fix
## warningnote when gm[] <- NULL
## documentation
## crop channels when merging in waterbodies - check same number  of inputs and outputs (unless wb added at top)


#' Function for assisting in the conversion of object to be suitable channel inputs to a dynatopGIS object
#'
#' @description Converts SpatVect of the channel network into SpatVect containing polgons suitable for dynatopGIS
#' @param chn a SpatVect object or a file which can read by terra::vect to create one
#' @param property_names a named vector of containing the columns of existing data properties required in the final SpatialPolygonsDataFrame
#' @param defaults default values used to replace missing widths, slopes and depths
#'
#' @return A SpatVect containing polygons of the channel network, with at least the following properties: uid, width, depth, slope, startNode and endNode.
#'
#' @details The processing follows the follwoing sequence:
#'   - The `property_names` input is used to rename the data associated with spatial objects. Remaining data is dropped
#'   - Varaibles are converted to the expected type then checked for missing values; in doing this
#'       - Missing `slope`, `width` and `depth`  values are populated by the default and a warning issued
#'   - If the spatial objects are lines
#'      - The spatial objects are buffered using the half `width`
#'
#' @examples
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' vect_lines <- terra::vect(channel_file)
#' property_names <- c(uid="identifier",endNode="endNode",startNode="startNode")
#' chn <- convert_channel(vect_lines,property_names)
#' @export
convert_channel <- function(chn,
                            property_names=c(uid = "uid",
                                             startNode = "startNode",
                                             endNode = "endNode",
                                             width = "width",
                                             slope = "slope",
                                             depth = "depth"),
                            defaults = c("width"=2,"slope"=0.001,"depth"=1)
                            ){

    ## read in the chn sp object is a character sting
    if(is.character(chn)){
        if(file.exists(chn)){
            chn <- terra::vect(chn)
        }else{
            stop("chn is a character string but the file specified does not exist")
        }
    }

    ## check the SptVect object and property names
    stopifnot("The channel network does not have a line or polygon geometry (even when read in)" =
                  terra::geomtype(chn) %in% c("lines","polygons"),
              "A field given in chn_property_names does not exist" = all( property_names %in% names(chn) ),
              "A new field name is duplicated" = !any(duplicated(names(property_names)))
              )
    is_polygon <- terra::geomtype(chn)=="polygons"

    ## checks the on defults
    stopifnot("defaults should be numeric" = is.numeric(defaults),
              "defaults should be finite" = all(is.finite(defaults)),
              "default values for width, slope and depth required" = all( c("width","slope","depth") %in% names(defaults) ))


    ## mutate the names so that they match those on the property_names
    chn <- chn[,property_names]
    names(chn) <- names(property_names)

    ## check required columns with defualts exist
    chn[,setdiff(c("width","depth","slope"),names(chn))] <- NA

    ## some type conversions
    chn$uid <- as.character(chn$uid)
    chn$width <- as.numeric(chn$width)
    chn$slope <- as.numeric(chn$slope)
    chn$startNode <- as.character(chn$startNode)
    chn$endNode <- as.character(chn$endNode)
    chn$depth <- as.numeric(chn$depth)

    ## populate the default missing values
    for(ii in c("slope","width","depth")){
        idx <- is.na(chn[[ii]])
        if(any(idx)){
            warning( paste("Replacing missing",ii,"values with default") )
            chn[[ii]][idx] <- defaults[ii]
        }
    }

    ## make into polygom if required
    if(!is_polygon){
        chn <- terra::buffer(chn, width=chn$width/2)
    }

    ## final check it is a channel...
    check_channel(chn)

    return(chn)
}

#' Trim channel objects
#'
#' @description Trim channel objects to those draining to named outlets
#'
#' @param chn the channel object
#' @param outlets names of the outlets to consider
#' @param removed should the removed objects be returned
#'
#' @details Starting with the named outlet objects the channel network is passed up using the
#' connectivity given by the start and end node, only those flaged as connected are returned.
#' @export
trim_channel <- function(chn,outlets,removed=FALSE){
       stop("Not updated")

    ## check chn is a channel
    if(is.character(chn)){ x <- terra::vect(chn) }
    check_channel(chn)

    ## check outlets are in the channel
    stopifnot("Not all outlets are in the channel object" = all( outlets %in% chn$name ))

    ## set up initial conditions
    idx <- chn$name %in% outlets
    not_eval <- !idx
    in_network <- idx

    ## for speed
    sN <- chn$startNode
    eN <- chn$endNode

    while( any(idx) ){
        idx <- (eN %in% sN[idx]) & not_eval
        not_eval[idx] <- FALSE
        in_network[idx] <- TRUE
    }

    if( removed ){ in_network <- !in_network }

    chn <- chn[in_network,]

    if( !removed ){ check_channel(chn,outlets) }

    return(chn)
}


#' Merge two channel objects
#'
#' @description Merge on channel object into another by burning them in
#'
#' @param x the channel object to burn into
#' @param y the channel object to be burnt in
#' @param outlets names of the channel lengths that are outlets, must be in x
#' @param verbose show progress
#'
#' @details Elements in x are either cropped, or fully removed if they lie under y. Attempts are made to ensure that the start and end Nodes are suitably replaced
#' @export
merge_channels <- function(x,y,outlets=NULL,verbose=FALSE,check_connectivity=TRUE){

    ## check x is a channel
    if(is.character(x)){ x <- terra::vect(x) }
    check_channel(x,outlets)

    ## check y is a channel
    if(is.character(y)){ x <- terra::vect(y) }
    check_channel(y)

    ## check names and node names are unique and projection
    stopifnot("names must be unique" = !any(x$uid %in% y$uid),
              "node names must be unique" = !any(c(x$startNode,x$endNode) %in% c(y$startNode,y$endNode)), ## this is prtty ineffiecent?
              "Projection of channel objects does not match" =  terra::crs(x, proj=TRUE)==terra::crs(y, proj=TRUE)
              )

    if(verbose){
        pb <- txtProgressBar(min = 0, max = nrow(y), initial = 0, char = "=",
                             width = NA, "progress...", "whoop", style = 3, file = "")
        on.exit( close(pb) )
    }

    keep_x <- rep(TRUE,nrow(x))
    keep_y <- rep(FALSE,nrow(y))

    x_sn <- x$startNode
    x_en <- x$endNode
    y_sn <- y$startNode
    y_en <- y$endNode

    for(ii in 1:nrow(y)){
        idx <- terra::is.related(x, y[ii,], "intersects")
        idx <- idx & keep_x ## drop any intersection with already dropped items
        if(any(idx)){
            keep_y[ii] <- TRUE

            if(sum(idx)==1){
                ## assume an outlet else there should be two intersections
                ## presumes no intersection at the outlet of x
                edx <- FALSE
                sdx <- TRUE
            }else{
                edx <- x_en[idx] %in% x_sn[idx] ## true is end of reach is start of another reach in subset
                sdx <- x_sn[idx] %in% x_en[idx] ## true is start is also and endNode for another reach i
                ## edx <- x$endNode[idx] %in% x$startNode[idx] ## true is end of reach is start of another reach in subset
                ## sdx <- x$startNode[idx] %in% x$endNode[idx] ## true is start is also and endNode for another reach in subet
            }

            ## redo end nodes
            tmp <- unique( x_en[idx][ edx ] ) ## end Nodes in the new object
            x_en[ x_en %in% tmp ] <- y_sn[ii]

            ## redo start nodes
            tmp <- unique( x_sn[idx][sdx] )
            x_sn[ x_sn %in% tmp ] <- y_en[ii]

            ## flag ones to remove
            keep_x[idx][ (edx & sdx) ] <- FALSE

        }

        if( verbose ){ setTxtProgressBar(pb, ii, title = NULL, label = NULL) }
    }

    x$endNode <- x_en
    x$startNode <- x_sn

    x <- rbind(x[keep_x,],y[keep_y,])

    terra::writeVector(x,"test_merge.gpkg",overwrite=TRUE)

    check_channel(x,outlets,check_connectivity=check_connectivity)

    return(x)
}



#' simplify and channel object
#'
#' @description Simplify a channel object by merging spatial objects
#'
#' @param chn the channel object
#' @param simplify_length minimum length of channel sections to remove
#' @param outlets specify outlet reaches to check the simplified channel on exit
#' @param strict_routing see details
#' @param verbose print a higher level of output
#'
#' @details strict_routing will maintain the flow pathways by only removing nodes with a single input and output, otherwise the routing may alter the effective flow lengths.
#' The is_wb feild in the data is used to flag objects that shouldn't be merged..
#' @export
simplify_channel <- function(chn, simplify_length=100, outlets=NULL, strict_routing = TRUE, verbose = FALSE){
   stop("Not updated")

    ## check chn is a channel
    if(is.character(chn)){ chn <- terra::vect(chn) }
    check_channel(chn,outlets)

    if( simplify_length <= 0 ){
        stop("simplify_length must be positive")
    }


    ## split geom and data - seems quicker...
    gm <- chn
    terra::values(gm) <- NULL
    gm <- split(gm,1:nrow(gm))
    chn <- as.data.frame(chn)

    loop_flag <- TRUE
    lp <- 0

    while( loop_flag ){
        lp <- lp + 1
        print(paste("starting loop",lp))


        to_keep <- rep(TRUE,nrow(chn))

        if("is_wb" %in% names(chn)){ not_wb <- !chn$is_wb }
        else{ not_wb <- rep(TRUE, nrow(chn)) }

        ## find channels to merge and order
        idx <- which( (chn$length < simplify_length) & not_wb)
        idx <- idx[order(chn$length[idx])]

        n <- length(idx)
        if( verbose ){
            pb <- txtProgressBar(min = 0, max = n, initial = 0, char = "=",
                                 width = NA, "progress...", "whoop", style = 3, file = "")
            on.exit( close(pb) )
            cnt <- 0
        }



        for(ii in idx){
            if( chn$length[ii] >= simplify_length ){
                ## catch incase a merge has already made it long
                cnt <- cnt+1
                setTxtProgressBar(pb, cnt, title = NULL, label = NULL)
                next
            }

            in_hn <-  which( chn$endNode == chn$startNode[ii] )
            out_hn <-  which( chn$startNode == chn$startNode[ii] )
            in_en <- which( chn$endNode == chn$endNode[ii] )
            out_en <- which( chn$startNode == chn$endNode[ii] )
            jj <- NA
            if( length(in_en)==1 &&
                length(out_en)==1 &&
                not_wb[out_en] ){
                ## endNode is a 1-to-1 join and can merge d/s
                jj <- out_en
                chn$startNode[jj] <- chn$startNode[ii]
            }
            if( length(in_hn)==1 &&
                not_wb[in_hn] &&
                length(out_hn)==1 &&
                is.na(jj) ){
                ## startNode is a 1-to-1 join and can merge u/s
                jj <- in_hn
                chn$endNode[jj] <- chn$endNode[ii]
            }
            if( length(out_hn)==1 &&
                is.na(jj) &&
                !strict_routing ){
                ## remove the channel entirely...
                jj <- NA
                if( length(out_en) > 0 &&
                    any(not_wb[out_en]) ){
                    ## merge with one downstream
                    jj <- out_en[not_wb[out_en]][1]
                }
                if( is.na(jj) &&
                    length(in_hn) > 0 &&
                    any(not_wb[in_hn]) ){
                    ## merge with one upstream
                    jj <- in_hn[not_wb[in_hn]][1]
                }
                if(!is.na(jj)){
                    chn$endNode[in_hn] <- chn$endNode[ii]
                }


                ## ## endNode is a single outlet so merge d/s but reroute other flows to
                ## ## new start of reach
                ## jj <- out_en
                ## chn$endNode[ in_en ] <- chn$startNode[ii]
                ## chn$startNode[jj] <- chn$startNode[ii]
            }

            if(!is.na(jj)){ ## can simplify

                ##print(paste(cnt,ii,jj))

                ## merge properties (default to those for jj)
                chn$length[jj] <- chn$length[jj] + chn$length[ii]
                chn$area[jj] <- chn$area[jj] + chn$area[ii]
                chn$width[jj] <- chn$area[jj] / chn$length[jj]
                chn$channelVol[jj] <- chn$channelVol[jj] + chn$channelVol[ii]
                chn$endNode[ii] <- chn$startNode[ii] <- NA # to stop matching on deleted segments

                ## merge geom
                gm[[jj]] <- terra::combineGeoms( gm[[jj]], gm[[ii]], minover=0 )
                to_keep[ii] <- FALSE
            }

            if( verbose ){
                cnt <- cnt+1
                setTxtProgressBar(pb, cnt, title = NULL, label = NULL)
            }
        }

        close(pb)

        if( all(to_keep) ){ loop_flag <- FALSE }
        else{
            gm <- gm[to_keep]
            chn <- chn[to_keep,]
        }

    }

    ## merge back into data.frame
    chn <- cbind(terra::vect(gm[to_keep]),chn[to_keep,])

    check_channel(chn,outlets)

    return(chn)
}

#' Function for checking channel networks are valid and suitable for dynatopGIS
#'
#' @description Performs chacks ont he values and properties of a channel then checks connnections to outlets if given
#' @param chn a SpatVect object representing the channel newtork as returned by convert_channel
#' @param outlets uid values of the outlets of the channel network
#'
#' @return TRUE if completes else failure message
#'
#' @details The processing follows the follwoing sequence:
#'   - Checks of the required properties of the channel netowrk for use in dynatopGIS
#'   - checks on the connectivity for loops
#'   - If outlets provided comparision of identified and known outlets
#'
#' and checks on the comparision fThe `property_names` input is used to rename the data associated with spatial objects. Remaining data is dropped
#'   - Varaibles are converted to the expected type then checked for missing values; in doing this
#'       - Missing `slope`, `width` and `depth`  values are populated by the default and a warning issued
#'   - If the spatial objects are lines
#'      - The spatial objects are buffered using the half `width`
#'
#' @examples
#' stop("example not updated")
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' vect_lines <- terra::vect(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(vect_lines,property_names)
#' @export
check_channel <- function(chn,outlets=NULL,chn_is_lines=FALSE,check_connectivity=TRUE){
    ## direct checks
    allowedGeom <- ifelse(chn_is_lines,"lines","polygons")
    stopifnot(
        "The channel network does not have a suitable geometry" = terra::geomtype(chn) == allowedGeom,
        ## names
        "uid property is missing" = "uid" %in% names(chn),
        "uids should be strings" = is.character(chn$uid),
        "uids cannot be missing" = !any(is.na(chn$uid)),
        "uids should be unique" = !any(duplicated(chn$uid)),
        ## startNode
        "startNode property is missing" = "startNode" %in% names(chn),
        "startNodes should be strings" = is.character(chn$startNode),
        "startNodes cannot be missing" = !any(is.na(chn$startNode)),
        ## endNode
        "endNode property is missing" = "endNode" %in% names(chn),
        "endNodes should be strings" = is.character(chn$endNode),
        "endNodes cannot be missing" = !any(is.na(chn$endNode)),
        ## depth
        "depth property is missing" = "depth" %in% names(chn),
        "depths should be numeric" = is.numeric(chn$depth),
        "depths cannot be missing" = !any(is.na(chn$depth)),
        "depths should be posistive" = all(chn$depth>0),
        ## slope
        "slope property is missing" = "slope" %in% names(chn),
        "slopes should be numeric" = is.numeric(chn$slope),
        "slopes cannot be missing" = !any(is.na(chn$slope)),
        "slopes should be positive" = all(chn$slope>0),
        ## width
        "width property is missing" = "width" %in% names(chn),
        "widths should be numeric" = is.numeric(chn$width),
        "widths cannot be missing" = !any(is.na(chn$width)),
        "widths should be positive" = all(chn$width>0)
    )
    ## check outlets
    if( !is.null(outlets) ){
        chn_outlets <- chn$uid[ !(chn$endNode %in% chn$startNode) ]
        if( !all(outlets %in% chn_outlets) ){
            stop(paste("Missing channel outlets:",paste( setdiff(outlets,chn_outlets), collapse=", ")))
        }
        if( !all(chn_outlets %in% outlets) ){
            stop(paste("Additional channel outlets:",paste( setdiff(chn_outlets,outlets), collapse=", ")))
        }
    }

    ## check connectivity
    if(check_connectivity){
        ## This is much quicker using vectors and not constantly accessing via the SpatVect object
        sN <- chn$startNode
        eN <- chn$endNode

        idx <- !(eN %in%sN) ## outlets are channel lengths whose outlet does not join another channel
        it <- 0
        cnt <- table(sN) ## we should never vist a node more times then it is a starting point

        while( any(idx) & it <= nrow(chn)){
            jdx <- sN[idx]
            tmp <- table(jdx)
            cnt[ names(tmp) ] <- cnt[ names(tmp) ] - tmp
            if( any(cnt<0) ){
                stop(paste("Failing loop involving nodes:", paste(names(cnt)[cnt<0],collapse=", ")))
            }
            jdx <- jdx[cnt[jdx]==0] ## only move up if it is the last visit to the startNode
            idx <- eN %in% jdx
            it <- it+1
        }
        ##browser()
        ##pjs <- 999
        stopifnot(
            "Error ingesting channel: problem with visiting all points" = all(cnt==0),
            "To many iterations" = it <= nrow(chn)
        )
    }

    invisible(TRUE)
}


#' Locate flow gauges on a channel network
#'
#' @description Locate flow gauges on a channel network
#'
#' @param chn the channel object
#' @param gauges either points of gauge locations of polygons of catchment boundaries
#' @param gauge_name feild name in gauges to use for name of gauge
#' @param max_dist maxium distance between gauge points and channel
#'
#' @details If gauges are points the nearest channel length within max_dist is chosen. If catchment polygons are given all channel lengths leaving the catchment are returned.
#' @export
locate_gauges <- function(chn,gauges,gauge_name="name",max_dist = 100){
    stop("Not updated")
    ## check chn is a channel
    if(is.character(chn)){ chn <- terra::vect(chn) }
    check_channel(chn)

    ## check y is a terra object of points of polygons
    if(is.character(gauges)){ gauges <- terra::vect(gauges) }
    stopifnot("The gauges do not have a point or polygon geometry (even when read in)" =
                  terra::geomtype(gauges) %in% c("points","polygons"),
              "The gauge_name field does not exist" = gauge_name %in% names(gauges)
              )

    if( terra::geomtype(gauges) == "points" ){
        g_buff <- terra::buffer(gauges, width=max_dist)
        idx <- terra::is.related(chn, g_buff, "intersects")
        chn <- chn[idx,]
        dst <- terra::distance(gauges,chn) ## distance between gauges(row) and channels(column)
        idx <- apply(dst,1,which.min) ## index of closest channel for each gauge
        out <- data.frame(name = terra::values(gauges)[[gauge_name]],
                          channel_name = chn$name[ apply(dst,1,which.min) ],
                          distance = apply(dst,1,min))
        out$channel_name[ out$distance > max_dist ] <- NA
    }else{
        out <- list()
        gnm <- terra::values(gauges)[[gauge_name]]
        for(ii in 1:nrow(gauges)){
            idx <- terra::is.related(chn, gauges[ii,], "intersects")
            if( any(idx) ){
                out[[ii]] <- data.frame(name = gnm[ii],
                                        channel_name = chn$name[idx][ !(chn$endNode[idx] %in% chn$startNode[idx]) ]
                                        )
            }else{
                out[[ii]] <- data.frame(name = gnm[ii],
                                        channel_name = NA_character_
                                        )
            }

        }
        out <- do.call(rbind,out)
    }
    return(out)
}

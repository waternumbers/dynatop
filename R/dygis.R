#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
#' @examples
#' ## The vignettes contains more examples of the method calls.
#'
#' ## create temporary directory for output
#' demo_dir <- tempfile("dygis")
#' dir.create(demo_dir)
#'
#' ## initialise processing
#' ctch <- dynatopGIS$new(file.path(demo_dir,"test.tif"))
#'
#' ## add a catchment outline based on the digital elevation model
#' dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
#' dem <- terra::rast(dem_file)
#' dem <- terra::extend(dem,1)
#' catchment_outline <- terra::ifel(is.finite(dem),1,NA)
#' ctch$add_catchment(catchment_outline)
#'
#' ## add digital elevation and channel data
#' ctch$add_dem(dem)
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' sp_lines <- terra::vect(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(sp_lines,property_names)
#' ctch$add_channel(chn)
#'
#' ## compute properties
#' ctch$sink_fill() ## fill sinks in the catchment and computes dem flow directions
#' \donttest{
#' ##ctch$compute_band()
#' ctch$compute_properties() # like topograpihc index and contour length
#' ##ctch$compute_flow_lengths()
#' }
#' ## classify and create a model
#' \donttest{
#' ctch$add_layer( ctch$classify("atb_20","atb",cuts=20) )# classify using the topographic index
#' ##ctch$get_method("atb_20") ## see the details of the classification
#' ctch$combine_classes("atb_20_band",c("atb_20","band")) ## combine classes
#' ctch$create_model(file.path(demo_dir,"new_model"),"atb_20") ## create a model
#' list.files(demo_dir,pattern="new_model*") ## look at the output files for the model
#' }
#' ## tidy up
#' unlink(demo_dir)
#' @export
dynatopGIS <- R6::R6Class(
    "dynatopGIS",
    public = list(
        #' @description Initialise a project, or reopen an existing project
        #'
        #' @param projectFile Either the tif file of spatial maps or the correspondingly names shapefile of channels
        #'
        #' @details This loads either the tif or shape file, looks for the other then does some basic consistancy checks. If neither exist nothing is done but the filename is recorded.
        #'
        #' @return A new `dynatopGIS` object
        initialize = function(projectFile){
            private$apply_initialize( projectFile )
            invisible(self)
        },
        #' @description Add a catchment outline to the `dynatopGIS` project
        #'
        #' @param catchment a \code{SpatRaster} object or the path to file containing one which contains a rasterised catchment map.
        #'
        #' @details If not a \code{SpatRaster} object the the catchment is read in using the terra package. Finite values in the raster indicate that the area is part of the catchment; with each subcatchment taking a unique finite value. Note that in the later processing it is assumed that outflow from the subcatchments can occur only through the channel network. The resolution and projection of the project is taken from the provided catchment
        #'
        #' @return \code{invisible(self)}
        add_catchment = function(catchment){
            if(!("SpatRaster" %in% class(catchment))){ catchment <- terra::rast(as.character(catchment)) }
            if(!("SpatRaster" %in% class(catchment))){ stop("catchment is not a SpatRaster") }
            private$apply_add_catchment(catchment)
            invisible(self)
        },
        #' @description Import a dem to the `dynatopGIS` object
        #'
        #' @param dem a \code{raster} layer object or the path to file containing one which is the DEM
        #' @param fill_na  should NA values in dem be filled. See details
        #' @param verbose Should additional progress information be printed
        #'
        #' @details If not a \code{raster} the DEM is read in using the terra package. If \code{fill_na} is \code{TRUE} all NA values other then those that link to the edge of the dem are filled so they can be identified as sinks.
        #'
        #' @return suitable for chaining
        add_dem = function(dem,fill_na=-9999){
            if(!("SpatRaster" %in% class(dem))){ dem <- terra::rast(as.character(dem)) }
            if(!("SpatRaster" %in% class(dem))){ stop("dem is not a SpatRaster") }
            private$apply_add_dem(dem,fill_na)
            invisible(self)
        },
        #' @description Import channel data to the `dynatopGIS` object
        #'
        #' @param channel a SpatVect object or file path that can be loaded as one containing the channel information
        #' @param verbose Should additional progress information be printed
        #' @details Takes the representation of the channel network as a SpatVect with properties name, length, area, startNode, endNode and overlaying it on the DEM. In doing this a variable called id is created (or overwritten) other variables in the data frame are passed through unaltered.
        #'
        #' @return suitable for chaining
        add_channel = function(channel,verbose=FALSE){
            if(!is(channel,"SpatVector")){ channel <- terra::vect( as.character(channel) ) }
            if(!is(channel,"SpatVector")){ stop("channel is not a SpatVector object") }

            private$apply_add_channel(channel,as.logical(verbose))
            invisible(self)
        },
        #' @description Add a layer of geographical information
        #'
        #' @param layer the raster layer to add (see details)
        #' @param layer_name name to give to the layer
        #'
        #' @details The layer should either be a raster layer or a file that can be read by the \code{raster} package. The projection, resolution and extent are checked against the existing project data. Only layer names not already in use (or reserved) are allowed. If successful the layer is added to the project tif file.
        #' @return suitable for chaining
        add_layer = function(layer,layer_name=names(layer)){
            layer_name <- as.character(layer_name)
            if(!("SpatRaster" %in% class(layer))){ layer <- terra::rast(as.character(layer)) }
            if(!("SpatRaster" %in% class(layer))){ stop("layer is not a SpatRaster") }
            if( length(layer_name) != terra::nlyr(layer) ){ stop("Length of names does not match number of layers") }
            private$apply_add_layer(layer,layer_name)
            invisible(self)
        },
        #' @description Get a layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a `raster` layer of the requested information if layer_name is given else a vector of layer names
        get_layer = function(layer_name=NULL){
            ## create a vector of available layers
            tmp <- names(private$brk)
            if( "channel" %in% tmp ){ tmp <- c(tmp,"channel_vect") }

            if(is.null(layer_name)){ return(tmp) }

            ## check layer name exists
            layer_name <- match.arg(layer_name,tmp,several.ok=TRUE)

            ## make raster and return
            if( "channel_vect" %in% layer_name){
                return( private$chn )
            }else{
                return( private$brk[[layer_name]] )
            }
        },
        #' @description Plot a layer
        #' @param layer_name the name of layer to plot
        #' @param add_channel should the channel be added to the plot
        #' @return a plot
        plot_layer = function(layer_name,add_channel=TRUE){
            layer_name <- layer_name[1]
            lyr <- self$get_layer(layer_name)
            terra::plot( lyr, main = layer_name)
            if( add_channel & length(private$chn) > 0){
                terra::plot(private$chn, add=TRUE )
            }
        },
        #' @description The sink filling algorithm of Planchona and Darboux (2001)
        #'
        #' @param min_grad Minimum gradient between cell centres
        #' @param max_it maximum number of replacement cycles
        #' @param verbose print out additional diagnostic information
        #' @param hot_start start from filled_dem if it exists
        #' @param flow_type The type of flow routing to apply see details
        #' @details The algorithm implemented is based on that described in Planchona and Darboux, "A fast, simple and versatile algorithm to fill the depressions in digital elevation models" Catena 46 (2001). A pdf can be found at (<https://horizon.documentation.ird.fr/exl-doc/pleins_textes/pleins_textes_7/sous_copyright/010031925.pdf>). The adaptations made are to ensure that all cells drain only within the subcatchments if provided.
        #'
        #' The flow_type can be either
        #' - "quinn" where flow is split across all downslope directions or
        #' - "d8" where all flow follows the steepest between cell gradient
        #'
        sink_fill = function(min_grad = 1e-4,max_it=1e6,verbose=FALSE, hot_start=FALSE, flow_type=c("quinn","d8")){
            flow_type <- match.arg(flow_type)
            private$apply_sink_fill(min_grad,max_it,verbose,hot_start,flow_type)
            invisible(self)
        },
        ## #' @description Computes the computational band of each cell
        ## #'
        ## #' @param type type of banding
        ## #' @param verbose print out additional diagnostic information
        ## #'
        ## #' @details Banding is used within the model to define the HRUs and control the order of the flow between them; HRUs can only pass flow to HRUs in a lower numbered band. Currently only a strict ordering of river channels and cells in the DEM is implemented. To compute this the algorithm passes first up the channel network (with outlets being in band 1) then through the cells of the DEM in increasing height.
        ## compute_band = function(type=c("strict"), verbose=FALSE){
        ##     type = match.arg(type)
        ##     private$apply_band(type,verbose)
        ##     invisible(self)
        ## },
        #' @description Computes statistics e.g. gradient, log(upslope area / gradient) for raster cells
        #'
        #' @param min_grad gradient that can be assigned to a pixel if it can't be computed
        #' @param verbose print out additional diagnostic information
        #'
        #' @details The algorithm passed through the cells in decreasing height. Min grad is applied to all cells. It is also used for missing gradients in pixels which are partially channel but have no upslope neighbours.
        compute_properties = function(min_grad = 1e-4,verbose=FALSE){
            private$apply_upward_pass(verbose)
            private$apply_downward_pass(min_grad,verbose)
            invisible(self)
        },
        ## #' @description Computes flow length for each pixel to the channel
        ## #'
        ## #' @param flow_routing TODO
        ## #' @param verbose print out additional diagnostic information
        ## #'
        ## #' @details The algorithm passes through the cells in the DEM in increasing height. Three measures of flow length to the channel are computed. The shortest length (minimum length to channel through any flow path), the dominant length (the length taking the flow direction with the highest fraction for each pixel on the path) and expected flow length (flow length based on sum of downslope flow lengths based on fraction of flow to each cell). By definition cells in the channel that have no land area have a length of NA.
        ## compute_flow_lengths = function(flow_routing=c("expected","dominant","shortest"), verbose=FALSE){
        ##     flow_routing = match.arg(flow_routing)
        ##     private$apply_flow_lengths(flow_routing,verbose)
        ##     invisible(self)
        ## },
        #' @description Create a catchment classification based cutting an existing layer into classes
        #' @param layer_name name of the new layer to create
        #' @param base_layer name of the layer to be cut into classes
        #' @param cuts values on which to cut into classes. These should be numeric and define either the number of bands (single value) or breaks between band (multiple values).
        #'
        #' @details This applies the given cuts to the supplied landscape layer to produce areal groupings of the catchment. Cuts are implement using \code{terra::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        classify = function(layer_name,base_layer,cuts){
            private$apply_classify(as.character(base_layer[1]), cuts, as.character(layer_name[1]))
            ##invisible(self)
        },
        #' @description Combine any number of classifications based on unique combinations and burns
        #' @param layer_name name of the new layer to create
        #' @param pairs a vector of layer names to combine into new classes through unique combinations. Names should correspond to raster layers in the project directory.
        #' @param burns a vector of layer names which are to be burnt on
        #'
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given, only positive burn values are applied.
        combine_classes = function(layer_name,pairs,burns=NULL){
            private$apply_combine_classes(as.character(layer_name[1]),
                                          as.character(pairs),
                                          as.character(burns))
        },
        #' @description Compute a Dynamic TOPMODEL
        #'
        #' @param layer_name name for the new model and layers
        #' @param class_layer the layer defining the topographic classes
        #' @param sf_opt Surface solution to use
        #' @param sz_opt transmissivity profile to use
        #' @param rain_layer the layer defining the rainfall inputs
        #' @param rain_label Prepended to rain_layer values to give rainfall series name
        #' @param pet_layer the layer defining the pet inputs
        #' @param pet_label Prepended to pet_layer values to give pet series name
        #' @param min_grad minimum gradient between cell centres (or channel reaches)
        #' @param verbose print more details of progress
        #'
        #' @details The \code{class_layer} is used to define the HRUs. Flow between HRUs is based on the ordering of the catchment (see the \code{compute_band} method). Flow from a HRU can only go to a HRU with a lower band.
        #' Setting the sf_opt and sz_opt options ensures the model is set up with the correct parameters present.
        #' The \code{rain_layer} (\code{pet_layer}) can contain the numeric id values of different rainfall (pet) series. If the value of \code{rain_layer} (\code{pet_layer}) is not \code{NULL} the weights used to compute an averaged input value for each HRU are computed, otherwise an input table for the models generated with the value "missing" used in place of the series name.
        create_model = function(layer_name,class_layer,
                                sf_opt = c("cnst","kin"),
                                sz_opt = c("exp","bexp","cnst","dexp"),
                                rain_layer=NULL, rain_label=character(0),
                                pet_layer=NULL, pet_label=character(0),
                                min_grad = 1e-4,
                                verbose=FALSE){

            ## check valid transmissivity and channel_solver
            sf_opt<- match.arg(sf_opt)
            sz_opt <- match.arg(sz_opt)

            private$apply_create_model(class_layer,
                                       rain_layer, rain_label,
                                       pet_layer, pet_label,
                                       layer_name,verbose,
                                       sf_opt, sz_opt,min_grad)

            invisible(self)
        },
        #' @description get the version number
        #' @return a numeric version number
        #' @details the version number indicates the version of the algorithms within the object
        get_version = function(){
            private$version
        }
        ## #' @description get the cuts and burns used to classify
        ## #' @param layer_name the name of layer whose classification method is returned
        ## #' @return a list with two elements, cuts and burns
        ## get_method = function(layer_name){
        ##     ## check layer name exists
        ##     layer_name <- match.arg(layer_name,names(private$brk))

        ##     jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[layer_name]])),".json")
        ##     if( !file.exists(jsonFile) ){
        ##         stop("No json file giving basis of the classifications")
        ##     }

        ##     return( jsonlite::fromJSON( jsonFile ) )
        ## }
    ),
    private = list(
        version = "0.4.0",
        projectFile = character(0),
        brk = NULL,
        chn = NULL,
        reserved_layers = c("catchment","dem","channel","channel_fraction","filled_dem",
                            "gradient","upslope_area","atb",
                            "band"),
        ## save and reload changes
        save_project = function(chn=FALSE){
            tmpFile <- paste0(private$projectFile,".tmp")
            terra::writeRaster(private$brk,tmpFile, filetype="Gtiff",datatype="FLT8S",overwrite=TRUE)
            file.copy(tmpFile,private$projectFile,overwrite=TRUE)
            private$brk <- terra::rast( private$projectFile )

            if(chn){
                channelFile <- paste0(tools::file_path_sans_ext(private$projectFile),".geojson")
                tmpFile <- paste0(channelFile,".tmp")
                terra::writeVector(private$chn,tmpFile,filetype="GeoJSON",overwrite=TRUE)
                file.copy(tmpFile,channelFile,overwrite=TRUE)
                private$chn <- terra::vect(channelFile)
            }
        },
        ## check and read the project files
        apply_initialize = function(projectFile){
            ## check extension
            fileType <- tools::file_ext(projectFile)
            stopifnot( "Incorrect file extension" = fileType == "tif" )

            ## see if file exists
            if( !file.exists(projectFile) ){
                message( "Creating a new project" )
                private$projectFile <- projectFile
                return()
            }

            ## read in raster data
            message( paste("Reading existing project at", projectFile) )
            brk <- terra::rast( projectFile )

            ## work out the channel file location
            channelFile <- paste0(tools::file_path_sans_ext(projectFile),".geojson")
            chn <- NULL

            ## initial error checks
            stopifnot("Raster data not valid: Processing currently only works on projected data with a square grid" =
                          all.equal(diff(terra::res(brk)),0) & !terra::is.lonlat(brk),
                      "No catchment layer in Raster data" = "catchment" %in% names(brk),
                      "Channel data file missmatch - missing file" = !("channel" %in% names(brk)) | ( ("channel" %in% names(brk)) & file.exists(channelFile) )
                      )

            ## read in shape file and check
            if( "channel" %in% names(brk) ){
                chn <- terra::vect( channelFile )

                if(  terra::crs(brk, proj=TRUE) != terra::crs(chn, proj=TRUE) ){
                    stop("Shape and Raster projections do not match")
                }
            }



            private$projectFile <- projectFile
            private$brk <- brk
            private$chn <- chn
        },
        ## add a catchment map
        apply_add_catchment = function(catchment){

            ## check is padded
            isPadded <- !(any( is.finite(catchment[1,][[1]]) ) |
                          any( is.finite(catchment[nrow(catchment),][[1]]) ) |
                          any( is.finite(catchment[,1][[1]]) ) |
                          any( is.finite(catchment[,ncol(catchment)][[1]]) ))

            ## check is projected
            isProjected <- all.equal(diff(terra::res(catchment)),0) & !terra::is.lonlat(catchment)

            stopifnot(
                "The catchment map already exists, start a new project" = !("catchment" %in% names(private$brk)),
                "The catchment must be one layer" = terra::nlyr(catchment)==1,
                "Processing currently only works on projected maps with a square grid" = isProjected,
                "Catchment must be padded with NA rows and columns" = isPadded
            )

            ## save
            names(catchment) <- "catchment"
            private$brk <- catchment
            private$save_project()
        },
        ## adding dem
        apply_add_dem = function(dem,fill_na){

            stopifnot(
                "The catchment map does not exists, try running add_catchment" =
                    "catchment" %in% names(private$brk),
                "The DEM already exists" = !("dem" %in% names(private$brk)),
                "New layer does not match resolution, extent or projection of project" =
                    terra::compareGeom(dem,private$brk,stopOnError=FALSE),
                "The dem must be one layer" = terra::nlyr(dem)==1
            )

            if( terra::global(is.na(dem) & !is.na(private$brk[["catchment"]]),max)>0 ){
                stop("dem has missing values in the catchment - try running fill_na first")
            }

            ## tidy up dem
            dem <- terra::mask(dem,private$brk[[ "catchment" ]])
            names(dem) <- "dem"

            ## save
            private$brk <- c(private$brk,dem)
            private$save_project()
        },
        ## add the channel
        apply_add_channel = function(chn,verbose){

            rq <- c("catchment")
            chn_variables <- c(
                "name" = "character",
                "length" = "numeric",
                "area" = "numeric",
                "startNode" = "character",
                "endNode" = "character",
                "slope" = "numeric"
            )

            stopifnot(
                "Not all required layers are available" = all(rq %in% names(private$brk)),
                "Channel already added" = is.null(private$chn),
                "A required property name is not specified" = all(names(chn_variables) %in% names(chn)),
                "Projection of channel object does not match that of project" =
                     terra::crs(private$brk, proj=TRUE) == terra::crs(chn, proj=TRUE)
            )

            ## check if there is an id feild which will be overwritten
            if( ("id" %in% names(chn)) ){
                warning("The name id is reserved and will be overwritten in the channel import")
            }

            ## ensure required properties are of correct type
            for(ii in names(chn_variables)){
                chn[[ii]][[ii]] <- as(chn[[ii]][[ii]],chn_variables[ii])
            }

            stopifnot(
                "Some non-finite values of length found!" = all(is.finite(chn$length)),
                "Some non-finite values of area found!" = all(is.finite(chn$area)),
                "Some non-finite values of slope found!" = all(is.finite(chn$slope)),
                "Some non-finite values of area found!" = all( is.finite(chn$area) ),
                "Some non-positive values of length found!" = all( chn$length > 0 ),
                "Some non-positive values of area found!" = all( chn$area > 0 ),
                "Some non-positive values of slope found!" = all( chn$slope > 0 ),
                "Some non-positive values of area found!" = all( chn$area > 0)
            )

            ## arrange id and band in order of flow direction - so lowest values at outlets of the network
            if( verbose ){ print("Computing channel id's and bands") }
            ## This is much quicker using vectors and not constantly accessing via the vect object
            id <- rep(as.integer(0),nrow(chn))
            bnd <- rep(as.integer(0),nrow(chn))
            sN <- chn$startNode
            eN <- chn$endNode

            idx <- !(eN %in%sN) ## outlets are channel lengths whose outlet does not join another channel
            it <- 1
            cnt <- table(sN) ## we should never vist a node more times then it is a starting point
            while(sum(idx)>0){
                id[idx] <- max(id) + 1:sum(idx)
                bnd[idx] <- it

                jdx <- sN[idx]
                for(ii in jdx){ cnt[ii] <- cnt[ii] - 1 } ## since sN might appear more then once..
                if( any(cnt<0) ){
                    stop(paste("Failing loop involving nodes:", paste(names(cnt)[cnt<0],collapse=", ")))
                }
                jdx <- jdx[cnt[jdx]==0] ## only move up if it is the last visit to the startNode
                idx <- eN %in% jdx
                it <- it+1
            }
            chn$id <- id
            chn$band <- bnd

            stopifnot(
                "Error ingesting channel: check connectivity" = all(chn$id > 0),
                "Error ingesting channel: problem with bands" = all(chn$band > 0),
                "Error ingesting channel: problem with visiting all points" = all(cnt==0)
            )

            chn <- chn[ order(chn$id),]

            ## TODO - possibly sort on length to try to identify bigger channels??
            ## create a raster of channel id numbers
            chn_rst <- terra::rasterize(chn,private$brk[["catchment"]],field = "id",touches=TRUE)
            chn_rst <- terra::mask(chn_rst,private$brk[["catchment"]])
            names(chn_rst) <- "channel"

            ## create a raster of channel coverage fractions
            chn_frac <- terra::rasterize(chn,private$brk[["catchment"]],background=0,cover=TRUE) ## fraction of cell covered by channel
            chn_frac <- terra::mask(chn_frac,private$brk[["catchment"]])
            names(chn_frac) <- "channel_fraction"
            terra::values(chn_frac) <- round(terra::values(chn_frac),2)## else get horrible rounding errors close to 1
            ## add a fraction to those cells with an ID but no fractions
            chn_frac[chn_frac==0 & !is.na(chn_rst)] <- 0.005

            ## rescale channel fractions which are <1 to match channel area
            cell_area <- prod(terra::res(chn_frac))
            chn_area <- sum(chn$area)
            tmp <- terra::values(chn_frac)
            idx <- is.finite(tmp) & tmp<1
            da <-  chn_area -  sum(tmp,na.rm=TRUE)*cell_area
            it <- 0
            while( abs(da) > 1e-6 & it<100){
                sc <- 1 + da / (sum( tmp[idx] )*cell_area)
                tmp[idx] <- pmin(1, tmp[idx]*sc)
                da <-  chn_area -  sum(tmp,na.rm=TRUE)*cell_area
                it <- it + 1
            }
            terra::values(chn_frac) <- tmp

            ## ## work out which cells each channel touches
            ## tmp <- terra::cells(private$brk[["catchment"]],chn)
            ## tmp <- split(tmp[,2],tmp[,1])
            ## tmp <- sapply(tmp,paste,collapse=",")
            ## chn$input_cells <- tmp

            ## save output
            private$brk <- c(private$brk,chn_rst,chn_frac)
            private$chn <- chn
            private$save_project(chn=TRUE)
        },
        ## Add a layer
        apply_add_layer=function(layer,layer_name){

            rq <- c("channel","catchment")
            stopifnot(
                "Missing required layers" = all(rq %in% names(private$brk)),
                "Name is reserved" = !any(layer_name %in% private$reserved_layers),
                "Name is already used" = !any(layer_name %in% names(private$brk)),
                "New layer does not match resolution, extent or projection of project" =
                    terra::compareGeom(layer,private$brk,stopOnError=FALSE)
            )

            ## tidy up the layers to be added
            layer <- terra::mask(layer,private$brk[[ "catchment" ]])
            names(layer) <- layer_name

            ## check for no NA values
            ## TODO check how global works on a stack, do we just need to call once...
            for(ii in names(layer)){
                if( terra::global( is.na(layer[[ii]]) &
                                   !is.na(private$brk[["catchment"]]) &
                                   is.na(private$brk[["channel"]]) ,max)>0 ){
                    stop(paste("Layer",ii,"has missing values in the catchment - try running fill_na first"))
                }
            }

            ## save output
            private$brk <- c(private$brk,layer)
            private$save_project()
        },
        ## Sink fill
        apply_sink_fill = function(min_grad,max_it,verbose,hot_start,flow_type){
            ## recall the catchments is padded with NA values

            d <- ifelse(hot_start,"filled_dem","dem")

            rq <- c(d,"channel","catchment")

            if(!all(rq %in% names(private$brk))){
                stop("Not all required layers are available")
            }

            d <- terra::as.matrix( private$brk[[d]] ,wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE )
            ctch <- terra::as.matrix( private$brk[["catchment"]] , wide=TRUE )

            ## values that should be valid
            to_be_valid <- !is.na(ctch) #!is.na(d) & is.na(ch)  # all values not NA should have a valid height
            is_valid <- !is.na(ch) # & !is.na(d) # TRUE if a channel cell for initialisation
            changed <- is_valid # cells changed at last iteration
            fd <- to_be_valid*Inf; fd[is_valid] <- d[is_valid]
            to_eval <- is_valid; to_eval[] <- FALSE

            ## dimensions
            nf <- sum(!is.na(ctch))

            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- matrix(sqrt(sum(rs^2)),3,3)
            dxy[1,2] <- dxy[3,2] <- rs[2]
            dxy[2,1] <- dxy[2,3] <- rs[1]
            dxy[2,2] <- 0
            dxy <- min_grad*dxy

            it <- 1
            ## start of iteration loop
            while(any(changed[]) & it<=max_it){

                ## work out all the the cells to evaluate
                ## should be next to those that are changed
                to_eval[] <- FALSE
                idx <- which(changed,arr.ind=TRUE) # index of changed cells
                sctch <- ctch[idx] # sub catchment number
                ##if( any(is.na(sctch)) ){ browser() }
                jdx <- idx
                for(ii in c(-1,0,1)){
                    for(jj in c(-1,0,1)){
                        if(ii==0 & jj==0){next}
                        ## adjust jdx
                        ##browser()
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        ## trim so only evaluate
                        ## (ii) cells within the same subcatchment
                        cjdx <- ctch[jdx]
                        kdx <- !is.na(cjdx) & (sctch == cjdx)
                        jjdx <- jdx[kdx,, drop=F]
                        to_eval[jjdx] <- !is_valid[jjdx] #TRUE
                    }
                }
                to_eval <- to_eval & to_be_valid #& !is_valid
                if(verbose){
                    cat("Iteration",it,"\n")
                    cat("\t","Cells to evaluate:",sum(to_eval),"\n")
                    cat("\t","Percentage Complete:",
                        round(100*sum(is_valid)/nf,1),"\n") #to_be_valid),1),"\n")
                }
                ## alter min value for the evaluation cells
                idx <- which(to_eval,arr.ind=TRUE) # index of changed cells
                sctch <- ctch[idx] # sub catchment number
                ##if( any(is.na(sctch)) ){ browser() }
                jdx <- idx
                mind <- rep(Inf,nrow(idx))
                for(ii in c(-1,0,1)){
                    for(jj in c(-1,0,1)){
                        if(ii==0 & jj==0){next}
                        ## adjust jdx
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        ## trim so only evaluate
                        ## (ii) cells within the same subcatchment
                        cjdx <- ctch[jdx]
                        kdx <- !is.na(cjdx) & (sctch == cjdx)
                        if( !any(kdx) ){ next }
                        jjdx <- jdx[kdx,,drop=F]
                        ##if( any(is.na(kdx)) | any(is.na(jjdx)) ){ browser() }
                        ##if( sum(kdx)==0 ){ browser() }##| sum(kdx) != nrow(jjdx) ){ browser() }
                        tmp <- pmin(mind[kdx],fd[jjdx] + dxy[ii+2,jj+2],na.rm=TRUE)
                        ##if( length(tmp) != length( mind[kdx] )){ browser() }
                        mind[kdx] <- pmin(mind[kdx],fd[jjdx] + dxy[ii+2,jj+2],na.rm=TRUE)
                    }
                }
                changed[] <- FALSE
                is_valid[idx] <- d[idx]>mind ## cells where the dem value is valid
                mind <- pmax(d[idx],mind) ## mind is now the replacemnt value
                changed[idx] <- mind < fd[idx]
                fd[idx] <- mind

                ## end of loop
                it <- it+1
            }

            if(hot_start){
                terra::values(private$brk[["filled_dem"]]) <- fd
            }else{
                rfd <- terra::rast( private$brk[["dem"]], names="filled_dem", vals=fd )
                private$brk <- c( private$brk, rfd )
            }

            private$save_project()

            if(it>max_it){ stop("Maximum number of iterations reached, sink filling not complete") }

        },
        ## function to do property calculations on an upwards pass (low to high DEM values)
        ## if we go up in height order then we are working from near the channel to the heighest point
        ## could add back in flow distances here
        apply_upward_pass = function(verbose){

            rq <- c("filled_dem","channel","channel_fraction")
            stopifnot(
                "Not all required input layers have been generated \n Try running sink_fill first" =
                    all( rq %in% names( private$brk) )
            )

            ## rasterize channel band to start
            rbnd <- terra::rasterize(private$chn, private$brk[["catchment"]],field = "band",touches=TRUE)
            rbnd <- terra::mask(rbnd,private$brk$catchment) ## ensure channel bands are within the catchment - else later code fails
            names(rbnd) <- "band"

            ## load raster layer
            d <- terra::as.matrix( private$brk[["filled_dem"]], wide=TRUE )
            chn <- terra::as.matrix( private$brk[["channel"]], wide=TRUE )
            chn_frc <- terra::as.matrix( private$brk[["channel_fraction"]], wide=TRUE )
            bnd <- terra::as.matrix( rbnd,  wide=TRUE )

            if( verbose ){ print("Computing upward pass") }

            idx <- order(d,na.last=NA) ## search order
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1) ## neighbours

            ## set up printing variables
            if(verbose){
                print_step <- c(1,rep(round(length(idx)/20,2)),length(idx)) # current, next print, step, total
            }

            ## main loop
            dz <- rep(NA,8)
            for(ii in idx){

                if( is.finite(chn[ii]) ){ ## then cell is a channel
                    if(chn_frc[ii]<1){
                        ## since it is mixed cell - partly landuse , partly channel
                        bnd[ii] <- bnd[ii] + 1
                    }
                }else{
                    ## it is not a channel
                    jdx <- ii+delta
                    dz[] <- d[ii] - d[jdx]
                    is_lower <- is.finite(dz) & dz>0
                    bnd[ii] <- max( bnd[jdx[is_lower]] ) + 1
                }

                if(verbose){
                    print_step[1] <- print_step[1] + 1
                    if( print_step[1] > print_step[2] ){
                        cat(round(100*print_step[1] / print_step[4],1),
                            "% complete","\n")
                        print_step[2] <- print_step[2] + print_step[3]
                    }
                }

            }

            ## save
            terra::values(rbnd) <- bnd
            private$brk <- c(private$brk,rbnd)
            private$save_project()


        },
        apply_downward_pass = function(min_grad,verbose){

            if( verbose ){ print("Loading data for downward pass") }
            rq <- c("filled_dem","channel")
            stopifnot(
                "Not all required input layers have been generated \n Try running sink_fill first" =
                    all( rq %in% names( private$brk) )
            )

            ## load raster layer
            d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE)
            ch_frc <- terra::as.matrix( private$brk[["channel_fraction"]] , wide=TRUE)

            if( verbose ){ print("Setting up computation") }

            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs) ## assumes square
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

            ## initialise output
            gr <- upa <- atb <- d*NA
            upa <- prod(rs)*(1-ch_frc) ## initialise upslope area from resolution

            idx <- order(d,na.last=NA,decreasing=TRUE) ## search order

            ## set up printing variables
            if(verbose){
                print_step <- c(1,rep(round(length(idx)/20,2)),length(idx)) # current, next print, step, total
            }

            if( verbose ){ print("Computing hillslope") }

            uA <- private$chn$area ## upsteam areas for channels

            ## loop downslope
            w <- rep(0,8)
            for(ii in idx){
                if( is.finite(ch[ii]) ){
                    ## pass on upslope area to channel
                    uA[ ch[ii] ] <- uA[ ch[ii] ] + upa[ii]
                    if( ch_frc[ii] < 1 ){
                        ## mixed cells
                        ## work out gradient from cells flowing in
                        jdx <- ii+delta
                        grd <- (d[ii]-d[jdx])/dxy
                        gcl <- grd*dcl
                        is_higher <- is.finite(gcl) & gcl<0 #& is.finite(cjdx) & cjdx==ctch[ii]
                        if( any(is_higher) ){
                            sum_gcl <- sum( gcl[is_higher] )
                            sum_dcl <- sum( dcl[is_higher] )
                            gr[ii] <- max(-sum_gcl / sum_dcl,min_grad)
                        }else{ ## nothing flows into the cell
                            gr[ii] <- min_grad
                        }
                        atb[ii] <- log(upa[ii]/gr[ii])
                    }else{ ## pure water cell
                        upa[ii] <- NA
                        gr[ii] <- NA
                        atb[ii] <- NA
                    }

                }else{

                    ## it is not a channel
                    w[] <- 0
                    jdx <- ii+delta
                    grd <- (d[ii]-d[jdx])/dxy
                    gcl <- grd*dcl
                    is_lower <- is.finite(gcl) & gcl>0 #& is.finite(cjdx) & cjdx==ctch[ii]
                    sum_gcl <- sum( gcl[is_lower] )
                    sum_dcl <- sum( dcl[is_lower] )
                    w[is_lower] <- gcl[is_lower] / sum_gcl

                    if( !any(w>0) ){ stop(paste("Cell",ii,"is a hillslope cell with no outflows")) }

                    ## gradient
                    gr[ii] <- max(sum_gcl / sum_dcl,min_grad)
                    ## topographic index
                    atb[ii] <- log(upa[ii]/gr[ii]) #log( upa[ii] / sum(gcl) )
                    ## propogate area downslope
                    upa[ jdx ]  <- upa[ jdx ] + w*upa[ii]
                }

                ## verbose output here
                if(verbose){
                    print_step[1] <- print_step[1] + 1
                    if( print_step[1] > print_step[2] ){
                        cat(round(100*print_step[1] / print_step[4],1),
                            "% complete","\n")
                        print_step[2] <- print_step[2] + print_step[3]
                    }
                }
            }

            if( verbose ){ print("Computing channel") }
            sN <- private$chn$startNode
            eN <- private$chn$endNode


            ## merge upslope areas into the channel object
            ##ch_upa <- tapply(upa,ch,sum)
            ##ch_upa <- ch_upa[setdiff(names(ch_upa),"NaN")]
            ##idx <- match(names(ch_upa),paste(private$chn$id)) #,names(ch_upa))
            ##uA[idx] <- as.numeric(ch_upa)

            ## remove channel area bit from hillslope upslope area
            ##upa[ch_frc==1] <- NA

            ## compute catchment area to each reach
            for(ii in length(sN):1){
                idx <- sN == eN[ii]
                uA[idx] <- uA[idx] + uA[ii] / sum(idx) ## TO CHECK not sure how channel routing fractions originally done
            }
            stopifnot(
                "All channel upstream areas should be finite" = all(is.finite(uA)),
                "All channel upstream areas should be non-negative" = all(uA>0)
            )

            private$chn$upstream_area <- uA

            ## save output
            private$brk <- c(private$brk,
                             terra::rast( private$brk[["dem"]], names="gradient", vals=gr ),
                             terra::rast( private$brk[["dem"]], names="upslope_area", vals=upa ),
                             terra::rast( private$brk[["dem"]], names="atb", vals=atb )
                             )

            private$save_project(chn=TRUE)
        },
        ## split_to_class
        apply_classify = function(base_layer,cuts,layer_name){

            rq <- c("channel",base_layer)
            stopifnot(
                "Missing channel or base layer to classify" = all( rq %in% names(private$brk) ),
                "layer_name has zero length" = nchar(layer_name)>0,
                "layer_name is already used" = !(layer_name %in% names(private$brk)),
                "layer_name is reserved" = !(layer_name %in% names(private$reserved_layers))
            )

            ## load base layer and mask out channel
            ##x <-  terra::mask( private$brk[[base_layer]], private$brk[["channel"]], inverse=TRUE)
            x <- private$brk[[base_layer]]

            ## work out breaks
            brk <- as.numeric(cuts)
            rng <- as.numeric( terra::global(x, fun="range",na.rm=TRUE) )
            if( length(brk)==1 ){
                ## this defines brks in the same way as cut would otherwise
                brk <- seq(rng[1],rng[2],length=brk+1)
            }else{
                brk <- sort(brk)
                if( brk[1] > rng[1]){ brk <- c(rng[1],brk) }
                if( tail(brk,1) < rng[2]){ brk <- c(brk,rng[2]) }
            }
            if( any(is.na(brk)) ){ stop("NA value in brk") }
            M <- cbind( head(brk,-1), tail(brk,-1), 1:(length(brk)-1) )

            ## cut the raster and save

            return( terra::classify(x,rcl=M,include.lowest=TRUE,names=layer_name) )
        },
        ## split_to_class
        apply_combine_classes = function(layer_name,pairs,burns){

            stopifnot(
                "Missing layers in pairs list" = all(pairs %in% names(private$brk)),
                "Missing layers in burns list" = all(burns %in% names(private$brk)),
                "There must be at lest one entry in pairs" = length(pairs)>0
            )

            x <- terra::as.matrix( private$brk[[ pairs ]] )
            idx <- rowSums(is.finite(x)) == ncol(x) ## thses are the valid cells
            xstr <- apply(x,1,function(r){paste(r,collapse="_")})

            ## add burns sequentally
            if(length(burns)>0){
                y <- terra::as.matrix( private$brk[[ burns ]] )
                y[y<=0] <- NA
                for(ii in 1:ncol(y)){
                    jdx <- is.finite(y[,ii])
                    xstr[jdx] <- paste("burn",ii,y[jdx,ii],sep="_")
                }
            }

            ## make numeric class
            uxstr <- unique(xstr[idx])
            ux <- setNames(1:length(uxstr),uxstr)
            z <- rep(NA,nrow(x))
            z[idx] <- ux[xstr[idx]]

            return( terra::rast( private$brk[["dem"]], names=layer_name, vals=z ) )
        },

        ## create a model
        apply_create_model = function(class_lyr,
                                      rain_lyr,rainfall_label,
                                      pet_lyr,pet_label,
                                      layer_name,verbose,
                                      sf_opt,
                                      sz_opt,min_grad){


            ## check layers
            rq <- c("filled_dem","channel","channel_fraction",
                    "band","gradient",class_lyr,
                    rain_lyr,pet_lyr)

            stopifnot(
                "Missing layers" = all(rq %in% names(private$brk))
            )
            ## work out some properties of the brick and channel
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs) ## assumes square
            nr <- terra::ncol(private$brk) ##ra::nrow(private$brk)
            delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            cell_area <- prod(rs)
            n_channel <- nrow(private$chn)

            ## get the raster data - risky with large maps
            hru_data <- terra::as.matrix(private$brk[[c("filled_dem","channel","channel_fraction",
                                                       "band","gradient")]])
            hru_class <- terra::as.matrix(private$brk[[class_lyr]])
            if(!is.null(rain_lyr)){
                cell_precip <- paste0(rainfall_label,terra::as.matrix(private$brk[[rain_lyr]]))
            }else{
                cell_precip <- rep("precip",nrow(hru_data))
            }
            if(!is.null(pet_lyr)){
                cell_pet <- paste0(pet_label,terra::as.matrix(private$brk[[pet_lyr]]))
            }else{
                cell_pet <- rep("pet",nrow(hru_data))
            }

            ## work out the number of HRUs
            id <- rep(NA,nrow(hru_data))
            idx <- order(hru_data[,"filled_dem"],na.last=NA) ## search order
            nhru <- n_channel + sum( hru_data[idx,"channel_fraction"]!=1 ) ## number of hrus is num of channels + number of cells with hillslopes

            ## construct template for HRU
            tmplate <- list(uid = c(id = NA_integer_, band = NA_integer_, cell = NA_integer_),
                            states = setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz")),
                            properties = setNames(rep(0,3), c("area","Dx","gradient")),
                            sf = list(),
                            rz = list(type="orig", parameters = c("s_rzmax" = 0.1)),
                            uz = list(type="orig", parameters = c("t_d" = 8*60*60)),
                            sz = list(),
                            sf_flow_direction = list(id = integer(0), fraction = numeric(0)),
                            sz_flow_direction = list(id = integer(0), fraction = numeric(0)),
                            initialisation = c("s_rz_0" = 0.75, "r_uz_sz_0" = 1e-7),
                            precip = numeric(0),
                            pet = numeric(0),
                            class = list()
                            )

            tmplate$sf <- switch(sf_opt,
                                 "cnst" = list(type = "cnst",
                                               parameters = c("v_sf" = 0.3,"s_raf" = 0.0, "t_raf" = 999.9)),
                                 "kin" = list(type = "kin",
                                              parameters = c("n" = 0.03,"s_raf" = 0.0, "t_raf" = 999.9)),
                                 stop("Unrecognised surface option")
                                 )

            tmplate$sz <- switch(sz_opt,
                                 "exp" = list(type = "exp",
                                              parameters = c( "t_0" = 0.135, "m" = 0.04 )),
                                 ## "bexp" = list(type = "bexp",
                                 ##               parameters = c( "t_0" = 0.135, "m" = 0.04 , "h_sz_max" = 5)),
                                 "dexp" = list(type = "dexp",
                                               parameters = c( "t_0" = 0.135, "m" = 0.04, "m2" = 0.1, "omega"=0.5)),
                                 ## "cnst" = list(type = "cnst",
                                 ##               parameters = c( "v_sz" = 0.1, "h_sz_max" = 5 )),
                                 stop("Unrecognised saturated zone option")
                                 )
            if(is.null(rain_lyr)){ tmplate$precip <- c("precip"=1) }#list(name="precip", fraction = 1) }
            if(is.null(pet_lyr)){ tmplate$pet <- c("pet"=1) }#list(name = "pet", fraction = 1) }

            ## initalise the hrus
            if(verbose){ cat("Initialise the HRUs","\n") }
            hru <- rep(list(tmplate), nhru )

            ## pass through all the cells...
            if( verbose ){ cat("Passing over cells","\n") }
            cnt <- n_channel
            for(ii in idx){
                chn_frc <- hru_data[ii,"channel_fraction"]
                if( chn_frc == 1 ){ next } ## totally handled in the channel part

                ## process the hillslope part of the cell
                cnt <- cnt + 1 ## get new id
                id[ii] <- cnt

                ##populate the uid
                hru[[cnt]]$uid["id"] <- as.integer(cnt)
                hru[[cnt]]$uid["band"] <- as.integer(hru_data[ii,"band"])
                hru[[cnt]]$uid["cell"] <- as.integer(ii)
                ## add class data
                hru[[cnt]]$class <- as.list(hru_class[ii,])

                ## work out flow direction and associated properties
                jdx <- ii+delta
                grd <- (hru_data[ii,"filled_dem"]-hru_data[jdx,"filled_dem"])/dxy
                gcl <- grd*dcl
                area <- cell_area * (1-chn_frc)
                if( is.finite(hru_data[ii,"channel"]) ){ ## cell is part channel - part hillslope
                    ## work out gradient from cells flowing in
                    is_higher <- is.finite(gcl) & gcl<0 #& is.finite(cjdx) & cjdx==ctch[ii]
                    if( any(is_higher) ){
                        sum_gcl <- sum( gcl[is_higher] )
                        sum_dcl <- sum( dcl[is_higher] )
                        gr <- max(-sum_gcl / sum_dcl,min_grad)
                    }else{ ## nothing flows into the cell
                        gr <- min_grad
                    }
                    sum_dcl <- mean(rs) ## to get correct width
                    hru[[cnt]]$sf_flow_direction = list(id = as.integer(hru_data[ii,"channel"]),
                                                        fraction=1)
                    wdth <- mean(rs)
                }else{
                    ## flow goes to lower hillslopes
                    is_lower <- is.finite(gcl) & gcl>0
                    kk <- jdx[is_lower]
                    stopifnot(
                        "All hillslope cells must flow to those with lower id's" =
                            all( cnt>id[kk] | hru_data[kk,"channel_fraction"]==1 ),
                        "All hillslope cells must flow to those with bands" =
                            all( hru_data[ii,"band"] > hru_data[kk,"band"] ),
                        "Hillslopes must drain down" = length(kk)>0
                    )
                    sum_gcl <- sum( gcl[is_lower] )
                    sum_dcl <- sum( dcl[is_lower] )
                    gr <- max(sum_gcl / sum_dcl,min_grad)
                    ## set flow directions
                    hru[[cnt]]$sf_flow_direction = list(id = as.integer(id[kk]),
                                                        fraction = gcl[is_lower]/sum(gcl[is_lower]))
                    wdth <- sum_dcl
                }
                ## set proerties
                hru[[cnt]]$properties[c("area","Dx","gradient")] <-
                    as.numeric(c(area, area/wdth, gr))

                ## work out precipitation and pet
                kk <- cell_precip[ii]
                hru[[cnt]]$precip[kk] <- 1
                kk <- cell_pet[ii]
                hru[[cnt]]$pet[kk] <- 1


                if( length(hru[[cnt]]$sf_flow_direction)==0 ){ stop("Hillslope HRU with no outflow") }
                if( hru[[cnt]]$properties["area"]==0 ){ stop("Hillslope HRU with no area") }

            }

            if( verbose ){ cat("Processing channel inputs","\n") }
            input_tbl <- terra::extract(private$brk[[c(rain_lyr,pet_lyr)]],private$chn) ## slow ish

            if( verbose ){ cat("Processing channel HRUs","\n") }
            shp <- as.data.frame(private$chn) ## copy channel data since quicker
            chn_class_names <- setdiff(names(shp), c("id","band","length","slope","area")) ## channel class info to copy
            outlets <- list() ## initialise list of outlets

            for(ii in 1:n_channel){
                ## it is a channel HRU
                hru[[ii]]$uid["id"] <- as.integer( shp$id[ii] )
                hru[[ii]]$uid["band"] <- as.integer( shp$band[ii] )
                hru[[ii]]$properties["Dx"] <- as.numeric( shp$length[ii] )
                hru[[ii]]$properties["gradient"] <- as.numeric( shp$slope[ii] )
                hru[[ii]]$properties["area"] <- as.numeric( shp$area[ii] )
                hru[[ii]]$class <- as.list( shp[ii,chn_class_names] )


                tbl <- table(input_tbl[input_tbl$ID==ii,2])
                hru[[ii]]$precip <- setNames(tbl/sum(tbl), paste0(rainfall_label,names(tbl)))
                tbl <- table(input_tbl[input_tbl$ID==ii,3])
                hru[[ii]]$pet <- setNames(tbl/sum(tbl), paste0(pet_label,names(tbl)))

                ## do downstream routing
                kk <- shp$id[ shp$startNode == shp$endNode[ii] ]
                if(length(kk)>0){
                    ## has downstream
                    hru[[ii]]$sf_flow_direction <- list(id = as.integer(kk), fraction = rep(1/length(kk),length(kk)))
                }else{
                    ## is an outlet
                    outlets[[length(outlets)+1]] <- data.frame(name = paste0("q_sf_",ii),
                                                               id = as.integer( ii ),
                                                               flux = "q_sf", scale = 1.0)
                }
            }

            if( verbose ){ cat("Passing through HRUs to sort indexing","\n") }
            for(ii in 1:nhru){
                ## for both - sort out so 0 indexed
                hru[[ii]]$uid["id"] <- as.integer(hru[[ii]]$uid["id"] - 1) ## since 0 indexed in dynatop
                hru[[ii]]$sf_flow_direction$id <- hru[[ii]]$sf_flow_direction$id - 1L
                hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
            }

            ## tidy up outlets
            outlets <- do.call(rbind,outlets)
            outlets$id <- as.integer( outlets$id - 1 )
            outlets$name <- paste0("q_sf_",outlets$id)

            ## correct maps etc to 0 index
            id[hru_data[,"channel_fraction"]==1] <- NA
            id <- id - 1

            ## make output
            if(verbose){ cat("Making output","\n") }
            rst <- terra::rast( private$brk[[ "channel" ]],
                               names = "hru", vals=id)
            model <- list(hru=hru, output_flux = outlets)
            model$map <- paste0(layer_name,".tif")
            terra::writeRaster(rst,model$map,overwrite=TRUE)
            saveRDS(model,paste0(layer_name,".rds"))
        }

    )
    )



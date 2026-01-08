#' R6 Class for processing a catchment to make a dynatop TOPMODEL
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
#' ## load the dem
#' dem_file <- system.file("extdata/gis", "SwindaleDTM40m.tif", package="dynatop", mustWork = TRUE)
#' dem <- terra::rast(dem_file)
#'
#' ## Build catchment outline - would normally provide this
#' catchment_outline <- terra::terra::ifel(is.finite(dem),1,NA)
#' catchment_outline <- as.polygons(catchment_outline, dissolve = TRUE, na.rm = TRUE)
#'
#' ## add the catchment and dem to the project
#' ctch$add_catchment(catchment_outline,dem)
#'
#' ## add channel to the project
#' channel_file <- system.file("extdata/gis", "SwindaleRiverNetwork.shp",
#' package="dynatop", mustWork = TRUE)
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
        #' @param projectFile is a tif file of spatial maps
        #'
        #' @details If it exists this loads the tif file or states a new project is being created.
        #'
        #' @return A new `dynatopGIS` object
        initialize = function(projectFile){

            private$apply_initialize( projectFile )
            invisible(self)
        },
        #' @description Add a catchment outline to the `dynatopGIS` project
        #'
        #' @param catchment a \code{SpatVect} object or the path to file containing one which contains the catchment outline
        #' @param dem a \code{SpatRast} object or the path to file containing one which is the DEM

        #' @details If teh inputs are not spatial objects they are read in using the terra package. The catchment is rasterised to the DEM resolution. If an \code{id} column is present in the catchment this is used so subcatchments can be delineated. The resolution and projection of the project is taken from the provided dem. all NA values in the dem within the catchment are replaced by the \code{fill_na} values so they can be identified as sinks.
        #'
        #' @return \code{invisible(self)}
        add_catchment = function(catchment, dem, fill_na=-9999){
            ## check catchment outline is a vector
            if(!("SpatVector" %in% class(catchment))){ catchment <- terra::vect(as.character(catchment)) }
            if(!("SpatVector" %in% class(catchment))){ stop("catchment is not a SpatVect") }
            ## check dem is a raster
            if(!("SpatRaster" %in% class(dem))){ dem <- terra::rast(as.character(dem)) }
            if(!("SpatRaster" %in% class(dem))){ stop("dem is not a SpatRaster") }

            private$apply_add_catchment(catchment,dem)
            invisible(self)
        },
        #' @description Import channel data to the `dynatopGIS` object
        #'
        #' @param channel a \code{SpatVect} object or file path that can be loaded as one containing the channel information
        #' @param verbose Should additional progress information be printed
        #' @details Takes the representation of the channel network as a SpatVect with properties name, length, area, startNode, endNode and overlaying it on the DEM. In doing this a variable called id is created (or overwritten) other variables in the data frame are passed through unaltered.
        #'
        #' @return suitable for chaining
        add_channel = function(channel,verbose=FALSE){
            if(!is(channel,"SpatVector")){ channel <- terra::vect( as.character(channel) ) }
            if(!is(channel,"SpatVector")){ stop("channel is not a SpatVector object") }
            ## TODO check is polygons

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
            stopifnot("layer is not a SpatRaster" = "SpatRaster" %in% class(layer),
                      "Length of names does not match number of layers" = length(layer_name) == terra::nlyr(layer)
                      )
            private$apply_add_layer(layer,layer_name)
            invisible(self)
        },
        #' @description Get a layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a `raster` layer of the requested information if layer_name is given else a vector of layer names
        get_layer = function(layer_name=NULL){
            ## create a vector of available layers
            tmp <- names(private$brk)
            if( "channel_id" %in% tmp ){ tmp <- c(tmp,"channel") }

            if(is.null(layer_name)){ return(tmp) }

            ## check layer name exists
            layer_name <- match.arg(layer_name,tmp,several.ok=TRUE)

            ## make raster and return
            if( "channel" %in% layer_name){
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
                if(layer_name %in% names(private$chn)){
                    terra::plot(private$chn,layer_name, add=TRUE )
                }else{
                    terra::plot(private$chn, add=TRUE )
                }
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
        sink_fill = function(min_grad = 1e-4,max_it=1e6,verbose=FALSE, hot_start=FALSE){
            private$apply_sink_fill(min_grad,max_it,verbose,hot_start)
            private$apply_downward_pass(verbose)
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
        ## compute_properties = function(min_grad = 1e-4,verbose=FALSE){
        ##     private$apply_upward_pass(verbose)
        ##     private$apply_downward_pass(min_grad,verbose)
        ##     invisible(self)
        ## },
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
        ## classify = function(layer_name,base_layer,cuts){
        ##     private$apply_classify(as.character(base_layer[1]), cuts, as.character(layer_name[1]))
        ##     ##invisible(self)
        ## },
        #' @description Combine any number of classifications based on unique combinations and burns
        #' @param layer_name name of the new layer to create
        #' @param pairs a vector of layer names to combine into new classes through unique combinations. Names should correspond to raster layers in the project directory.
        #' @param burns a vector of layer names which are to be burnt on
        #'
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given, only positive burn values are applied.
        ## combine_classes = function(layer_name,pairs,burns=NULL){
        ##     private$apply_combine_classes(as.character(layer_name[1]),
        ##                                   as.character(pairs),
        ##                                   as.character(burns))
        ## },
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
        create_model = function(model_name,class_layer = NULL,
                                rain_layer=NULL,# rain_label="precip",
                                pet_layer=NULL,# pet_label="pet",
                                verbose=FALSE){

            ## check valid transmissivity and channel_solver
            ##sf_opt<- match.arg(sf_opt)
            ##sz_opt <- match.arg(sz_opt)
            private$apply_create_model(model_name,class_layer,
                                       rain_layer,# rain_label,
                                       pet_layer,# pet_label,
                                       verbose
                                       )

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
        version = "0.5.0",
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

            ## check is padded
            isPadded <- !(any( is.finite(brk[1,][[1]]) ) |
                          any( is.finite(brk[nrow(brk),][[1]]) ) |
                          any( is.finite(brk[,1][[1]]) ) |
                          any( is.finite(brk[,ncol(brk)][[1]]) ))

            ## initial error checks
            stopifnot("Raster data not valid: Processing currently only works on projected data with a square grid" =
                          all.equal(diff(terra::res(brk)),0) & !terra::is.lonlat(brk),
                      "No catchment layer in Raster data" = "catchment" %in% names(brk),
                      "No DEM layer in Raster data" = "dem" %in% names(brk),
                      "Raster data not valid: Processing currently only works with padded grids" = isPadded)

            if( "channel_id" %in% names(brk) ){
                channelFile <- paste0(tools::file_path_sans_ext(projectFile),".geojson")
                stopifnot( "No channel file available" = file.exists(channelFile) )
                shp <- terra::vect( channelFile )
            }else{
                shp <- NULL
            }

            private$projectFile <- projectFile
            private$brk <- brk
            private$chn <- shp
        },
        ## add a catchment map
        apply_add_catchment = function(catchment,dem){

            ## inital stops
            stopifnot(
                "The catchment map already exists, start a new project" = !("catchment" %in% names(private$brk)),
                "DEM data not valid: Processing currently only works on projected data with a square grid" =
                    all.equal(diff(terra::res(dem)),0) & !terra::is.lonlat(dem),
                "The DEM must be one layer" = terra::nlyr(dem)==1
            )

            ## reproject the catchment
            catchment <- terra::project(catchment,dem) ## convert catchment to dem projection

            ## check the extents
            extDiff <- ( as.vector(terra::ext(dem)) -
                         as.vector(terra::ext(catchment)) ) * c(-1,1,-1,1)## xmin, xmax, ymin, ymax
            stopifnot("Catchment extends beyound raster" = all(extDiff>=0))

            ## crop the dem to the catchment
            dem <- terra::crop(dem,catchment,mask=TRUE)
            dem <- terra::extend(dem,1) ## pad with NA
            names(dem) <- "dem"

            ## rasterise the catchment
            if("id" %in% names(catchment)){
                ctch <- terra::rasterize(catchment,dem,field="id")
            }else{
                ctch <- terra::rasterize(catchment,dem,field=1L)
            }
            names(ctch) <- "catchment"

            ## check is padded - should not be needed but safer to check again..
            isPadded <- !(any( is.finite(ctch[1,][[1]]) ) |
                          any( is.finite(ctch[nrow(ctch),][[1]]) ) |
                          any( is.finite(ctch[,1][[1]]) ) |
                          any( is.finite(ctch[,ncol(ctch)][[1]]) ))

            ## save
            private$brk <- c(ctch,dem)
            private$save_project()
        },
        ## add the channel
        apply_add_channel = function(chn,verbose){
            rq <- c("catchment","dem")

            stopifnot(
                "Not all required layers are available" = all(rq %in% names(private$brk)),
                "Channel already added" = !("channel" %in% names(private$brk)),
                "Projection of channel object does not match that of project" =
                    terra::crs(private$brk, proj=TRUE) == terra::crs(chn, proj=TRUE)
            )

            if(verbose){ print("Checking channel") }
            check_channel(chn)

            ## pass down channel adding id and bands
            if("id" %in% names(chn)){ warning("Overwriting id variable") }
            chn$id <- 1:nrow(chn)

            ## channel raster of id
            chn_rst <- terra::rasterize(chn,private$brk[["catchment"]],field = "id",touches=TRUE)
            names(chn_rst) <- "channel_id"
            ## create a raster of channel coverage fractions
            chn_frac <- terra::rasterize(chn,private$brk[["catchment"]],background=0,cover=TRUE) ## fraction of cell covered by channel
            chn_frac <- terra::mask(chn_frac,private$brk[["catchment"]])
            names(chn_frac) <- "channel_fraction"
            terra::values(chn_frac) <- round(terra::values(chn_frac),2) ## else get horrible rounding errors close to 1
            chn_frac[chn_frac==0 & !is.na(chn_rst)] <- 0.005 ## add a fraction to those cells with an ID but no fractions

            ## save output
            private$brk <- c(private$brk,chn_rst,chn_frac) #,chn_depth)
            private$chn <- chn
            private$save_project(chn=TRUE)

            ## TODO remove channels that aren't in the raster - since these will be small and cause grief
        },
        ## Add a layer
        apply_add_layer=function(layer,layer_name){

            rq <- c("channel_id","catchment")
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
                                   is.na(private$brk[["channel_id"]]) ,max)>0 ){
                    stop(paste("Layer",ii,"has missing values in the catchment - try running fill_na first"))
                }
            }

            ## save output
            private$brk <- c(private$brk,layer)
            private$save_project()
        },
        ## Sink fill
        apply_sink_fill = function(min_grad,max_it,verbose,hot_start){ #,flow_type){

            ## recall the catchments is padded with NA values
            if(!hot_start & ("filled_dem" %in% names(private$brk))){
                warning("Filled DEM already created")
                return()
            }

            d <- ifelse(hot_start,"filled_dem","dem")

            rq <- c(d,"channel_id","catchment")

            if(!all(rq %in% names(private$brk))){
                stop("Not all required layers are available")
            }


            d <- terra::as.matrix( private$brk[[d]] ,wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel_id"]] , wide=TRUE )
            ctch <- terra::as.matrix( private$brk[["catchment"]] , wide=TRUE )

            ## values that should be valid
            to_be_valid <- !is.na(ctch) # all values not NA should have a valid height
            is_valid <- !is.na(ch) # TRUE if a channel cell for initialisation
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
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        to_eval[jdx] <- !is_valid[jdx]
                    }
                }
                to_eval <- to_eval & to_be_valid
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
                        ## adjust mind
                        mind <- pmin(mind,fd[jdx] + dxy[ii+2,jj+2],na.rm=TRUE)
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
        ## function to do property calculations on an upwards pass (high to low DEM values)
        apply_downward_pass = function(verbose){

            rq <- c("filled_dem","channel_id","channel_fraction")
            stopifnot(
                "Not all required input layers have been generated \n Try running sink_fill first" =
                    all( rq %in% names( private$brk) )
            )

            ## load dem
            d <- terra::values( private$brk[["filled_dem"]] )
            chn_id <- terra::values( private$brk[["channel_id"]] )
            chn_frc <- terra::values( private$brk[["channel_fraction"]] )

            ## initialise the bands
            bnd <- is.finite(d)*1 ## initialise everything
            chn_bnd <- rep(1,nrow(private$chn)) ## we will populate this by id so will need merging back into private$chn carefully

            if( verbose ){ print("Computing downward pass of hillslope") }

            idx <- order(d,na.last=NA,decreasing=TRUE) ## search order

            nr <- terra::ncol(private$brk); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1) ## neighbours

            ## set up printing variables
            if(verbose){
                print_step <- c(1,rep(round(length(idx)/20,2)),length(idx)) # current, next print, step, total
            }

            ## main loop through hillslope cells
            dz <- rep(NA,8)
            for(ii in idx){
                jj <- chn_id[ii]
                if( is.na(jj) ){
                    jdx <- ii+delta
                    dz[] <- d[ii] - d[jdx]
                    is_lower <- is.finite(dz) & dz>0
                    bnd[ jdx[is_lower] ] <-  pmax( bnd[ jdx[is_lower] ], bnd[ii]+1 )
                }else{
                    chn_bnd[jj] <- max( bnd[ii], chn_bnd[jj] )
                    if(chn_frc[ii]==1){bnd[ii] <- NA} ## this is important else these cells get an HRU id
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

            if( verbose ){ print("Computing downward pass of channels") }
            ## assume channel id values are ordered 1:nrow(chn)
            stopifnot("Channel id are not correct" = all( 1:nrow(chn) == private$chn$id ))
            ## This is much quicker using vectors and not constantly accessing via the SpatVect object
            sN <- private$chn$startNode
            eN <- private$chn$endNode

            idx <- !(sN %in%eN) ## initial channel reaches
            it <- 0
            cnt <- table(eN) ## we should never vist a node more times then it is an endpoint
            while(sum(idx)>0){
                chn_bnd[idx] <- pmax( chn_bnd[idx], it )
                jdx <- eN[idx]
                for(ii in jdx){ cnt[ii] <- cnt[ii] - 1 } ## since eN might appear more then once..
                if( any(cnt<0) ){
                    stop(paste("Failing loop involving nodes:", paste(names(cnt)[cnt<0],collapse=", ")))
                }
                jdx <- jdx[cnt[jdx]==0] ## only move down if it is the last visit to the endNode
                idx <- sN %in% jdx

                it <- it+1
            }
            stopifnot(
                "Error processing channel: problem with visiting all points" = all(cnt==0)
            )
            private$chn$band <- chn_bnd

            if( verbose ){ print("Computing HRU IDs") }
            ## work out hru id
            max_hru <- -1
            hru <- rep(NA,length(bnd))
            chn_hru <- rep(NA,length(chn_bnd))
            for(ii in 1:max(chn_bnd)){
                idx <- bnd==ii
                nidx <- sum(idx)
                hru[idx] <- max_hru + 1:nidx
                max_hru <- max_hru + nidx
                idx <- chn_bnd == ii
                nidx <- sum(idx)
                chn_hru[idx] <- max_hru + 1:nidx
            }
            private$chn$hru <- chn_hru

            ## save
            rbnd <- private$brk["filled_dem"]
            names(rbnd) <- "band"
            terra::values(rbnd) <- bnd
            rhru <- private$brk["filled_dem"]
            names(rhru) <- "hru"
            terra::values(rhru) <- hru
            private$brk <- c(private$brk,rbnd,rhru)
            private$save_project(chn=TRUE)

        },
        ## create a model
        apply_create_model = function(model_name,class_lyr,
                                      rain_lyr,
                                      pet_lyr,
                                      verbose){

            ## check layers
            rq <- c("hru","band",
                    "filled_dem",
                    "channel_id","channel_fraction",
                    rain_lyr,
                    pet_lyr,
                    class_lyr)

            stopifnot(
                "Missing layers" = all(rq %in% names(private$brk))
            )

            ## work out some properties of the brick
            rs <- terra::res( private$brk )[1]
            nr <- terra::ncol(private$brk)
            delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            sc <- c(sqrt(2),1,sqrt(2),1,1,sqrt(2),1,sqrt(2))
            dxy <- sc*rs
            dcl <- (0.5/sc)*rs
            cell_area <- rs*rs

            ## read in the data required for hillslope
            mdl <- terra::values(private$brk[[rq]], dataframe=TRUE)

            ## rename and sort out some variables
            if(is.null(rain_lyr)){
                mdl$precip <- 2L
            }else{
                mdl$precip <- mdl[[rain_lyr]]
                mdl[[rain_lyr]] <- NULL
            }

            if(is.null(pet_lyr)){
                mdl$pet <- 3L
            }else{
                mdl$pet <- mdl[[pet_lyr]]
                mdl[[pet_lyr]] <- NULL
            }

            mdl$z <- mdl$filled_dem
            mdl$filled_dem <- NULL

            mdl$is_channel <- FALSE
            mdl$edges <- NA_character_

            ## compute area and set value to NA if no hillslope
            mdl$area <- cell_area * (1 - mdl$channel_fraction)

            ## initialise the channel HRUs
            chn <- as.data.frame(private$chn)
            chn$is_channel <- TRUE
            chn$edges <- NA_character_
            chn$z <- NA

            ## compute the index of cellls to evaluate
            cell_idx <- which( is.finite(mdl$z) )

            ## initialise the channel storage to compute as we loop
            chn_n <- nrow(chn)
            total_chn_frac <- rep(0,chn_n)
            chn_precip <- rep(list(NULL),chn_n)
            chn_pet <- rep(list(NULL),chn_n)

            ## Loop hillslope cells
            slp <- rep(NA,8)
            for(ii in cell_idx){
                jj <- mdl$channel_id[ii]
                if( is.finite(jj) ){
                    total_chn_frac[jj] <- total_chn_frac[jj] + mdl$channel_fraction[ii]
                    ## edges
                    mdl$edges[ii] <- sprintf('{"hru":%i,"width":%f,slope:%f}',
                                             chn$hru[jj], rs, chn$depth[jj] / rs)
                    ## add precip
                    str <- paste(mdl$precip[ii])
                    if( !(str %in% names(chn_precip[[jj]])) ){ chn_precip[[jj]][str] <- 0 }
                    chn_precip[[jj]][str] <- chn_precip[[jj]][str] + mdl$channel_fraction[ii]
                    ## add pet
                    str <- paste(mdl$pet[ii])
                    if( !(str %in% names(chn_pet[[jj]])) ){ chn_pet[[jj]][str] <- 0 }
                    chn_pet[[jj]][str] <- chn_pet[[jj]][str] + mdl$channel_fraction[ii]

                }else{
                    jdx <- ii+delta
                    slp[] <- (mdl$z[ii] - mdl$z[jdx])/dxy
                    is_lower <- is.finite(slp) & (slp > 0)
                    if(!any(is_lower)){
                        stop(paste("No lower cells for",ii))
                    }
                    mdl$edges[ii] <- sprintf('{"hru":%s,"width":%s,slope:%s}',
                                             paste0("[", paste(mdl$hru[jdx[is_lower]],collapse=","),"]"),
                                             paste0("[", paste(slp[is_lower],collapse=","),"]"),
                                             paste0("[", paste(dcl[is_lower],collapse=","),"]"))
                }
            }

            ## apply some more channel HRU calculations
            chn$area <- total_chn_frac * cell_area
            fin <- function(x){
                out <- if(length(x)==0){ NA_integer_ }else{ as.integer(names(x)[which.max(x)]) }
                return(out)
            }
            chn$precip <- sapply(chn_precip,fin)
            chn$pet <- sapply(chn_pet,fin)

            ## loop channels to make edges
            outlets <- NULL
            for(ii in 1:nrow(chn)){
                ## edges
                idx <- chn$startNode == chn$endNode[ii]
                if(any(idx)){
                    ne <- sum(idx)
                    chn$edges[ii] <- sprintf('{"hru":%s,"width":%s,slope:%s}',
                                             paste0("[", paste(chn$hru[idx],collapse=","),"]"),
                                             paste0("[", paste(rep(chn$slope[ii],ne), collapse=","),"]"),
                                             paste0("[", paste(rep(chn$width[ii]/ne,ne), collapse=","),"]"))
                }else{
                    outlets <- c(outlets,chn$hru[ii])
                }
            }

            ## add geometries and parameters to the channel then remove unwanted variables
            chn$geom <- terra::geom(private$chn,wkt=TRUE)
            chn$sf_type <- 2; chn$sf_param <- paste0('["n":0.1,"depth":',chn$depth,']')
            chn$rz_type <- 1; chn$rz_param <- '["s_rz_max":0.1]'
            chn$uz_type <- 1; chn$uz_param <- '["t_uz":0.001]'
            chn$sz_type <- 1; chn$sz_param <- '["t_0":0.0,"m":0.02]'
            chn[,c("startNode","endNode","depth","slope","width","id")] <- NULL

            mdl <- mdl[is.finite(mdl$hru),]
            mdl$geom <- terra::geom(terra::as.polygons(private$brk[["hru"]],aggregate=FALSE,values=FALSE),
                                    wkt=TRUE)
            mdl$sf_type <- 1; mdl$sf_param <- '["n":0.3]'
            mdl$rz_type <- 1; mdl$rz_param <- '["s_rz_max":0.1]'
            mdl$uz_type <- 1; mdl$uz_param <- '["t_uz":0.001]'
            mdl$sz_type <- 1; mdl$sz_param <- '["t_0":0.001,"m":0.02]'
            mdl$channel_id <- mdl$channel_fraction <- NULL

            mdl <- merge(mdl,chn,all=TRUE)
            mdl$s_sf <- NA
            mdl$s_rz <- NA
            mdl$s_uz <- NA
            mdl$s_sz <- NA

            ## sort outlets
            output <- data.frame(names = paste0("q_sf_",outlets),
                                 hru = outlets,
                                 flux = "q_sf",
                                 scale = 1)

            ## save the model
            write.table(mdl, paste0(model_name,".csv"),
                        sep=";",row.names=FALSE,
                        na="",quote=FALSE,
                        fileEncoding = "UTF-8")
            write.csv(mdl,paste0(model_name,"_outlets.csv"),row.names=FALSE)

            saveRDS(list(hru=mdl, output = output),paste0(model_name,".rds"))

            mdl <- terra::vect(mdl,geom="geom")
            terra::writeVector(mdl,paste0(model_name,".geojson"),overwrite=TRUE)
        }
    )
    )

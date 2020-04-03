#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
dynatopGIS <- R6::R6Class(
    "dynatop",
    public = list(
        #' @description Sets up the \code{dynatop} object from a model description
        #'
        #' @param model a model objectRasterLayer containing the dem
        #' @param delta error term in checking redistribution sums
        #'
        #' @details This function makes some basic consistency checks on a list representing a dynamic TOPMODEL model.
        #' @return A new `dynatop` object
        initialize = function(mdl, verbose=FALSE, delta = 1e-13){

            private$check_model(mdl,verbose,delta)
            private$create_prop()
            ## return so can be chained
            invisible(self)
        },
    ),
    private = list(
        ## stores of data
        version = 0.1,
        scope = list(components <- c("hillslope","channel","param",
                                     "gauge","point_inflow")),
        model = list() # store for the model
        prop = list() # store for the properties extracted from the model and states
        states = list() # store for the states as saved
        input = list() # store for the input series
        output = list() # store for the output series

        ## functions
        # Function to check the model to be used in a dynamic TOPMODEL run
        check_model = function(model, verbose=FALSE, delta=1e-13){

            ## check all components of the model exist
            components <- private$scope$components
            idx <- components %in% names(model)
            if( !all(idx) ){
                stop(paste("Missing componets:",paste(components[!idx],collapse=",")))
            }

            ## check components that should be data.frames of given structure

            ## check the HRU table properties
            req_names <- list(output_names = list(),
                              parameter = list(),
                              data_series = list())
            for(ii in setdiff(components,"param")){
                ## what should the properties of each column be
                prop <- private$definition(ii)
                
                if(!is.data.frame(model[[ii]])){
                    stop(paste("Table",ii,"should be a data.frame"))
                }

                idx <- df_prop[[ii]]$name %in% names(model[[ii]])

                if( !all( idx ) ){# check it has required columns
                    stop( paste("Table",ii,"is missing columns:",
                                paste(df_prop[[ii]]$name[!idx],collapse=",")) )
                }

                ## check data types
                tmp <- sapply(model[[ii]],class) # types of the columns
                idx <- tmp[ df_prop[[ii]]$name ] != df_prop[[ii]]$type
                if( any( idx ) ){
                    stop( paste("Incorrect types in table",ii,"columns:",
                                paste(df_prop[[ii]]$name[idx],collapse=",")) )
                }

                ## take the required names
                for(jj in names(req_names)){
                    tmp <- prop$name[prop$role==jj]
                    req_names[[jj]][[ii]] <- unlist(model[[ii]][,tmp])
                }
            }

            ## unpack the required names to vectors
            for(jj in names(req_names)){
                req_names[[jj]] <- do.call(c,req_names[[jj]])
            }
            

            ## parameter vector should be named numeric vector and contain all required names
            if( !all(is.vector(model$param), is.numeric(model$param)) ){
                stop("param should be a numeric vector")
            }
            if( length(unique(names(model$param))) != length(model$param) ){
                stop("All values in param should have a unique name")
            }
            idx  <- req_names$parameter %in% names(model$param)
            if(!all(idx)){
                stop(paste("The following parameters are not specified:",
                           paste(req_names[!idx],collapse=",")))
            }
            idx  <- names(model$param) %in% req_names$parameter
            if(!all(idx)){
                stop(paste("The following parameters are not used:",
                           paste(names(model$param)[!idx],collapse=",")))
            }

            ## check all output series have unique names
            if( length(req_names$output_names) != length(unique(req_names$output_names)) ){
                stop("All output series should have a unique name")
            }

            ## checks on hillslope and channel HSU ids
            all_hsu <- c(model$hillslope$id,model$channel$id)
            if( length(all_hsu) != length(unique(all_hsu)) ){
                stop("HSU id values should be unique") }
            if( !all(is.finite(all_hsu)) ){ stop("HSU id values should be finite") }
            if( !all(range(all_hsu)==c(1,length(tmo))) ){
                stop("HSU id's should be numbered consecuativly from 1")
            }
            
            ## all points_inflows and gauges should be on a channel
            ## with fractions between 0 & 1
            for(jj in c("gauge","point_inflow")){
                if(nrow(model[[jj]]) == 0){next}
                idx <- (model[[jj]]$id %in% model[['channel']]$id) &
                    (model[[jj]]$fraction >= 0) &
                    (model[[jj]]$fraction <= 1)
                if( any(!idx) ){
                    stop(paste("The following", ii , "are incorrectly specified:",
                               paste(model[[jj]]$name[!idx],collapse=" ")))
                }
            }

            
            ## checks on redistribution
            fcheck <- function(x){
                all(x$idx %in% all_hsu) & (abs(sum(x$frc)-1) < delta)
            }
            for(jj in c("hillslope","channel")){
                idx <- sapply(model[[jj]]$flow_direction,fcheck)
                if( any(!idx) ){
                    stop(paste("Flow redistribution is not valid for HSUs:",
                               paste(model[[jj]]$id[!idx],sep=" ")))
                }
            }
                 
 
            ## specific checks on channel network connectivity
            chn_con <- lapply(model$channel$flow_direction,function(x){x$idx})
            if( any(sapply(chn_con,length)>1) ){
                stop("Only channels routing to single HSUs are supported")
            }
            chn_con <- do.call(c,chn_con)
            chn_con <- chn_con[!is.na(chn_con)]
            if( !all(chn_con %in% model$channel$id) ){
                stop("Channels routing to non channel HSUs, set next_id to NA to represent an outflow")
            }
            
            to_outlet <- rep(FALSE,length(model$channel$id))
            to_outlet[is.na(chn_con)] <- TRUE # set outlets to true
            ## loop channels at top of network
            for(ii in setdiff(model$channel$id,chn_con)){
                ## set up a record of place in search down tree
                in_search <- rep(FALSE,length(model$channel$id))
                jj <- ii
                while( !in_search[model$channel$id==jj] & # fails if loop
                       !to_outlet[model$channel$id==jj] ){ # fails at outlet
                           in_search[model$channel$id==jj] <- TRUE
                           jj <- chn_con[model$channel$id==jj]
                       }
                if(to_outlet[model$channel$id==jj]){
                    to_outlet[in_search] <- TRUE
                }
            }
            if( any(!to_outlet) ){
                stop(paste("The following channels do not drain to an outlet:",
                           paste(model$channel$id[!to_outlet],collapse=" ")))
            }
            
            ## verbose printing of head and tail channels
            if(verbose){
                ## print out head channels
                message(paste("The head channels are:",
                              paste(setdiff(model$channel$id,chn_con),
                                    collapse=", "),
                              sep="\n"))
                ## print out tail channels
                message(paste("The channels with outfalls:",
                              paste(model$channel$id[is.na(chn_con)]
                                    collapse=", "),
                              sep="\n"))
            }

            ## if here we have passed all test so add model
        },
        ## function to mutate from model to internal structure
        create_prop = function(){
            components <- private$scope$components
            for(jj in setdiff(components,"param")){
                ## what should the properties of each column be
                prop <- private$definition(jj,TRUE,TRUE)
                out <- list()

                for(ii in 1:nrow(prop)){
                    nm <- prop$name[ii]
                    out[[ nm ]] <- switch(
                        prop$role[ii],
                        attribute = unname( provate$model[[tbl]][[nm]] ),
                        data_series = unname( provate$model[[tbl]][[nm]] ),
                        parameter = unname( provate$model$param[ provate$model[[tbl]][[nm]] ]),
                        state = rep(NA,nrow(provate$model[[tbl]])),
                        tmp = rep(NA, nrow(provate$model[[tbl]]))
                    )
                    
                }
                private$prop[[jj]] <- out
            }
        },
        ## Function for defining the elements of a model object
        ##
        ## @description The function returns a description, in terms of variable names, names and role, for the selected table within a Dynamic TOPMODEL object.
        ##
        ## @param type the names of the data frame within the Dynamic TOPMODEL object
        ## @param include_states shoould states be returned
        ## @param include_tmp should tempory variables used in simulation be returned
        ##
        ## @return A data frame descriping the properties of one of the tables found within the model
        ##
        model_description = function(type=c("hillslope","channel","point_inflow","gauge"),
                                     include_states=FALSE,include_tmp=FALSE){
            type <- match.arg(type)
            
            if(type=="hillslope"){
                out <- data.frame(name = c("id","atb_bar","s_bar","area","delta_x","sz_band","sf_band", # attributes associated with catchment HSU
                                           "precip","pet", # names of input series
                                           "q_sfmax","s_rzmax","s_rz0","ln_t0","m","t_d","t_sf", # parameter names
                                           "s_sf","s_rz","s_uz","s_sz","sum_l_sz_in","l_sz","l_szmax", # states
                                           "p","e_p","e_t","l_sf","q_sf_rz","q_rz_sf","q_rz_uz","q_uz_sz","q_uz_sf","q_sz_sf","e_t","l_sf","sum_l_sz_in_t","l_sz_t","Q_minus_t","Q_plus_t","Q_minus_tDt"), ## tempory stores not needed for next timestep
                                  role = c(rep("attribute",7),
                                           rep("data_series",2),
                                           rep("parameter",7),
                                           rep("state",7),
                                           rep("tmp",17)),
                                  type = c("integer",rep("numeric",6),
                                           rep("character",2),
                                           rep("character",7),
                                           rep("numeric",7),
                                           rep("numeric",17)
                                           ),
                                  stringsAsFactors=FALSE)
                
            }
            if(type=="channel"){
                out <- data.frame(name= c("id","area","sz_band","sf_band","length","next_id", # states
                                          "precip","pet", # inputs
                                          "v_ch", # parameters
                                          "sum_l_sz_in", #state
                                          "p","e_p","l_sf","s_ch","sum_l_sz_in_t"),
                                  role = c(rep("attribute",6),
                                           rep("data_series",2),
                                           rep("parameter",1),
                                           rep("state",1),
                                           rep("tmp",5)),
                                  type = c("integer",rep("numeric",4),"integer",
                                           rep("character",2),
                                           rep("character",1),
                                           rep("numeric",1),
                                           rep("numeric",5)),
                                  stringsAsFactors=FALSE)
            }

            if(type=="point_inflow"){
                out <- data.frame(
                    name = c("name","id","fraction"),
                    type=c("character","integer","numeric"),
                    role = c("data_series",rep("property",2)),
                    stringsAsFactors=FALSE)
            }
            
            if(type=="gauge"){
                out <- data.frame(
                    name = c("name","id","fraction"),
                    type=c("character","integer","numeric"),
                    role = c("output_label",rep("property",2)),
                    stringsAsFactors=FALSE)
            }
            
            if(!include_states){
                out <- out[out$role!="state",,drop=FALSE]
            }
            if(!include_tmp){
                out <- out[out$role!="tmp",,drop=FALSE]
            }
            
            return(out)
        },
        
        
        
        
    )
)

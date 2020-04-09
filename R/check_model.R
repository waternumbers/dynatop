#' Function to check the model to be used in a dynamic TOPMODEL run
#'
#' @description This function makes some basic consistency checks on a list representing a dynamic TOPMODEL model.
#'
#' @param model a dynamic TOPMODEL list object
#' @param use_states logical if states are to be checked
#' @param verbose if set prints out further information
#' @param delta error term in checking redistribution matrix sums
#'
#' @return a vector of the names of the expected input series
#'
#' @details The checks performed and basic 'sanity' checks. They do not check for the logic of the parameter values nor the consistncy of states and parameters. Sums of the redistribution matrices are checked to be less then 1+delta.
#' @export
check_model <- function(model, verbose=FALSE, use_states=FALSE,delta=1e-13){

    ## check all components of the model exist
    components <- c("hillslope","channel","param","gauge","point_inflow")
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
        prop <- model_description(ii)
        
        if(!is.data.frame(model[[ii]])){
            stop(paste("Table",ii,"should be a data.frame"))
        }
        
        idx <- prop$name %in% names(model[[ii]])
        
        if( !all( idx ) ){# check it has required columns
            stop( paste("Table",ii,"is missing columns:",
                        paste(prop$name[!idx],collapse=",")) )
        }
        
        ## check data types
        
        tmp <- sapply(model[[ii]],class) # types of the columns
        idx <- tmp[ prop$name ] != prop$type
        if( any( idx ) ){
            stop( paste("Incorrect types in table",ii,"columns:",
                        paste(prop$name[idx],collapse=",")) )
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
    if( !all(range(all_hsu)==c(1,length(all_hsu))) ){
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
    ## TODO ad check that going down band?
    ## TODO add check on bound parameter
    fcheck <- function(x){
        all(x$idx %in% all_hsu) & (abs(sum(x$frc)-1) < delta)
    }
    idx <- sapply(model$hillslope$sz_dir,fcheck)
    if( any(!idx) ){
            stop(paste("Saturated flow redistribution is not valid for HSUs:",
                       paste(model$hillslope$id[!idx],collapse=" ")))
    }
    idx <- sapply(model$hillslope$sf_dir,fcheck)
    if( any(!idx) ){
        stop(paste("Surface flow redistribution is not valid for HSUs:",
                   paste(model$hillslope$id[!idx],collapse=" ")))
    }
    idx <- sapply(model$channel$flow_dir,fcheck)
    if( any(!idx) ){
        stop(paste("Channel flow redistribution is not valid for HSUs:",
                   paste(model$channel$id[!idx],collapse=" ")))
    }
    
    ## specific checks on channel network connectivity
    chn_con <- lapply(model$channel$flow_dir,function(x){x$idx})
    
    if( any(sapply(chn_con,length)>1) ){
        stop("Only channels routing to single HSUs are supported")
    }
    is_outlet <- sapply(chn_con,is.null) # identify outlets
    if( !all(do.call(c,chn_con) %in% model$channel$id) ){
        stop("Channels routing to non channel HSUs, set next_id to NA to represent an outflow")
    }

    
    to_outlet <- is_outlet
    ## loop channels at top of network
    for(ii in setdiff(model$channel$id,do.call(c,chn_con))){
        ## set up a record of place in search down tree
        in_search <- rep(FALSE,length(model$channel$id))
        jj <- which(model$channel$id==ii)
        while( !in_search[jj] & # fails if loop
               !to_outlet[jj] ){ # fails at outlet
                   in_search[jj] <- TRUE
                   jj <- which(model$channel$id==chn_con[[jj]])
               }
        if(to_outlet[jj]){
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
                      paste(model$channel$id[is.na(chn_con)],
                            collapse=", "),
                      sep="\n"))
    }
    
    ## if here we have passed all tests
    return(unique(req_names[["data_series"]]))
}

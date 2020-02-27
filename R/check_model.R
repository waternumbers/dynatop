#' Function to check the model to be used in a dynamic TOPMODEL run
#'
#' @description This function makes some basic consistency checks on a list representing a dynamic TOPMODEL model.
#'
#' @param model a dynamic TOPMODEL list object
#' @param verbose if set prints out further information
#' @param delta error term in checking redistribution matrix sums
#' 
#' @return a vector of the names of the expected input series
#'
#' @details The checks performed and basic 'sanity' checks. They do not check for the logic of the parameter values nor the consistncy of states and parameters. Sums of the redistribution matrices are checked to be less then 1+delta.
#' @export
check_model <- function(model, verbose=FALSE, delta=1e-13){

    ## check all components of the model exist
    components <- c("hillslope","channel","param","Fsf","Fsz","gauge","point_inflow")
    idx <- components %in% names(model)
    if( !all(idx) ){
        stop(paste("Missing componets:",paste(components[!idx],collapse=",")))
    }

    ## check components that should be data.frames
    
    ## create a list of data.frames describing the data frames in the model
    df_prop <- list(
        hillslope = data.frame(name = c("id","area","s_bar","delta_x",
                                        "precip","pet",
                                        "q_sfmax","s_rzmax","s_rz0","ln_t0","m","t_d","t_sf"),
                               type=c(rep("numeric",4),rep("character",9)),
                               role = c(rep("property",4),rep("data_series",2),rep("parameter",7)),
                               stringsAsFactors=FALSE),
        channel = data.frame(name = c("id","area","length","next_id","precip","pet","v_ch"),
                             type=c(rep("numeric",4),rep("character",3)),
                             role = c(rep("property",4),rep("data_series",2),"parameter"),
                             stringsAsFactors=FALSE),
        point_inflow = data.frame(
            name = c("name","id","fraction"),
            type=c("character",rep("numeric",2)),
            role = c("data_series",rep("property",2)),
            stringsAsFactors=FALSE),
        gauge = data.frame(
            name = c("name","id","fraction"),
            type=c("character",rep("numeric",2)),
            role = c("output_label",rep("property",2)),
            stringsAsFactors=FALSE)
    )
    
    ## check the HRU table properties
    for(ii in names(df_prop)){
        if(!is.data.frame(model[[ii]])){
            stop(paste("Table",ii,"should be a data.frame"))
        }
        
        idx <- df_prop[[ii]]$name %in% names(model[[ii]])
        
        if( !all( idx ) ){# check it has required columns
            stop( paste("Table",ii,"is missing columns:",
                        paste(df_prop[[ii]]$name[!idx],collapse=",")) )
        }else{ ## check data types
            tmp <- sapply(model[[ii]],class) # types of the columns
            idx <- tmp[ df_prop[[ii]]$names ] != df_prop[[ii]]$type
            if( any( idx ) ){
                stop( paste("Incorrect types in table",ii,"columns:",
                            paste(df_prop[[ii]]$name[idx],collapse=",")) )
            }
        }
    }

    ## parameter vector should be named numeric vector and contain all required names
    if( !all(is.vector(model$param), is.numeric(model$param)) ){
        stop("param should be a numeric vector")
    }
    if( length(unique(names(model$param))) != length(model$param) ){
        stop("All values in param should have a unique name")
    }

    ## check paraemter names
    req_names <- NULL
    for(jj in names(df_prop)){
        tmp <- df_prop[[jj]]$name[df_prop[[jj]]$role=="parameter"]
        tmp <- unique(unlist(model[[jj]][,tmp]))
        req_names <- unique(c(req_names,tmp))
    }
    
    idx  <- req_names %in% names(model$param)
    if(!all(idx)){ stop(paste("The following parameters are not specified:",paste(req_names[!idx],collapse=","))) }
    idx  <- names(model$param) %in% req_names
    if(!all(idx)){ stop(paste("The following parameters are not used:",paste(names(model$param)[!idx],collapse=","))) }

    ## ids should be numeric, sequential and start from 1
    ids <- c(model[['hillslope']]$id,model[['channel']]$id)
    rng_ids <- range(ids)
    if( length(ids) != length(unique(ids)) ){ stop("id should be unique") }
    if( !all(is.finite(ids)) ){ stop("id values should be finite") }
    if( !all(range(ids)==c(1,length(ids))) ){ stop("id's should be numbered consecuativly from 1") }

    ## all points_inflows and gauges should be on a channel with fractions between 0 & 1
    for(jj in c("gauge","point_inflow")){
        if(nrow(model[[jj]]) == 0){next}
        
        for(ii in 1:nrow(model[[jj]])){
            if( !all( c(model[[jj]][ii,'id'] %in% model[['channel']]$id,
                        model[[jj]][ii,'fraction'] >=0,
                        model[[jj]][ii,'fraction'] <=1) ) ){
                stop(paste(jj,model[[jj]][ii,'name'],"is incorrectly specified"))

            }
        }
    }


    
    ## checks on redistribution matrices
    for(ii in c("Fsz","Fsf")){
        if(!(attr(class(model[[ii]]),"package")=="Matrix")){#is(model[[ii]],"Matrix")){
            stop( ii," should be a numeric Matrix" )
        }
        if( !all(dim(model[[ii]])==length(ids)) ){
            stop(paste(ii,"has incorrect dimensions"))
        }
        if( any( model[[ii]] < 0)){
            stop(ii," should have values greater or equal to 0")
        }
        if( any( model[[ii]] > 1)){
            stop(ii," should have values less or equal to 1")
        }
        if( any( colSums(model[[ii]]) > 1+delta)){
            stop(ii," should have column sums that are less or equal to 1")
        }
    }

    ## all output series must have unique names
    req_names <- NULL
    for(jj in names(df_prop)){
        tmp <- df_prop[[jj]]$name[df_prop[[jj]]$role=="output_label"]
        tmp <- unique(unlist(model[[jj]][,tmp]))
        req_names <- unique(c(req_names,tmp))
    }
    if( length(req_names) != length(unique(req_names)) ){
        stop("All output series should have a unique name")
    }

    ## check all next_id are valid
    idx <- !is.na(model$channel$next_id)
    if( !all( model$channel$next_id[idx] %in%  model$channel$id ) ){
        stop("Channels routing to channels not in network, set next_id to NA to represent an outflow")
    }
    ## check for loops - all reaches at top of system must drain to outfall with all reaches being passed through
    to_outlet <- NULL
    for(ii in setdiff(model$channel$id,model$channel$next_id)){    
        id <- ii
        tmp_rec <- ii
                
        cnt <- 0
        while(cnt <= nrow(model$channel)){
            ## work out next id
            id <- model$channel$next_id[model$channel$id==id]
            tmp_rec <- c(tmp_rec,id)
            cnt <- cnt + 1
            ## break if get to outlet
            if(is.na(id)){break}
            if(id %in% to_outlet){break}
        }
        to_outlet <- unique(c(to_outlet,tmp_rec))
    }
    idx <- model$channel$id %in% to_outlet
    if( !all(idx) ){
        stop(paste("Loop in channel network. Check reaches:",paste(model$channel$id[!idx],collapse=TRUE)))
    }
        
    ## verbose printing of head and tail channels
    if(verbose){
        ## print out head channels
        message(paste("The head channels are:",
                      paste(setdiff(model$channel$id,model$channel$next_id),
                            collapse=", "),
                      sep="\n"))
        ## print out tail channels
        message(paste("The channels with outfalls:",
                      paste(model$channel$id[is.na(model$channel$next_id)],
                            collapse=", "),
                      sep="\n"))
    }

    ## work out required input series
    req_names <- NULL
    for(jj in names(df_prop)){
        tmp <- df_prop[[jj]]$name[df_prop[[jj]]$role=="data_series"]
        tmp <- unique(unlist(model[[jj]][,tmp]))
        req_names <- unique(c(req_names,tmp))
    }
    return(req_names)
}


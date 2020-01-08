#' Function to check the model to be used in a dynamic TOPMODEL run
#'
#' @description This function makes some basic consistency checks on a list representing a dynamic TOPMODEL model.
#'
#' @param model a dynamic TOPMODEL list object
#' @param verbose if set prints out further information
#' 
#' @return NULL if tests passed, otherwise fails
#'
#' @details The checks performed and basic 'sanity' checks. they do not check for the logic of the parameter values nor the consistncy of states and parameters.
#' @export
check_model <- function(model, verbose=FALSE){

    ## check all components of the model exist
    components <- c("hillslope","channel","param","Dex","Dsz","gauge","point_inflow")
    idx <- components %in% names(model)
    if( !all(idx) ){
        stop(paste("Missing componets:",paste(components[!idx],collapse=",")))
    }

    ## check components that should be data.frames
    
    ## create a list of data.frames describing the data frames in the model
    df_prop <- list(
        hillslope = data.frame(name = c("id","area","s_bar",
                                        "precip_input","pet_input",
                                        "qex_max","srz_max","srz_0","ln_t0","m","td","tex"),
                               type=c(rep("numeric",3),rep("character",9)),
                               role = c(rep("property",3),rep("data_series",2),rep("parameter",6)),
                               stringsAsFactors=FALSE),
        channel = data.frame(name = c("id","area","length","next_id","precip_input","pet_input","v_ch"),
                             type=c(rep("numeric",4),rep("character",3)),
                             role = c(rep("property",4),rep("data_series",2),"parameter"),
                             stringsAsFactors=FALSE),
        point_inflow = data.frame(
            name = c("name","id","fraction"),
            type=c("character",rep("numeric",2)),
            role = c("data_series",rep("property",2)),
            stringsAsFactors=FALSE),
        gauge <- data.frame(
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
            stop( paste(c("Table",ii,"is missing columns:",
                          hru_properties[[ii]]$names[!idx]),collapse=",") )
        }else{ ## check data types
            tmp <- sapply(model[[ii]],class) # types of the columns
            idx <- tmp[ df_prop[[ii]]$names ] != df_prop[[ii]]$type
            if( any( idx ) ){
                stop( paste(c("Incorrect types in table",ii,"columns:",
                              df_prop[[ii]]$names[idx]),collapse=",") )
            }
        }
    }

    ## ids should be numeric, sequential and start from 1
    ids <- c(model[['hillslope']]$id,model[['channel']]$id)
    rng_ids <- range(ids)
    if( length(idx) != length(unique(idx)) ){ stop("id should be unique") }
    if( !all(is.finite(ids)) ){ stop("ids should be finite") }
    if( !all(range(ids)==c(1,length(ids))) ){ stop("id's should be numbered consecuativly from 1") }

    ## all points should be on a channel with fractions between 0 & 1
    for(jj in c("gauge","point_inflow")){
        for(ii in 1:nrow(model[[jj]])){
            if( !all( c(model[[jj]][ii,'id'] %in% model[['channel']]$id,
                        model[[jj]][ii,'fraction'] >=0,
                        model[[jj]][ii,'fraction'] =<1) ) ){
                stop(paste(jj,model[[jj]][ii,'name'],"is incorrectly specified"))

            }
        }
    }
    
    ## checks on only redistribution matrices
    for(ii in c("Dsz","Dex")){
        if(!(attr(class(model[[ii]]),"package")=="Matrix")){#is(model[[ii]],"Matrix")){
            stop( ii," should be a numeric Matrix" )
        }
        if( !all(dim(model[[ii]])==length(idx)) ){
            stop(paste(ii,"has incorrect dimensions"))
        }
        if( any( model[[ii]] < 0)){
            stop(ii," should have values greater or equal to 0")
        }
        
    }

    ## checks only on param
    if(!is.vector(model$param) | !is.numeric(model$param)){
        stop("param should be a numeric vector")
    }

    ## checks on channel
    if( check_channel ){
        ## check all next_id are valid
        idx <- !is.na(model$channel$next_id)
        if( !all( model$channel$next_id[idx] %in%  model$channel$id ) ){
            stop("Channels routing to channels not in network, set next_id to NA to represent an outflow")
        }
        ## check for loops - all reaches at top of system must drain to outfall
        for(ii in setdiff(model$channel$id,model$channel$next_id)){
            id <- ii
            cnt <- 0
            while(cnt <= nrow(model$channel)){
                ## work out next id
                id <- model$channel$next_id[model$channel$id==id]
                cnt <- cnt + 1
                ## break if get to outlet
                if(is.na(id)){break}
                ## check loop
                if( id==ii ){
                    stop(paste("Loop of river flow incorporating channel id",ii))
                }
            }
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
    }
    
    ## END of checks on individual objects 

    
    ## check unique ids
    for(ii in c("hillslope","channel")){
        if( length(unique(model[[ii]]$id)) < nrow(model[[ii]]) ){
            stop(paste(ii,"HRU ids should be unique"))
        }
    }

    ## check parameters
    pnames <- NULL
    for(ii in c("hillslope","channel")){
        tmp <- hru_properties[[ii]][ hru_properties[[ii]]$is_param , 'name']
        pnames <- c(pnames,
                    unlist( model[[ii]][,tmp] ))
    }
    pnames <- unique(pnames)
    ## check the parameters that are set
    cnames <- names(model$param)
    if(!check_channel & ("v_ch" %in% names(model$channel))){
        ## remove the names of the velocity parameters if they exists
        cnames <- setdiff(cnames,model$channel$v_ch)
    }
    
    if( !all( pnames %in% cnames ) ){
        stop("Not all parameters are specified")
    }
    
    if( !all( cnames %in% pnames ) ){
        stop("Some parameters are not used")
    }

    ## check Wex and Wsat
    hillslope_hru <- as.character( model$hillslope[,'id'] )
    for(ii in c("Wsz","Dex")){
        if( any( colnames(model[[ii]]) != rownames(model[[ii]]) ) ){
            stop(ii," should have identically ordered column and row names")
        }
        if( any( colnames(model[[ii]]) != hillslope_hru ) |
            any( rownames(model[[ii]]) != hillslope_hru ) ){
            stop(ii," column and row names should match the order of the hillslopes in the HRU table")
        }
    }

    ## check Fex and Fsat
    channel_hru <- as.character( model$channel[,'id'] )
    for(ii in c("Fsat","Fex")){
        if( any(colnames(model[[ii]]) != hillslope_hru) ){
            stop(ii," column names should match order of hillslopes in the HRU table")
        }
        if( any( rownames(model[[ii]]) != channel_hru ) ){
            stop(ii," row names should match the order of the channels in the HRU table")
        }
    }
    
    ## check redistribution sums for saturated zone
    tmp <- colSums( rbind(model$Wsz,model$Fsz) )
    if( any(tmp>1) ){
        stop("Saturated flow redistribution fractions sum to greater then 1 for HRUs: ",
             paste(colnames(model$Wsat)[tmp>1],collapse=', '))
    }

    ## check redistribution sums for surface excess
    tmp <- colSums( rbind(model$Wex,model$Fex) )
    if( any(tmp>1) ){
        stop("Surface excess redistribution fractions sum to greater then 1 for HRUs: ",
             paste(colnames(model$Wex)[tmp>1],collapse=', '))
    }

    ## check point tables
    if( check_channel ){
        for(ii in c("gauge","point_inflow")){
            if( !all( model[[ii]]$channel_id %in%  model$channel$id ) ){
                stop(paste("Not all",ii,"are on model channels"))
            }
            if( !all( (model[[ii]]$fraction <= 1) &
                      (model[[ii]]$fraction >= 0 ))){
                stop(paste("Not all",ii,"reach fractions are in [0,1]"))
            }
            if( !( length(unique(model[[ii]]$name)) == nrow(model[[ii]])) ){
                stop(paste(ii,"should have unique names"))
            }
        }
    }
    
    return(NULL)
}


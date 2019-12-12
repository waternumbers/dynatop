#' Function to check the model to be used in a dynamic TOPMODEL run
#'
#' @description This function makes some basic consistency checks on a list representing a dynamic TOPMODEL model.
#'
#' @param model a dynamic TOPMODEL list object
#' @param check_channel flag indicating is additional checks for channel routing are to be carried out
#' @param verbose if set prints out further information
#' 
#' @return NULL if tests passed, otherwise fails
#'
#' @details The checks performed and basic 'sanity' checks. they do not check for the logic of the parameter values nor the consistncy of states and parameters.
#' @export
check_hrus <- function(hru_list){
    if(!is.list(hru_list)){
        stop("HRU list should be a list")
    }

    nms <- names(hru_list)
    if( length(nms)!=length(hru_list) || length(nms)!=length(unique(nms)) ){
        stop("hru_list does not have unique names")
    }
    
    ## function to return properties of each hru type
    f_hru_prop <- function(type){
        if(type=="hillslope"){
            out <- list(properties = c("mapid","area","atb_bar","s_bar"),
                        inputs = c("precip","pet"),
                        parameters = c("srz_max","srz_0","ln_t0","m","td","tex")
                        )
        }
        if(type=="channel"){
            out <- list(properties = c("mapid","area","length"),
                        inputs = c("precip","pet"),
                        parameters = c("v_ch")
                        )
        }
        return(out)
    }

    f_check_prop <- function(x,nms){
        if(!(x$type %in% c("hillslope","channel"))){
            stop("Invlaid type of hru")
        }
        rtype <- switch(x$type,
                        hillslope=c("wex","wsz"),
                        channel=c("wch"))
        p <- f_hru_prop(x$type)
        ## properties
        if(!setequal(p$properties,names(x$properties)) || !is.numeric(x$properties)){
            stop("Incorrect properties are specified")
        }
        ## inputs
        if(!setequal(p$inputs,names(x$inputs)) || !is.character(x$inputs)){
            stop("Incorrect properties are specified")
        }
        ## parameters
        if(!setequal(p$parameters,names(x$parameters)) || !is.character(x$parameters)){
            stop("Incorrect properties are specified")
        }
        ## check redistribution
        for(ii in rtype){
            if(!setequal(nms,names(x[[ii]])) || !is.numeric(x[[ii]]) || sum(x[[ii]])>1 || any(x[[ii]]<0)){
                stop(paste("Error in specification of",ii))
            }
        }
        x$type
    }
            
    ## First perform checks that rely only on individual components of the model
    type <- sapply(hru_list,f_check_prop,nms=nms)

    ## check channel network for loops and outlets (biurfurcastions??)- TO DO
    ## for(ii in setdiff(model$channel$id,model$channel$next_id)){
    ##     id <- ii
    ##     cnt <- 0
    ##     while(cnt <= nrow(model$channel)){
    ##         ## work out next id
    ##         id <- model$channel$next_id[model$channel$id==id]
    ##         cnt <- cnt + 1
    ##         ## break if get to outlet
    ##         if(is.na(id)){break}
    ##         ## check loop
    ##         if( id==ii ){
    ##             stop(paste("Loop of river flow incorporating channel id",ii))
    ##         }
    ##     }
    
    
        
    ##     ## verbose printing of head and tail channels
    ##     if(verbose){
    ##         ## print out head channels
    ##         message(paste("The head channels are:",
    ##                       paste(setdiff(model$channel$id,model$channel$next_id),
    ##                             collapse=", "),
    ##                       sep="\n"))
    ##         ## print out tail channels
    ##         message(paste("The channels with outfalls:",
    ##                       paste(model$channel$id[is.na(model$channel$next_id)],
    ##                             collapse=", "),
    ##                       sep="\n"))
    ##     }
    ## }

    ## TO DO recognise and add points
    ## ## check point tables
    ## if( check_channel ){
    ##     for(ii in c("gauge","point_inflow")){
    ##         if( !all( model[[ii]]$channel_id %in%  model$channel$id ) ){
    ##             stop(paste("Not all",ii,"are on model channels"))
    ##         }
    ##         if( !all( (model[[ii]]$fraction <= 1) &
    ##                   (model[[ii]]$fraction >= 0 ))){
    ##             stop(paste("Not all",ii,"reach fractions are in [0,1]"))
    ##         }
    ##         if( !( length(unique(model[[ii]]$name)) == nrow(model[[ii]])) ){
    ##             stop(paste(ii,"should have unique names"))
    ##         }
    ##     }
    ## }
    
    return(TRUE)
}


#' Function to check the model to be used in a dynamic TOPMODEL run
#'
#' @description This function makes some basic consistency checks on a list representing a dynamic TOPMODEL model.
#'
#' @param model a dynamic TOPMODEL list object
#' @return NULL if tests passed, otherwise fails
#'
#' @details The checks performed and basic 'sanity' checks. they do not check for the logic of the parameter values nor the consistncy of states and parameters.
#' @export
check_model <- function(model){

    ## First perform checks that rely only on individual components of the model
    
    ## checks only on the HRU table
    if(!is.data.frame(model$hru)){ # check HRU is a data frame
        stop("HRU table should be a data.frame")
    }else{
        ## columns of parameters and the type of the contents 
        hru_parameter_columns <- c(srz_max="character",
                                   srz_0="character",
                                   ln_t0="character",
                                   m="character",
                                   td="character",
                                   tex="character")
        ## all required columns in the HRU table and type of the contents
        hru_req_columns <- c(id="numeric",
                             type="character",
                             area="numeric",
                             atb_bar='numeric',
                             precip_input="character",
                             pet_input="character",
                             hru_parameter_columns)
        
        if(!all(names(hru_req_columns) %in% names(model$hru))){ # check it has required columns
            stop( paste(c("HRU table missing columns:",setdiff(names(hru_req_columns), names(model$hru))),collapse=" ") )
        }else{ ## check data types
            tmp <- sapply(model$hru,class) # types of the columns
            if( any( tmp[names(hru_req_columns)] != hru_req_columns ) ){
                stop( paste(c("Incorrect types in HRU table columns:",names(hru_req_columns)[which(tmp[names(hru_req_columns)] != hru_req_columns)]),collapse=" ") )
            }
        }
    }
            
    ## checks on only redistribution matrices
    for(ii in c("Wsat","Fsat","Wex","Fex")){
        if(!is.matrix(model[[ii]]) | !is.numeric(model[[ii]])){
            stop( ii," should be a numeric matrix" )
        }
        if( is.null(colnames(model[[ii]])) | is.null(rownames(model[[ii]])) ){
            stop(ii," should be have column and row names")
        }
        if( any( model[[ii]] < 0)){
            stop(ii," should have values greater or equal to 0")
        }
    }

    ## checks only on param
    if(!is.vector(model$param) | !is.numeric(model$param)){
        stop("param should be a numeric matrix")
    }

    ## END of checks on individual objects return if failed

    ## check unique ids
    if( length(unique(model$hru[,'id'])) < nrow(model$hru) ){
        stop("HRU ids should be unique")
    }

    ## check HRUs are recognised types
    if( !all(model$hru[,'type'] %in% c('hillslope','channel')) ){
        stop("HRU types should be hillslope or channel (case sensitive)")
    }

    ## check parameters
    tmp <- unique( unlist(model$hru[,names(hru_parameter_columns)]) ) ## all unique parameter names
    if( !all( tmp %in% names(model$param) ) ){
        stop("Not all parameters are specified")
    }
    if( !all( names(model$param) %in% tmp ) ){
        stop("Some parameters are not used")
    }

    ## check Wex and Wsat
    hillslope_hru <- as.character( model$hru[ model$hru[,'id']=='hillslope' ,'id'] )
    for(ii in c("Wsat","Wex")){
        if( any( colnames(model[[ii]]) != rownames(model[[ii]]) ) ){
            stop(ii," should have identically ordered column and row names")
        }
        if( any( colnames(model[[ii]]) != hillslope_hru ) |
            any( rownames(model[[ii]]) != hillslope_hru ) ){
            stop(ii," column and row names should match the order of the hillslopes in the HRU table")
        }
    }

    ## check Fex and Fsat
    channel_hru <- as.character( model$hru[ model$hru[,'id']=='channel' ,'id'] )
    for(ii in c("Wsat","Wex")){
        if( any(colnames(model[[ii]]) != hillslope_hru) ){
            stop(ii," column names should match order of hillslopes in the HRU table")
        }
        if( any( rownames(model[[ii]]) != channel_hru ) ){
            stop(ii," row names should match the order of the channels in the HRU table")
        }
    }
    
    ## check redistribution sums for saturated zone
    tmp <- colSums( rbind(model$Wsat,model$Fsat) )
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

    return(NULL)
}


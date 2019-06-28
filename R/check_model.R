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
    
    ## checks only on the HRU tables
    if(!is.data.frame(model$hillslope) | !is.data.frame(model$channel)){ # check HRU is a data frame
        stop("HRU tables should be a data frames")
    }else{
        ## make data frames of columns expected in each table, the type of value expected and if they are a parameter
        hru_properties <- list(
            hillslope = data.frame(name = c("id","area","atb_bar","s_bar",
                                            "precip_input","pet_input",
                                            "srz_max","srz_0","ln_t0","m","td","tex"),
                                   type=c(rep("numeric",4),rep("character",8)),
                                   is_param = c(rep(FALSE,6),rep(TRUE,6)),
                                   stringsAsFactors=FALSE),
            channel = data.frame(name = c("id","area","precip_input","pet_input"),
                                 type=c(rep("numeric",2),rep("character",2)),
                                 is_param = rep(FALSE,4),
                                 stringsAsFactors=FALSE)
        )

        ## check just the HRU tables
        for(ii in names(hru_properties)){
            idx <- hru_properties[[ii]]$name %in% names(model[[ii]])
            browser()
            if( !all( idx ) ){# check it has required columns
                stop( paste(c("HRU table",ii,"missing columns:",
                              hru_properties[[ii]]$names[!idx]),collapse=" ") )
            }else{ ## check data types
                tmp <- sapply(model[[ii]],class) # types of the columns
                idx <- tmp[ hru_properties[[ii]]$names ] != hru_properties[[ii]]$type
                if( any( idx ) ){
                    stop( paste(c("Incorrect types in HRU table",ii,"columns:",
                                  hru_properties[[ii]]$names[idx]),collapse=" ") )
                }
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
        stop("param should be a numeric vector")
    }

    ## END of checks on individual objects return if failed

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
    browser()
    if( !all( pnames %in% names(model$param) ) ){
        stop("Not all parameters are specified")
    }
    
    if( !all( names(model$param) %in% pnames ) ){
        stop("Some parameters are not used")
    }

    ## check Wex and Wsat
    hillslope_hru <- as.character( model$hillslope[,'id'] )
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


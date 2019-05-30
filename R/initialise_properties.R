#' Initialises the hru properties
#'
#' @description The properties of the hru are the non timevarying characteristics that describe the hru. Also contains properties of the combined set of hru
#'
#' @param hru TODO
#' @param param TODO
#' @param type TODO
#'
#' @return a list of properties, each property a named vector
#' @export
initialise_properties <- function(model,param,type=c("hillslope","channel")){

    hru <- model$hru
    
    ## check the input
    type <- match.arg(type)

    if(type=="hillslope"){
        hru_var <- c("id","area","atb_bar","precip_input","pet_input")
        param_var <- c("srz_max","srz_0","ln_t0","m","td","Tex")
    }
    if(type=="channel"){
        hru_var <- c("id","area","precip_input","pet_input")
        param_var <- NULL
    }
    out <- list()
    idx <- which(hru[,'type']==type)
    for(ii in hru_var){
        out[[ii]] <- hru[idx,ii]
        names(out[[ii]]) <- NULL
    }
    for(ii in param_var){
        out[[ii]] <- param[ hru[idx,ii] ]
        names(out[[ii]]) <- NULL
    }

    return(out)
}



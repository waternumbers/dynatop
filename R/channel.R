#' Functions for a channel HRU
#'
#' @description Functions for creating, checking, initialising and evolving a hillslope HRU.
#'
#' @param h a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#'
#' @rdname channel_hru

#' @name channel_hru
#' @export
create_channel <- function(){
    prop_names <- c("area","precip_series","pet_series")
    input_names <- c("precip","lex","lsz")
    output_names <- c("qch")
    param_names <- character(0)
    state_names <- character(0)
    
    out <- list(id=numeric(0),type="channel")
    for(ii in prop_names){
        out$prop[[ii]] <- numeric(0)
    }
    for(ii in input_names){
        out$input[[ii]] <- numeric(0)
    }
    for(ii in output_names){
        out$output[[ii]] <- list(id=numeric(0),weight=numeric(0),val=numeric(0))
    }
    for(ii in param_names){
        out$param[[ii]] <- numeric(0)
    }
    for(ii in state_names){
        out$state[[ii]] <- numeric(0)
    }

    return(out)
}

#' @name channel_hru
#' @export
check_channel <- function(h){
    #warning("hillslope check not implimented")
    return(h)
}

#' @name channel_hru
#' @export
initialise_channel <- function(h){
    return(h)
}

#' @name channel_hru
#' @export
evolve_channel <- function(h,delta_t){

    h$output$qch <- (h$input$lsz + h$input$lex + delta_t*h$input$precip)

    return(h)

}

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
create_channel <- function(id){

    flst <- function(x){
        list(id=x,
             type="channel",
             external_input==c(precip=0,pet=0),
             conectivity=list(lex=list(id=numeric(0),
                                       w=numeric(0)),
                              lsz=list(id=numeric(0),
                                       w=numeric(0)))
             internal_input=c(lex=0,lsz=0),
             output=list("qch"=list(id=numeric(0),weight=numeric(0),val=numeric(0))),
             state=list()
             )
    }
    lapply(id,flst)
}

             
##                          list(
##     lapply(
##     prop_names <- c("area","band")
##     input_names <- c("precip","lex","lsz")
##     output_names <- c("qch")
##     param_names <- c("precip_series","pet_series")
##     state_names <- character(0)
    
##     out <- list(id=numeric(0),type="channel")
##     for(ii in prop_names){
##         out$prop[[ii]] <- numeric(0)
##     }
##     for(ii in input_names){
##         out$input[[ii]] <- numeric(0)
##     }
##     for(ii in output_names){
##         out$output[[ii]] <- list(id=numeric(0),weight=numeric(0),val=numeric(0))
##     }
##     for(ii in param_names){
##         out$param[[ii]] <- numeric(0)
##     }
##     for(ii in state_names){
##         out$state[[ii]] <- numeric(0)
##     }

##     return(out)
## }

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

    ## TO DO the id should be set in the initialisation
    #browser()
    h$output$qch$id <- h$id
    h$output$qch$val <- (h$input$lsz + h$input$lex + delta_t*h$input$precip)

    return(h)

}

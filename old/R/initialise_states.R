#' Initialises the hru states for a dynatop simulation
#'
#' @description Initialise the states of the different hru types. these can be time varying such as store of fluxes or fixed properties such as parmeter values and areas
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#'
#' @name initialise_dynatop
#' @export
initialise_dynatop <- function(model,initial_recharge){

    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }

    hru <- create_hru(model)

    
    for(ii in 1:length(hru)){
#        print(ii)
        
        hru[[ii]] <- switch(hru[[ii]]$type,
                            hillslope = initialise_hillslope(hru[[ii]],initial_recharge),
                            channel = initialise_channel(hru[[ii]])
                            )
    }
    
    return(hru)
}


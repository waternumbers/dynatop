#' Initialises the hru states for a dynatop simulation
#'
#' @description Initialise the states of the different hru types. these can be time varying such as store of fluxes or fixed properties such as parmeter values and areas
#'
#' @param hru a list of one or more dynamic TOPMODEL HRUs
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return An altered hru list with states added
#'
#' @name initialise_dynatop
#' @export
initialise_cpp <- function(hru,initial_recharge){

    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }

    check_hru(hru)

    output <- template_output()

    for(ii in 1:length(hru)){
        ## only the hillsope has states
        if( hru[[ii]]$type == "hillslope" ){
            ## set so zero surface excess
            hru[[ii]]$state['ex'] <- 0

            ## initialise root zone based on fraction constrained within 0,1
            hru[[ii]]$state['rz'] <- max( min( hru[[ii]]$param['srz_0'],1) ,0 ) *
                hru[[ii]]$param['srz_max']

            ## initialise based on steady state of saturated zone given constant recharge from unsturated zone

            ## maximum lateral flow from saturated zone
            hru[[ii]]$state['lsz_max'] <- exp(hru[[ii]]$param['ln_t0'])*
                hru[[ii]]$prop['s_bar'] # for unit width
            hru[[ii]]$state['lsz_max'] <- hru[[ii]]$state['lsz_max'] / hru[[ii]]$prop['delta_x'] ## for specified width

            ## initialise
            hru[[ii]]$state['lsz'] <- hru[[ii]]$state['lsz_in'] <- q0
            hru[[ii]]$state['lsz'] <- min(hru[[ii]]$state['lsz'],
                                          hru[[ii]]$state['lsz_max'])

            ## compute the deficit
            gamma <- hru[[ii]]$prop['atb_bar'] - hru[[ii]]$param['ln_t0']
            hru[[ii]]$state['sz'] <- hru[[ii]]$param['m']*
                (gamma + log(hru[[ii]]$state['lsz']))
            hru[[ii]]$state['sz'] <- max(0,hru[[ii]]$state['sz'])

            ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
            hru[[ii]]$state['uz'] <- max(0,q0 * hru[[ii]]$param['td'] * hru[[ii]]$state['sz'], na.rm=TRUE)
        }
        hru[[ii]]$output <- output
    }

    return(hru)
}


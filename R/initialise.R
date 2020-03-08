#' Initialisation of the hillslope HSU store
#'
#' @description Initialises a Dynamic TOPMODEL int eh simpliest way possible.
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/s
#'
#' @details Initialises a Dynamic TOPMODEL int eh simpliest way possible.
#'
#' @return \code{model} with additional variable representing the initial state values added to the model object.
#'
#' @export
initialise <- function(model,initial_recharge){

    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }

    ## ###########################################
    ## Initialise the states
    #browser()
    ## maximum lateral flow from saturated zone
    model$hillslope$l_szmax <- exp( model$param[model$hillslope$ln_t0] )*model$hillslope$s_bar

    ## initialise the root zone
    model$hillslope$s_rz <- pmax( pmin( model$param[ model$hillslope$s_rz0 ], 1) ,0 ) * model$param[ model$hillslope$s_rzmax ]

    model$hillslope$s_uz  <- model$hillslope$s_sf <- 0

    ## initialise the saturated zone
    model$hillslope$sum_l_sz_in <- model$hillslope$l_sz <- pmin(model$hillslope$l_szmax,initial_recharge)
    model$hillslope$s_sz <- pmax(0, model$param[ model$hillslope$m ]*( log(model$hillslope$l_szmax) - log(model$hillslope$l_sz)))

    model$channel$sum_l_sz_in <- initial_recharge
    ## initialise the unsaturated zone based on recharge
    ##model$hillslope$state$s_uz <- pmax(0, initial_recharge * model$hillslope$param$t_d * model$hillslope$state$s_sz)

    ##saturated zone
    ## initialise the saturated zone
    ##model$hillslope$state$l_sz_in <- model$hillslope$state$l_sz <- pmin(model$hillslope$state$l_szmax,initial_recharge)

    ## compute the deficit
    ##gamma <- sum(model$hillslope$attr$area*(model$hillslope$attr$atb_bar - model$hillslope$param$ln_t0))  / sum(model$hillslope$attr$area)
    ##model$hillslope$state$s_sz <- pmax(0, model$hillslope$param$m*(gamma + log(model$hillslope$state$l_sz)))

    ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983

    ## model$hillslope$state$s_uz <- pmax(0, initial_recharge * model$hillslope$param$t_d * model$hillslope$state$s_sz)


    return(model)
}

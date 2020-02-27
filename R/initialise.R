#' Contains functions for initialising and checking the hillslope HSU store
#'
#' @description TODO
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @details TODO
#' 
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#'
#' @name initialise
NULL

#' @rdname initialise
## returns a list of varaibles names by list elements expected in a hillslope HSU store
required_vrbl <- function(type=c("hillslope","channel")){
    type <- match.arg(type)
    
    if(type=="hillslope"){
        out <- list(attr = c("id","atb_bar","s_bar","area","delta_x","sz_band","sf_band","precip","pet"), # values associated with catchment HSU
                    param = c("q_sfmax","s_rzmax","s_rz0","ln_t0","m","t_d","t_sf"), # parameter names - also attributes but seperate since handled differently
                    input = c("p","e_p","l_sf"), # tempory store of input, not needed for next timestep
                    state = c("s_sf","s_rz","s_uz","s_sz","l_sz_in","l_sz","l_szmax"), # states - required for next timestep
                    flux = c("q_ex_rz","q_rz_ex","q_rz_uz","q_uz_sz","q_uz_ex","q_sz_ex","e_t","l_sf","l_sz_in_t","Q_minus_t","Q_plus_t") # fluxes within and from the HSU not needed for the next time step
                    )
    }
    if(type=="channel"){
        out <- list(attr = c("id","area","sz_band","sf_band","precip","pet"),
                    input = c("p","e_p","l_sf"),
                    param = character(0),
                    state = c("l_sz_in"),
                    flux = c("s_ch","l_sz_in_t"))
    }
    return(out)
}

#' @rdname initialise
#' @export
initialise <- function(model,initial_recharge, return_sim=FALSE){
        
    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }
    
    ## convert the model to the correct list type
    mdl <- store_to_sim(model)

    ## ###########################################
    ## Initialise the states
    #browser()
    ## maximum lateral flow from saturated zone
    mdl$hillslope$state$l_szmax <- exp(mdl$hillslope$param$ln_t0)*mdl$hillslope$attr$s_bar
    
    ## initialise the root zone
    mdl$hillslope$state$s_rz <- pmax( pmin( mdl$hillslope$param$s_rz0, 1) ,0 ) * mdl$hillslope$param$s_rzmax

    mdl$hillslope$state$s_uz[] <- 0
    
    ## initialise the saturated zone
    mdl$hillslope$state$l_sz_in <- mdl$hillslope$state$l_sz <- pmin(mdl$hillslope$state$l_szmax,initial_recharge)
    mdl$hillslope$state$s_sz <- pmax(0, mdl$hillslope$param$m*(log(mdl$hillslope$state$l_szmax) - log(mdl$hillslope$state$l_sz)))
    
    ## initialise the unsaturated zone based on recharge
    ##mdl$hillslope$state$s_uz <- pmax(0, initial_recharge * mdl$hillslope$param$t_d * mdl$hillslope$state$s_sz)
    
    ##saturated zone
    ## initialise the saturated zone
    ##mdl$hillslope$state$l_sz_in <- mdl$hillslope$state$l_sz <- pmin(mdl$hillslope$state$l_szmax,initial_recharge)
    
    ## compute the deficit
    ##gamma <- sum(mdl$hillslope$attr$area*(mdl$hillslope$attr$atb_bar - mdl$hillslope$param$ln_t0))  / sum(mdl$hillslope$attr$area)
    ##mdl$hillslope$state$s_sz <- pmax(0, mdl$hillslope$param$m*(gamma + log(mdl$hillslope$state$l_sz)))
    
    ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
   
    ## mdl$hillslope$state$s_uz <- pmax(0, initial_recharge * mdl$hillslope$param$t_d * mdl$hillslope$state$s_sz)
    

    ## convert form back if required
    if( !return_sim ){
        mdl <- sim_to_store(mdl)
    }
    
    return(mdl)
}

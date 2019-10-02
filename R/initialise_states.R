#' Initialises the hru states
#'
#' @description Initialise the states of the different hru types. these can be time varying such as store of fluxes or fixed properties such as parmeter values and areas
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of vectors of states and properties for use in computation
#' 
#' @name initialise_hru
#' @rdname initialise_hru
NULL

#' Hillslope initialisation
#' @export
#' @rdname initialise_hru
initialise_hillslope <- function(model,initial_recharge){

    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }
    
    ## create a list and fill with fixed properties
    hru_var <- c("id","area","atb_bar","s_bar","precip_input","pet_input")
    param_var <- c("srz_max","srz_0","ln_t0","m","td","tex")

    hillslope <- list()
    for(ii in hru_var){
        hillslope[[ii]] <- model$hillslope[,ii]
        names(hillslope[[ii]]) <- NULL
    }
    for(ii in param_var){
        hillslope[[ii]] <- model$param[ model$hillslope[,ii] ]
        names(hillslope[[ii]]) <- NULL
    }

    ## Now compute values based on initialisation
    
    ## compute the max flow from saturated area
    hillslope$lsz_max <- exp(hillslope$ln_t0)*hillslope$s_bar

    ## initialise based on steady state of saturated zone given constant recharge from unsturated zone
    q_0 <- initial_recharge # q_0 is specific area recharge in m/hr

    hillslope$quz <- rep(q_0,length(hillslope$id))

    ## a least squares solution with total outflow equals total inflow to saturated zone
    ## Assume all outflow to river
    X_init <- rbind( diag(nrow(model$Wsat)) - model$Wsat , colSums(model$Fsat) ) %*% diag(hillslope$area)
    y_init <- c( hillslope$area*q_0, sum(hillslope$area)*q_0 )
    
    hillslope$lsz <- q_0 # assume constant across catchment
    try({
        hillslope$lsz <- as.numeric( solve( t(X_init)%*%X_init, t(X_init)%*%y_init ) )
        #        names(hillslope$lsz) <- NULL
        if(any( hillslope$lsz <= 0 )){
            stop("Solution for saturated zone initialisation produced negative or zero flows.\n Using constant recharge across catchment")
        }
    })
        
    
    ## K_init <- diag(1/hillslope$area) %*% (diag(nrow(model$Wsat)) - model$Wsat) %*%
    ##     diag(hillslope$area)
    ## hillslope$lsz <- q_0 # assume constant across catchment
    ## try({
    ##     hillslope$lsz <- solve( K_init , hillslope$quz )
    ##     names(hillslope$lsz) <- NULL
    ##     if(any( hillslope$lsz <= 0 )){
    ##         stop("Solution for saturated zone initialisation produced negative or zero flows.\n Using constant recharge across catchment")
    ##     }
    ## })
    hillslope$lsz <- pmin( hillslope$lsz , hillslope$lsz_max )
    
    ## compute the saturated zone storage deficit based on the flow
    hillslope$ssz <- -hillslope$m * log(hillslope$lsz/exp(hillslope$ln_t0-hillslope$atb_bar))
    hillslope$ssz <- pmax(  hillslope$ssz, 0 )
    
    ## assume full recharge from saturated zone
    #hillslope$quz <- rep(q_0,length(hillslope$id))

    ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
    hillslope$suz <- hillslope$quz * hillslope$td * hillslope$ssz

    ## initialise root zone based on fraction constrained within 0,1
    hillslope$srz <- pmax( pmin( hillslope$srz_0, 1) ,0 ) * hillslope$srz_max

    ## initialise surface storage to 0
    hillslope$ex <- rep(0,length(hillslope$id))

    return(hillslope)
}


#' Channel initialisation
#' @export
#' @rdname initialise_hru
initialise_channel <- function(model){

    hru_var <- c("id","area","precip_input","pet_input")

    channel <- list()
    for(ii in hru_var){
        channel[[ii]] <- model$channel[,ii]
        names(channel[[ii]]) <- NULL
    }

    ## and set store
    channel$store <- rep(0,length(channel$id))
    
    return(channel)
}

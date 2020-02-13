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
#' @rdname hillslope
#' @export
#'

#' @name hillslope
#' returns a list of varaibles names by list elements expected in a hillslope HSU store
required_vrbl <- function(type=c("hillslope","channel")){
    type <- match.arg(type)
    
    if(type=="hillslope"){
        out <- list(attr = c("id","atb_bar","s_bar","area","delta_x","sz_band","sf_band"),
                    series = c("precip","pet"),
                    input = c("precip","pet","lsf","lsz"),
                    output = c("et","lsf","lsz"),
                    param = c("qsf_max","srz_max","srz_0","ln_t0","m","td","tsf"),
                    state = c("sf","rz","uz","sz","lsz_max"), # stages - not dependent upon time step
                    flux = c("qex_rz","qrz_ex","qrz_uz","quz_sz","quz_ex","qsz_ex") # fluxes depend on timestep
                    )
    }
    if(type=="channel"){
        out <- list(attr = c("id","area","sz_band","sf_band"),
                    series = c("precip","pet"),
                    input = c("precip","pet","lsf","lsz"),
                    output = c("lch"),
                    param = character(0),
                    state = c("sch"),
                    flux = character(0))
    }
    return(out)
}

#' @name hillslope
#' create an empty store
empty_store <- function(model,type=c("hillslope","channel")){
    #browser()
    type <- match.arg(type)

    ## get a list of all variables that should exist
    vrbl <- required_vrbl(type)

    tmp <- model[[type]]
    
    ## initialise depending upon type
    out <- list(attr=list(),
                series=list(),
                input=list(),
                param=list(),
                state=list(),
                flux=list())
    
    for(ii in vrbl$attr){
        out$attr[[ii]] <- unname( tmp[,ii] )
    }
    for(ii in vrbl$series){
        out$series[[ii]] <- unname( tmp[,ii] )
    }
    for(ii in vrbl$input){
        out$input[[ii]] <-  rep(0,nrow(tmp))
    }
    for(ii in vrbl$output){
        out$output[[ii]] <-  rep(0,nrow(tmp))
    }
    for(ii in vrbl$param){
        out$param[[ii]] <- unname( model$param[ tmp[,ii] ] )
    }
    for(ii in vrbl$state){
        out$state[[ii]] <- rep(0,nrow(tmp)) ## initialise all stores and fluxes as 0
    }
    for(ii in vrbl$flux){
        out$flux[[ii]] <- rep(0,nrow(tmp)) ## initialise all stores and fluxes as 0
    }
    
    return(out)
}


#' @name hillslope
#' @export
initialise <- function(model,initial_recharge){

    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }

    ## #########################################
    ## Initialise the output
    mdl <- list()

    
    ## collapsed matrices as lists
    ## ith element of list contains the parents (and fractions) for the
    ## hsu with id = i
    #browser()
    tmp <- as.matrix(summary(model$Fsf))
    tmp <- lapply(by(tmp[,2:3],tmp[,1],unique,simplify=FALSE),as.list)
    mdl$sf <- rep(list(NULL),max(as.numeric(names(tmp))))
    mdl$sf[as.numeric(names(tmp))] <- tmp
    
    tmp <- as.matrix(summary(model$Fsz))
    tmp <- lapply(by(tmp[,2:3],tmp[,1],unique,simplify=FALSE),as.list)
    mdl$sz <- rep(list(NULL),max(as.numeric(names(tmp))))
    mdl$sz[as.numeric(names(tmp))] <- tmp
    
    ## add the channel hsu
    mdl$channel <- empty_store(model,"channel")

    ## add the hillslope hsu
    hs <- empty_store(model,"hillslope")

    ## maximum lateral flow from saturated zone
    hs$state$lsz_max <- exp(hs$param$ln_t0)*hs$attr$s_bar

    ## initialise the root zone
    hs$state$rz <- pmax( pmin( hs$param$srz_0, 1) ,0 ) * hs$param$srz_max

    ## initialise the saturated zone
    hs$output$lsz <- hs$input$lsz <- pmin(hs$state$lsz_max,initial_recharge)

    ## compute the deficit
    #browser()
    gamma <- sum(hs$attr$area*(hs$attr$atb_bar - hs$param$ln_t0))  / sum(hs$attr$area)
    hs$state$sz <- pmax(0, hs$param$m*(gamma + log(hs$output$lsz)))

    ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
    hs$state$uz <- pmax(0, initial_recharge * hs$param$td * hs$state$sz)

    mdl$hillslope <- hs

    return(mdl)
}


##     ## #########################################
##     ## Hill slope initialisation
    
##     hs <- empty_store(model,"hillslope")
    
##     ## initialise root zone based on fraction constrained within 0,1
##     hs$rz <- pmax( pmin( hs$srz_0, 1) ,0 ) * hs$srz_max

##     ## maximum lateral flow from saturated zone
##     hs$lsz_max <- exp(hs$ln_t0)*hs$s_bar

    
##     ## initialise based on steady state of saturated zone given constant recharge from unsturated zone
##     q_0 <- initial_recharge

##     ## a least squares solution with total outflow equals total inflow
##     ## to saturated zone
##     tryCatch({
##         hs$lsz <- as.numeric( solve( Diagonal(length(hs$id)) - model$Dsz[hs$id,hs$id],
##                                     rep(q_0,length(hs$id))) )
##     },error = function(e){
##         warning("Solution for saturated zone initialisation produced negative or zero flows. Using constant recharge across catchment")
##         hs$lsz <- rep(q_0,length(hs$id))
##     })

##     ## threshold at max
##     hs$lsz <- pmin( hs$lsz , hs$lsz_max )

##     ## compute the deficit
##     gamma <- sum(hs$area*(hs$atb_bar - hs$ln_t0))  / sum(hs$area)
##     hs$sz <- hs$m*(gamma + log(hs$lsz))
##     hs$sz <- pmax(0,hs$sz)

##     ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
##     hs$uz <- pmax(0, rep(q_0,length(hs$id)) * hs$td * hs$sz, na.rm=TRUE)

##     cmdl$hillslope <- hs
    
##     ## ########################################
##     ## initialise the channel
##     cmdl$channel <- empty_store(model,"channel")
    
##     return(cmdl)
## }

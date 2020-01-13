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

    ## function for initialising the hillslope
    fhillslope <- function(model,initial_recharge){

        ## sort by order
        model$hillslope <- model$hillslope[order(model$hillslope$band),]

        hru_var <- c("id","area","s_bar","precip_input","pet_input","delta_x","band")

        lst <- list()
        for(jj in nrow(model$hillslope)){
            h <- list(type="hillslope")
            for(ii in hru_var){
                h$properties[[ii]] <- unname( model$hillslope[jj,ii] )
            }
            for(ii in param_var){
                h$parameters[[ii]] <- unname( model$param[ model$hillslope[jj,ii] ] )
            }

            ## initialise surface storage to 0
            h$states$ex <- 0

            ## initialise root zone based on fraction constrained within 0,1
            h$states$rz <- max( min( h$parameters$srz_0, 1) ,0 ) *
                h$parameters$srz_max

            ## initialise based on steady state of saturated zone given constant recharge from unsturated zone

            ## maximum lateral flow from saturated zone
            h$states$lsz_max <- exp(h$parameters$ln_t0)*h$properties$s_bar # for unit width
            h$states$lsz_max <- h$states$lsz_max / h$properties$delta_x # for specified width

            ## initialise
            h$states$lsz <- h$states$lsz_in <- initial_recharge
            h$states$lsz <- min(h$states$lsz,h$states$lsz_max)

            ## compute the deficit
            gamma <- sum(h$properties$area*(h$properties$atb_bar - h$parameters$ln_t0))  / sum(h$properties$area)
            h$states$sz <- h$parameters$m*(gamma + log(h$states$lsz))
            h$states$sz <- max(0,h$sz)

            ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
            h$states$uz <- pmax(0, intial_recharge * h$parameters$td * h$parameters$sz, na.rm=TRUE)

            ## take only what we need
            lst[[jj]] <- h
        }
        return(lst)
    }

    ## function for initialising the channel
    fchannel <- function(model){
        lst <- list()
        for(jj in 1:nrow(model$channel)){
            h <- list()
            hru_var <- c("id","area","precip_input","pet_input")

            for(ii in hru_var){
                h$properties$[[ii]] <- unname( model$channel[jj,ii] )
            }
            lst[[jj]] <- h
        }
        return(lst)
    }

    hru <- c(fhillslope(model,initial_recharge),
             fchannel(model))

    ## add redistribution
    tbl_ex <- as.matrix(summary(model$Dex))
    tbl_sz <- as.matrix(summary(model$Dsz))

    area <- sapply(hru,FUN=function(x){x$properties$area})
    for(ii in 1:length(hru)){
        ## do ex redistribution
        idx <- tbl_ex[,2] == hru[[ii]]$id
        hru[[ii]]$properties$rex$idx <- tbl_ex[idx,1]
        hru[[ii]]$properties$rex$w <- tbl_ex[idx,3]*area[ii]/area[idx]
        ## do ex redistribution
        idx <- tbl_sz[,2] == hru[[ii]]$id
        hru[[ii]]$properties$rsz$idx <- tbl_sz[idx,1]
        hru[[ii]]$properties$rsz$w <- tbl_sz[idx,3]*area[ii]/area[idx]
    }

    return(hru)
}


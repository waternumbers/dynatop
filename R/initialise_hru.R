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
#' @export
blank_hsu <- function(type=c("hillslope","channel")){
    type <- match.arg(type)
    
    if(type=="hillslope"){
        out <- list(type=type,
                    prop = setNames(numeric(5),c("id","atb_bar","s_bar",'area','delta_x')),
                    input = setNames(numeric(4),c('precip','pet','lsf','lsz')),
                    output = setNames(numeric(3),c('et','lsf','lsz')),
                    series = setNames(character(2),c('precip','pet')),
                    param = setNames(numeric(7),
                                     c("qsf_max","srz_max","srz_0","ln_t0","m","td","tex")),
                    state = setNames(numeric(5),
                                      c("ssf","srz","suz","ssz","lsz_max")), # stores and constants
                    flux = setNames(numeric(6),
                                    c("qex_rz","qrz_ex","qrz_uz","quz_sz","quz_ex","qsz_ex")) # vertical fluxes
                    )
    }
    if(type=="channel"){
        out <- list(type=type,
                    prop = setNames(numeric(2),c("id","area")),
                    input = setNames(numeric(4),c('precip','pet','lsf','lsz')),
                    output = setNames(numeric(3),c('et','lsz','lsf')),
                    series = setNames(character(2),c('precip','pet')),
                    states = c("sch"=0)
                    )
    }
    return(out)
}

#' @name hillslope
#' @export
initialise_hsu <- function(model,initial_recharge){
    
    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }

    ## get common things used in the ordering
    id <- c(model$hillslope$id,model$channel$id)
    sz_band <- c(model$hillslope$sz_band,model$channel$sz_band)
    sf_band <- c(model$hillslope$sf_band,model$channel$sf_band)
    
    sq <- list(sf = id[order(sf_band)],
                sz = id[order(sz_band)])
    ## #########################################
    ## Initialise the output
    hsu <- rep(list(NULL),max(id))

    print("initialising hillslopes")
    ## go through the hillslope table
    for(ii in 1:nrow(model$hillslope)){
        
        jj <- model$hillslope$id[ii]
        hsu[[jj]] <- blank_hsu("hillslope")
        tmp <- names(hsu[[jj]]$prop)
        hsu[[jj]]$prop <- setNames(as.numeric(model$hillslope[ii,tmp]),tmp)
        hsu[[jj]]$series <- setNames(as.character(model$hillslope[ii,c('precip_series','pet_series')]),
                                     c("precip","pet"))
        
        tmp <- names(hsu[[jj]]$param)
        hsu[[jj]]$param <- setNames(as.numeric(model$param[ as.character(model$hillslope[ii,tmp])]),tmp)
        ## initialise root zone based on fraction constrained within 0,1
        hsu[[jj]]$state['srz'] <- max( min( hsu[[jj]]$param['srz_0'], 1) ,0 ) * hsu[[jj]]$param['srz_max']
        ## set max saturated zone flux
        hsu[[jj]]$state['lsz_max'] <- exp(hsu[[jj]]$param['ln_t0'])*hsu[[jj]]$prop['s_bar']
        hsu[[jj]]$output['lsz'] <- min(initial_recharge,hsu[[jj]]$state['lsz_max'])
        hsu[[jj]]$input['lsz'] <- hsu[[jj]]$output['lsz'] ## TODO this is a  bodge
    }
    ## initialise the storage deficits
    gamma <- c(0,0)
    for(ii in model$hillslope$id){
        ## compute the deficit
        gamma[1] <- gamma[1] + hsu[[ii]]$prop['area']*(hsu[[ii]]$prop['atb_bar'] - hsu[[ii]]$param['ln_t0'])
        gamma[2] <- gamma[2] + hsu[[ii]]$prop['area']
    }
    gamma[1] <- gamma[1]/gamma[2]
    for(ii in model$hillslope$id){
        hsu[[ii]]$state['ssz'] <- hsu[[ii]]$param['m']*(gamma[1] + log(hsu[[ii]]$output['lsz']))
        hsu[[ii]]$state['ssz'] <- max(0,hsu[[ii]]$state['ssz'])
    }
    
    print("initialising channels")
    ## go through the channel table
    
    for(ii in 1:nrow(model$channel)){
        jj <- model$channel$id[ii]
        hsu[[jj]] <- blank_hsu("channel")
        tmp <- names(hsu[[jj]]$prop)
        hsu[[jj]]$prop <- setNames(as.numeric(model$channel[ii,tmp]),tmp)
        hsu[[jj]]$series <- setNames(as.character(model$channel[ii,c('precip_series','pet_series')]),
                                     c("precip","pet"))
    }


    print("initialising parents")
    #browser()
    ## collapsed matrices as lists
    tmp <- as.matrix(summary(model$Fsf))
    sf_prnt <- by(tmp[,2:3],tmp[,1],unique,simplify=FALSE)
    for(ii in names(sf_prnt)){
        jj <- as.numeric(ii)
        hsu[[jj]]$parent$sf <- sf_prnt[[ii]]
    }
    
    tmp <- as.matrix(summary(model$Fsz))
    sz_prnt <- by(tmp[,2:3],tmp[,1],unique,simplify=FALSE)
    for(ii in names(sz_prnt)){
        jj <- as.numeric(ii)
        hsu[[jj]]$parent$sz <- sz_prnt[[ii]]
    }
    
    return(list(hsu=hsu,
                sq=sq))
}

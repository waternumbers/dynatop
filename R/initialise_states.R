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

    
    ## extract redistribution tables
    redist <- list()
    redist$lex <- as.matrix(summary(model$Dex))
    redist$lsz <- as.matrix(summary(model$Dsz))

    ## sort by order - this is the basic sequencing at the moment
    model$hillslope <- model$hillslope[order(model$hillslope$band),]
    
    ## list for output
    hru <- list()

    ## do hillslopes first
    cnt <- 1
    ## convert parameters to values
    h0 <- create_hillslope()
    hillslope <- model$hillslope
    for(ii in names(h0$param)){
        hillslope[,ii] <- unname( model$param[hillslope[,ii]] )
    }
    hillslope <- split(hillslope, seq(nrow(hillslope)),drop=TRUE)
    for(hs in hillslope){
        h <- h0 #create_hillslope()
        h$id <- hs$id
        prop <- intersect(names(h$prop),names(hs))
        h$prop <- as.list(hs[prop])
        param <- intersect(names(h$param),names(hs))
        h$param <- as.list(hs[param])
        hru[[cnt]] <- initialise_hillslope(h,initial_recharge)
        cnt <- cnt+1
    }
    
    ## do channels second
    h <- create_channel()
    channel <- model$channel
    for(ii in names(h$param)){
        channel[,ii] <- unname( model$param[channel[,ii]] )
    }
    channel <- split(channel, seq(nrow(channel)),drop=TRUE)
    for(hs in channel){
        h <- create_channel()
        h$id <- hs$id
        prop <- intersect(names(h$prop),names(hs))
        h$prop <- as.list(hs[prop])
        param <- intersect(names(h$param),names(hs))
        h$param <- as.list(hs[param])
        hru[[cnt]] <- initialise_channel(h)
        cnt <- cnt+1
    }
    
    ## work out the reweighting
    area <- sapply(hru,FUN=function(x){x$prop$area})
    for(jj in names(redist)){
        idx <- split(redist$lex[,1],redist$lex[,2])
        w <- split(redist$lex[,3],redist$lex[,2])
        kdx <- sapply(split(redist$lex[,2],redist$lex[,2]),unique)
        for(ii in 1:length(kdx)){
            k <- kdx[ii]
            hru[[k]]$output[[jj]]$id <- idx[[ii]]
            hru[[k]]$output[[jj]]$w <- w[[ii]]*area[k]/area[idx[[ii]]]
        }
    }
    
    return(hru)
}


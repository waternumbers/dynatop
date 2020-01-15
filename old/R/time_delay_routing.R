#' Function to apply variable time delay routing
#'
#' @description This function makes a model along with the specification of a point and inputs to the channel. It applies a time delay historgram routing based on the channel length and velocities.
#'
#' @param model a dynamic TOPMODEL list object
#' @param channel_inflows an xts object of flows to the channels
#' @param diffuse_inflows an xts object of diffuse inflows
#' @param point_inflows and xts object of point inflows
#' @param initial_conditions a list of initial conditions
#'
#' @return An xts object of flows at the guage locations
#'
#' @details Timings are taken from the channel_inflows, other inputs are ignored unless there timings match
#' @export
time_delay_routing <- function(model,channel_inflows,
                               point_inflows=NULL,
                               diffuse_inflows=NULL,
                               initial_conditions=list()){

    
    ## check model including channel
    check_model(model,check_channel=TRUE,verbose=FALSE)

    ## start default initial conditions
    ic <- list()

    ## check channel_inflows - object type and match model
    ts <- check_obs(channel_inflows,as.character(model$channel[['id']]))
    ic$channel_inflows <- setNames(as.numeric(channel_inflows[1,]),
                                   names(channel_inflows))
    
    ## check diffuse inflows - time step, variables
    use_diffuse <- is.xts(diffuse_inflows)
    if( use_diffuse ){
        tmp <- setdiff(names(diffuse_inflows),paste(model$channel[['id']]))
        if(length(tmp)>0){
            warning(paste("The following diffuse inflows are not used:",
                          paste(tmp,collapse=", "),
                          sep="\n"))
        }
        if( !all(index(diffuse_inflows) == index(channel_inflows)) ){
            stop("Diffuse inflow data time index does not match that of channel inflows")
        }
        ic$diffuse_inflows <- setNames(as.numeric(diffuse_inflows[1,]),
                                       names(diffuse_inflows))
    }

    ## check point inflow locations and time series
    if( nrow(model$point_inflow)>0 ){
        if( !is.xts(point_inflows) ){
            stop("Point inflow series is required as an xts object")
        }
        tmp <- setdiff(model$point_inflow$name,names(point_inflows))
        if(length(tmp)>0){
            stop(paste("The following point inflows are not specified:",
                       paste(tmp,collapse=", "),
                       sep="\n"))
        }
        if( !all( index(point_inflows) == index(channel_inflows)) ){
            stop("Diffuse inflow data time index does not match that of channel inflows")
        }
        ic$point_inflows <- setNames(as.numeric(point_inflows[1,]),
                                     names(point_inflows))  
    }

    ## add alternative initial conditions
    for(ii in intersect(names(ic),names(initial_conditions))){
        jj <- intersect(names(ic[[ii]]),names(initial_conditions[[ii]]))
        ic[[ii]][jj] <- initial_conditions[[ii]][jj]
    }
    
    ## add diffuse inflows to channel_inflows to speed up computations
    if( use_diffuse ){
        for(ii in intersect(names(channel_inflows), names(diffuse_inflows))){
            channel_inflows[,ii] <- channel_inflows[,ii] + diffuse_inflows[,ii]
            ic$channel_inflows[ii] <- ic$channel_inflows[ii] + ic$diffuse_inflows[ii]
        }
    }
    
    ## initialise the output
    out <- reclass( matrix(NA,nrow(channel_inflows),
                           nrow(model$gauge)), match.to=channel_inflows)
    names(out) <- model$gauge$name

    ## Compute the time delays
    time_of_travel <- compute_time_delay(model)

    ## function to make polynonial
    fpoly <- function(h,f){
        if(f>h){stop("Interval for polynomial compution if wrong")}            
        dmax <- floor(h)+1
        dmin <- floor(f)
        ply <- rep(0,dmax)
        idx <- (dmin+1):dmax
        ply[idx] <- 1
        if(h!=f){
            ## take away potentially partial value around foot
            ply[dmin+1] <- ply[dmin+1] - (f - floor(f))
            ## take away potentially partial value around head
            ply[dmax] <- ply[dmax] - (1 - (h - floor(h)))
        }
        return(ply/sum(ply))
    }
    
    ## Loop gauges
    for(ii in model$gauge$name){
        ## initialise the point - set to 0
        out[,ii] <- 0

        
        ## loop channel+diffuse upstream
        upstream_channels <- colnames(time_of_travel$head_to_gauge)
        upstream_channels <- upstream_channels[is.finite(time_of_travel$head_to_gauge[ii,])]
        for(jj in upstream_channels){
            ## time in timesteps to foot and head of channel length adjusted for
            ## in reach gauge
            h2g <- time_of_travel$head_to_gauge[ii,jj] / ts$step
            rt <- time_of_travel$reach_time[jj] / ts$step
            f2g <- pmax(0,h2g-rt)
            qfrac <- (h2g-f2g)/rt
            ## generate polynomial
            ply <- fpoly(h2g,f2g)*qfrac
            ## evaluate filter
            npad <- length(ply)-1
            x <- c(rep(ic$channel_inflows[jj],npad),as.numeric(channel_inflows[,jj]))
            #browser()
            tmp <- filter(x,ply,method="conv",sides=1)
            x <- as.numeric( filter(x,ply,method="conv",sides=1) )
            if(npad>0){x <- x[-(1:npad)]}
            out[,ii] <- out[,ii] + x
            
        }

        ## apply point inputs
        upstream_points <- colnames(time_of_travel$point_to_gauge)
        upstream_points <- upstream_points[is.finite(time_of_travel$point_to_gauge[ii,])]
        ## loop points
        for(jj in upstream_points){
            p2g <- time_of_travel$point_to_gauge[ii,jj]/ts$step
            f2g <- floor(p2g)
            h2g <- floor(p2g)
            ## generate polynomial
            ply <- fpoly(h2g,f2g)
            ## evaluate filter - pad for initial conditions           
            npad <- length(ply)-1
            x <- c(rep(ic$point_inflows[jj],npad),as.numeric(point_inflows[,jj]))
            x <- as.numeric( filter(x,ply,method="conv",sides=1) )
            if(npad>0){x <- x[-(1:npad)]}
            out[,ii] <- out[,ii] + x
            
        }
    }
    
    #out <- reclass( out, match.to=channel_inflows)
    return(out)
}

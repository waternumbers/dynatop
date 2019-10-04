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

    ## check channel_inflows - object type and match model
    ts <- check_obs(channel_inflows,model$channel[['id']])

    ## check diffuse inflows - time step and variables
    use_diffuse <- is.xts(diffuse_inflows)
    if( use_diffuse ){
        tmp <- setdiff(names(diffuse_inflows),paste(model$channel[['id']]))
        if(length(tmp)>0){
            warning(paste("The following diffuse inflows are not used:",
                          paste(tmp,collapse=", "),
                          sep="\n"))
        }
        if( index(diffuse_inflows) != index(channel_inflows) ){
            stop("Diffuse inflow data time index does not match that of channel inflows")
        }
    }

    ## check point inflow locations and time series
    if( nrow(model$point_inflow)>0 ){
        if( !is.xts(point_inflows) ){
            stop("Point inflow series is not an xts object")
        }
        tmp <- setdiff(model$point_inflow$name,names(point_inflows))
        if(length(tmp)>0){
            stop(paste("The following point inflows are not specified:",
                       paste(tmp,collapse=", "),
                       sep="\n"))
        }
    }

    ## check initial conditions - match model and point inflows - fix if not
    ic <- list(channel_inflows = setNames(rep(0,ncol(channel_inflows)),
                                          names(channel_inflows)),
               diffuse_inflows = setNames(rep(0,ncol(diffuse_inflows)),
                                          names(diffuse_inflows)),
               point_inflows = setNames(rep(0,ncol(point_inflows)),
                                        names(point_inflows))
               )
    if(is.list(initial_conditions)){
        for(ii in intersect(names(ic), names(initial_conditions))){
            idx <- intersect(names(ic[[ii]]),names(initial_conditions[[ii]]))
            ic[[ii]][idx] <- initial_conditions[[ii]][idx]
        }
    }

    ## add diffuse inflows to channel_inflows
    for(ii in intersect(names(channel_inflows), names(diffuse_inflows))){
        channel_inflows[,ii] <- channel_inflows[,ii] + diffuse_inflows[,ii]
        ic$channel_inflows[ii] <- ic$channel_inflows[ii] + ic$diffuse_inflows[ii]
    }

    ## initialise the output
    out <- reclass( matrix(NA,nrow(channel_inflows),
                           nrow(model$gauges)), match.to=channel_inflows)
    names(out) <- model$gauges$name

    ## Compute the time delays
    time_of_travel <- compute_time_delay(model)

    ## Loop gauges
    for(ii in model$gauge$name){
        ## initialise the point - set to 0
        out[,ii] <- 0

        ## loop channel+diffuse upstream
        upstream_channels <- colnames(time_of_travel$reach_to_gauge)
        upstream_channels <- upstream_channels[is.finite(time_of_travel$reach_to_gauge[ii,])]
        for(jj in upstream_channels){
            ## time in timesteps to foot and head of channel length
            tw <- (time_of_travel$reach_to_gauge[ii,jj] +
                   c(0,time_of_travel$reach_time[jj]))/(3600*ts$time_step)
            ## time window can be negative if gauge is within reach
            ## fix this and work out fraction of flow with >0 time
            twz <- pmax(tw,0)
            qfrac <- diff(twz)/time_of_travel$reach_time[jj]
            ## create polynomial
            ## remember first element is zero time offset in filter
            ply <- rep(0,floor(twz[2])+1)
            ## fill with ones
            idx <- (floor(twz[1])+1):(floor(twz[2])+1)
            ply[idx] <- 1
            ## take away potentially partial value around foot
            idx <- floor(twz[1])+1
            ply[idx] <- ply[idx] - (1 - (twz[1] - floor(twz[1])))
            ## take away potentially partial value around head
            idx <- floor(twz[1])+1
            ply[idx] <- ply[idx] - (1 - (twz[2] - floor(twz[2])))
            ## scale includeing  flow fraction
            ply <- ply*qfrac/sum(ply)
            ## run filter
            out[,ii] <- out[,ii] + filter(channel_inflows[,jj],ply,
                                          method="conv",sides=1,
                                          init=ic$channel_inflows[jj])
        }

        ## apply point inputs
        upstream_points <- colnames(time_of_travel$point_to_gauge)
        upstream_points <- upstream_points[is.finite(time_of_travel$point_to_gauge[ii,])]
        ## loop points
        for(jj in upstream_points){
            ## time in timesteps to foot and head of channel length
            tw <- (time_of_travel$point_to_gauge[ii,jj] +
                   c(0,ts$time_step))/(3600*ts$time_step)
            ## time window can be negative if gauge is within reach
            ## fix this and work out fraction of flow with >0 time
            twz <- pmax(tw,0)
            qfrac <- diff(twz)/diff(tw)
            ## create polynomial
            ## remember first element is zero time offset in filter
            ply <- rep(0,floor(twz[2])+1)
            ## fill with ones
            idx <- (floor(twz[1])+1):(floor(twz[2])+1)
            ply[idx] <- 1
            ## take away potentially partial value around foot
            idx <- floor(twz[1])+1
            ply[idx] <- ply[idx] - (1 - (twz[1] - floor(twz[1])))
            ## take away potentially partial value around head
            idx <- floor(twz[1])+1
            ply[idx] <- ply[idx] - (1 - (twz[2] - floor(twz[2])))
            ## scale includeing  flow fraction
            ply <- ply*qfrac/sum(ply)
            ## run filter
            out[,ii] <- out[,ii] + filter(channel_inflows[,jj],ply,
                                          method="conv",sides=1,
                                          init=ic$channel_inflows[jj])
        }
    }

    return(out)
}

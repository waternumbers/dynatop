#' Function to compute the time delays between gauges and points in the network
#'
#' @description This function computes the time taken for between point in the network TODO
#'
#' @param model a dynamic TOPMODEL list object
#' 
#' @return A list with two matrices one of time delays between the foot of reaches and the gauges, one with time differences between point inflow and gauges
#'
#' @details This is done by a simple search - not very efficent for large networ#' @export
compute_time_delay <- function(model){
    
    ## check the model including the channel
    check_model(model,check_channel=TRUE,verbose=FALSE)
    
    ## compute the time to travel down each reach
    reach_time <- setNames(model$channel$length/model$param[model$channel$v_ch],
                           model$channel$id)
    
    ## get the next reach as a string
    next_reach <- setNames(paste(model$channel$next_id),
                           model$channel$id)
    
    ## compute times for foot of one rech to foot of another
    foot_to_foot <- matrix(NA,nrow(model$channel),nrow(model$channel),
                           dimnames=list(model$channel$id,
                                         model$channel$id))
    for(ii in row.names(foot_to_foot)){
        id <- ii
        ts <- 0
        while(id %in% names(reach_time)){
            foot_to_foot[ii,id] <- ts
            ts <- ts + reach_time[id]
            id <- next_reach[id]
        }
    }
    
    ## work out time to gauges
    reach_to_gauge <- matrix(NA,nrow(model$gauge),nrow(model$channel),
                             dimnames=list(model$gauge$name,
                                           model$channel$id))
    for(ii in 1:nrow(model$gauge)){
        gr <- paste(model$gauge$channel_id[ii])
        ## time to foot of the gauged reach
        reach_to_gauge[ii,] <- foot_to_foot[colnames(reach_to_gauge),gr]
        ## correct for time above foot
        reach_to_gauge[ii,] <- reach_to_gauge[ii,] -
            (1-model$gauge$fraction[ii])*reach_time[gr]
    }
    ## times between point inflows and gauges
    point_to_gauge <- matrix(NA,nrow(model$gauges),nrow(model$point_inflow))
    pr <- paste(model$point_inflow$channel_id[ii])
    for(ii in 1:nrow(model$gauge)){
        gr <- paste(model$gauge$channel_id[ii])
        reach_to_gauge[ii,] <- foot_to_foot[colnames(reach_to_gauge),gr]
        ## correct for time above foot of gauge
        reach_to_gauge[ii,] <- reach_to_gauge[ii,] -
            (1-model$gauge$fraction[ii])*reach_time[gr]
        ## correct for time above foot of point
        reach_to_gauge[ii,] <- reach_to_gauge[ii,] +
            (1-model$point_inflow$fraction)*reach_time[pr]
    }
    
    return(list(foot_to_foot=foot_to_foot,
                point_to_gauge=point_to_gauge,
                reach_to_gauge=reach_to_gauge,
                reach_time=reach_time))
}

        

#' Function to compute the time delays between gauges and points in the network
#'
#' @description This function computes the time taken for between point in the network TODO
#'
#' @param model a dynamic TOPMODEL list object
#' 
#' @return A list with two matrices one of time delays between the foot of reaches and the gauges, one with time differences between point inflow and gauges
#'
#' @details This is done by a simple search - not very efficent for large networks. The time returned is in h
#' @export
compute_time_delay <- function(model){
    
    ## check the model including the channel
    check_model(model,check_channel=TRUE,verbose=FALSE)

    
    ## compute the time to travel down each reach
    reach_time <- setNames(model$channel$length/model$param[model$channel$v_ch],
                           model$channel$id)
    
    ## get the next reach as a string
    next_reach <- setNames(paste(model$channel$next_id),
                           model$channel$id)

    ## compute times from the head of one reach to the head of another
    head_to_head <- matrix(NA,nrow(model$channel),nrow(model$channel),
                           dimnames=list(model$channel$id,
                                         model$channel$id))
    ## one row for one starting channel
    for(ii in row.names(head_to_head)){
        id <- ii
        ts <- 0
        while(id %in% names(reach_time)){
            head_to_head[ii,id] <- ts
            ts <- ts + reach_time[id]
            id <- next_reach[id]
        }
    }
    #browser()
    ## time from head of reach to gauge
    head_to_gauge <- matrix(NA,nrow(model$gauge),nrow(model$channel),
                            dimnames=list(model$gauge$name,
                                          model$channel$id))
    for(ii in 1:nrow(model$gauge)){
        gr <- paste(model$gauge$channel_id[ii])
        head_to_gauge[ii,] <- head_to_head[colnames(head_to_gauge),gr] +
            (model$gauge$fraction[ii]*reach_time[gr])
    }

    
    point_to_gauge <- matrix(NA,nrow(model$gauge),nrow(model$point_inflow),
                             dimnames=list(model$gauge$name,
                                           model$point_inflow$name))
    for(ii in 1:nrow(model$gauge)){
        gr <- paste(model$gauge$channel_id[ii])
        for(jj in 1:nrow(model$point_inflow)){      
            ir <- paste(model$point_inflow$channel_id[jj])
            point_to_gauge[ii,jj] <- head_to_head[ir,gr] +
                (model$gauge$fraction[ii]*reach_time[gr]) -
                (model$point_inflow$fraction[jj]*reach_time[ir])
        }
    }
    point_to_gauge[point_to_gauge<0] <- NA
    
    return(list(point_to_gauge=point_to_gauge/3600,
                head_to_gauge=head_to_gauge/3600,
                reach_time=reach_time/3600))
}


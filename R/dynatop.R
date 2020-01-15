#' Run dynamic topmodel
#' @param model TODO
#' @param obs_data TODO
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#' @param sim_time_step simulation timestep in hours, default value of NULL results in data time step
#' @param use_states should the states in the model be used (default FALSE)
#'
#' @details use_states, currently does not impose any checks
#'
#' @export
dynatop <- function(hru,obs_data,
                    initial_recharge=NULL,
                    sim_time_step=NULL,
                    sz_opt=list(omega=0.7,
                                theta=0.7)){

    ## check the hru
    check_hru(hru)
    input_series <- unique( unlist( sapply( hru, function(x){x$series}) ) )
    
    ## check input and get model timestep
    ts <- check_obs(obs_data,input_series,
                    sim_time_step)

    if(!is.null(initial_recharge)){
        hru <- initialise(hru,initial_recharge)
    }

    ## initialise the output
    channel_names <- unlist( sapply(hru,function(x){ifelse(x$type=="channel",
                                                           x$id,NULL)}) )
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(channel_names)),
                              match.to=obs_data )
    names(channel_inflow) <- channel_names

    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    #browser()
    for(it in 1:nrow(obs_data)){
        print(it)
        
        ## set the external inputs
        obs_it <- obs[it,]
        
        ## set the channel inflow to zero to be added to
        channel_inflow[it,] <- 0
        
        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            rcpp_dynatop(hru,obs_it)

            ## add channel_inflow
            tmp <- sapply(hru[channel_names],function(x){x$output['qch']})
            channel_inflow[it,] <- channel_inflow[it,] + tmp
        } ## end of inner loop
        ## correct channel_inflow to be averge of substeps
        channel_inflow[it,] <- channel_inflow[it,] / ts$n_step
        
    } ## end of timestep loop

    return( list(hru=hru,
                 channel_input = channel_inflow ) )
}

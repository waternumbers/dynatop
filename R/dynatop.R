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
dynatop <- function(model,obs_data,
                    initial_recharge=NULL,
                    sim_time_step=NULL,
                    use_states=FALSE,
                    sz_opt=list(omega=0.7,
                                theta=0.7,
                                max_iter=100,
                                tol=1e-6)){

   
    ## check the model
    input_series <- check_model(model)
    ## check input and get model timestep
    ts <- check_obs(obs_data,input_series,
                    sim_time_step)

    if(use_states){ # then just take states from the model object
        hru <- model$states
    }else{ # initialise the model
        hru <- initialise_dynatop(model,initial_recharge)
    }
      
    ## initialise the redistributed fluxes
    redist_flux <- list(
        precip = rep(0,length(hru)),
        pet = rep(0,length(hru)),
        lex = rep(0,length(hru)),
        lsz = rep(0,length(hru)),
        qch = rep(0,length(hru))
    )

    
    ## compute the summaries that are used in the simulation
    comp_summ <- list(
        precip_names  = sapply(hru,function(x){x$prop$precip_series}),
        pet_names = sapply(hru,function(x){x$prop$pet_series}),
        clear_redist = setdiff(names(redist_flux),c("precip","pet")),
        channel_idx = sapply(hru,function(x){ifelse(x$type=="channel",x$id,NA)})
    )
    comp_summ$channel_idx[is.finite(comp_summ$channel_idx)]
                          
    ## initialise the output
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(comp_summ$channel_idx)), match.to=obs_data )
    names(channel_inflow) <- comp_summ$channel_idx
    
    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    
    for(it in 1:nrow(obs_data)){
        
        ## set the precip and pet inputs
        redist_flux$precip <- obs_data[it,comp_summ$precip_names]
        redist_flux$pet <- obs_data[it,comp_summ$pet_names]

        ## set the channel inflow to zero to be added to
        channel_inflow[it,] <- 0

        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## clear all lateral fluxes
            for(ii in comp_summ$clear_redist){
                redist_flux[[ii]][] <- 0
            }

            for(ii in 1:length(hru)){

                h <- hru[[ii]]
                
                ## set inputs
                for(jj in names(hru[[ii]]$input)){
                    hru[[ii]]$input[[jj]] <- redist_flux[[jj]][h$id]
                }

                ## evolve hru
                hru[[ii]] <- switch(h$type,
                             "hillslope" = evolve_hillslope(hru[[ii]],ts$sub_step,sz_opt),
                             "channel" = evolve_channel(hru[[ii]],ts$sub_step))
                
                ## copy fluxes back
                for(jj in names(h$ouput)){
                    redist_flux[[jj]][hru[[ii]]$out[[jj]]$id] <-
                        redist_flux[[jj]][hru[[ii]]$out[[jj]]$id] + hru[[ii]]$output[[jj]]$val
                }
                
            } ## end of hru loop
            ## copy to model output
            channel_inflow[it,] <- channel_inflow[it,] + redist_flux$qch[comp_summ$channel_idx]
        } ## end of inner loop
        channel_inflow[it,] <- channel_inflow[it,]/(3600*ts$step)
    } ## end of timestep loop

    model$states <- hru

    return( list(model=model,
                 channel_input = channel_inflow ) )
}

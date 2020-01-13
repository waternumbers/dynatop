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

    #browser()
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


    ## initialise the vertical flux stores
    lateral_flux <- list(
        lex = rep(0,length(hru)),
        lsz = rep(0,length(hru))
    )

    ## simple function for solving ODE
    fode <- function(a,b,x0,t){
        b <- pmax(b,1e-10)
        unname( x0*exp(-b*t) + (a/b)*(1-exp(-b*t)) )
    }

    ## initialise the output
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(model$channel$id)), match.to=obs_data )
    names(channel_inflow) <- model$channel$id

    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){

        ## set the inputs to the hillslope
        precip <- obs_data[it,]
        pet <- obs_data[it,]

        ## set the channel inflow to zero to be added to
        channel_inflow[it,] <- 0

        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## clear all lateral fluxes
            for(ii in names(lateral_flux)){
                lateral_flux[[ii]][] <- 0
            }

            for(ii in 1:length(hru)){

                h <- hru[[ii]]

                ## set inputs
                for(jj in names(h$input)){
                    if(jj == precip){
                        h$input$precip$val <- precip[h$input$precip$id]
                        next
                    }
                    if(jj == pet){
                        h$input$pet$val <- pet[h$input$precip$id]
                        next
                    }
                    h$input[[jj]]$val <- lateral_flux[[jj]][h$input[[jj]]$id]
                }


                if(h$type == "hillslope"){
                    ##evolve states
                    h <- evolve_hillslope(h,ts$time_step,sz_opt)
                    ## copy fluxes back
                    for(jj in names(h$ouput)){
                        lateral_flux[[jj]][h$out[[jj]]$id] <-
                            lateral_flux[[jj]][h$out[[jj]]$id] + h$output[[jj]]$val
                }

                if(h$type==h$channel){
                    tmp <- lateral_flux$lex[ii] +
                        lateral_flux$lex[ii] +
                        precip[h$input$precip$id]*ts$time_step
                    channel_inflow[it,ii] <- channel_inflow[it,ii]/(3600*ts$time_step)
                }

                ## put hru back into list
                hru[[ii]] <- h
            } ## end of hru loop
        } ## end of inner loop
    } ## end of timestep loop

    model$states <- hru

    return( list(model=model,
                 channel_input = channel_inflow ) )
}

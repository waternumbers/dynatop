#' Run dynamic topmodel
#' @param model TODO
#' @param obs_data TODO
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#' @param sim_time_step simulation timestep in hours, default value of NULL results in data time step
#' @param use_states should the states in the model be used (default FALSE)
#' @param sz_opt TODO
#' @param mass_check return time series of mass balance errors
#'
#' @details use_states, currently does not impose any checks
#'
#' @export
dynatop <- function(model,obs_data,
                    initial_recharge=NULL,
                    sim_time_step=NULL,
                    use_states=FALSE,
                    mass_check=FALSE,
                    return_states=NULL,
                    sz_opt=list(omega=0.7,
                                theta=0.7)){

    ## check the model
    input_series <- check_model(model,use_states=use_states)

    ## check input and get model timestep
    ts <- check_obs(obs_data,input_series,
                    sim_time_step)

    ## initialise the model
    if( use_states ){
        model <- initialise(model,initial_recharge)
    }

    ## convert to lists for simulation
    ## this creates hillslope, channel, sqnc and lateral_flux
    list2env(convert_form(model),as.environment(-1))


    ## simple function for solving ODE
    fode <- function(a,b,x0,t){
        #b <- pmax(b,1e-10)
        ebt <- exp(-b*t)
        kappa <- (1-ebt)/b
        kappa[b==0] <- t
        ## kappa <- pmin(t,(1-ebt)/b)
        #x <- unname( x0*ebt + (a/b)*(1-ebt) )
        x <- unname( x0*ebt + a*kappa )
        return(x)
    }

    ## initialise the channel inflow output
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(model$channel$attr$id)), match.to=obs_data )
    names(channel_inflow) <- model$channel$attr$id

    ## initialise the mass check output if required
    if( mass_check ){
        mass_errors <- matrix(NA,nrow(obs_data)*ts$n_sub_step,6)
        colnames(mass_errors) <- c("DateTime","step","s_sf","s_rz","s_uz","s_sz")
    }

    ## check and initialise the state outputs
    if( length(return_states)>0 ){
        if( !("POSIXct" %in% class(return_states)) ){
            stop("Times for returning states should be POSIXct object")
        }
        idx <- index(obs_data) %in% return_states
        return_states <- rep(list(NULL),sum(idx))
        names(return_states) <- index(obs_data)[idx]
        return_states[['idx']] <- idx
        return_state$flag <- TRUE
    }else{
        return_state <- list(flag=FALSE)
    }


    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){
        ## set the inputs to the hillslope and channel
        ## set as rate m/s
        hillslope$p <- obs_data[it,hillslope$precip]/ts$step
        hillslope$e_p <- obs_data[it,hillslope$pet]/ts$step
        channel$p <- obs_data[it,channel$precip]/ts$step
        channel$e_p <- obs_data[it,channel$pet]/ts$step

        ## set the state of the channel to 0 since wany to accumulate over the timestep
        channel$s_ch[] <- 0

        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## stor previous version of model if need mass check
            if( mass_check ){
                hs0 <- hillslope
                ch0 <- channel
            }

            ## Step 1: Distribute any surface storage downslope
            lateral_flux$sf[] <- 0 # remove fluxes from previous time step since not needed - should be over written but...

            for(bnd in sqnc$sf_band){ ## loop all bands of surface
                ## hillslope
                idx <- bnd$hillslope ## index of hillslope HSUs
                if(any(idx)){
                    for(ii in idx){
                        hillslope$l_sf[ii] <- sum( hillslope$Fsf[[ii]]$x *
                                                   lateral_flux$sf[ hillslope$Fsz[[ii]]$j ])
                    }

                    hillslope$l_sf[idx] <- hillslope$l_sf[idx] / hillslope$area[idx] # inflow in m depth accrued over sub step

                    ## compute the new state value
                    tilde_sf <- fode( hillslope$l_sf[idx]/ts$sub_step, 1/hillslope$t_sf[idx],
                                     hillslope$s_sf[idx],ts$sub_step)
                    ## work out out flow
                    hillslope$l_sf[idx] <- hillslope$s_sf[idx] + hillslope$l_sf[idx] - tilde_sf
                    hillslope$s_sf[idx] <- tilde_sf
                    lateral_flux$sf[ hillslope$id[idx] ]  <- hillslope$l_sf[idx]*hillslope$area[idx]
                }



                ## channel
                idx <- bnd$channel #which(model$channel$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        channel$l_sf[ii] <- sum( channel$Fsf[[ii]]$x *
                                                 lateral_flux$sf[ channel$Fsf[[ii]]$j ] )
                    }
                    channel$l_sf[idx] <- channel$l_sf[idx]/channel$area[idx]
                    ## compute the new store of channel input in depth
                    channel$s_ch[idx] <- channel$s_ch[idx] + channel$l_sf[idx]
                    lateral_flux$sf[ channel$id[idx] ]  <- 0
                }
            }

            if( mass_check ){ mc_s_ch_sf <- channel$s_ch }
            ## print(paste("s_sf",it,inner,
            ##             sum(model$hillslope$state$s_sf*model$hillslope$attr$area) +
            ##             sum(model$channel$flux$s_ch*model$channel$attr$area) -
            ##             sum(model0$hillslope$state$s_sf*model0$hillslope$attr$area)-
            ##             sum(model0$channel$flux$s_ch*model0$channel$attr$area),
            ##             sum(model0$hillslope$state$s_sf*model0$hillslope$attr$area)+
            ##             sum(model0$channel$flux$s_ch*model0$channel$attr$area),
            ##             sum(model$hillslope$state$s_sf*model$hillslope$attr$area)+
            ##             sum(model$channel$flux$s_ch*model$channel$attr$area)
            ##             ))

            ## Step 2: solve the root zone for hillslope elements
            #browser()
            ## evaluate max integral of flow to rootzone
            hillslope$q_sf_rz <- pmin( hillslope$q_sfmax*ts$sub_step,hillslope$s_sf )
            hillslope$s_sf <- hillslope$s_sf - hillslope$q_sf_rz

            ## solve ODE
            tilde_rz <- fode( hillslope$p + (hillslope$q_sf_rz/ts$sub_step),
                             hillslope$e_p/hillslope$s_rzmax,
                             hillslope$s_rz,ts$sub_step )

            ## work out actual evapotranspiration by mass balance
            hillslope$e_t <- hillslope$s_rz + hillslope$p*ts$sub_step - tilde_rz
            ## new storage value
            hillslope$s_rz <- pmin(tilde_rz,hillslope$s_rzmax)

            ## split root zone flow
            tmp <- tilde_rz - hillslope$s_rz
            saturated_index <- hillslope$s_sz <= 0 # which areas are saturated
            hillslope$q_rz_sf <- tmp * saturated_index
            hillslope$q_rz_uz <- tmp * !saturated_index

            ## print(paste("s_rz",it,inner,
            ##             sum((model0$hillslope$state$s_rz +
            ##                 model$hillslope$input$p*ts$sub_step +
            ##                 model$hillslope$flux$q_sf_rz -
            ##                 tmp - model$hillslope$state$s_rz)*
            ##                 model0$hillslope$attr$area),
            ##             sum(model0$hillslope$state$s_rz*model0$hillslope$attr$area),
            ##             sum((model$hillslope$state$s_rz -
            ##                  model$hillslope$input$p*ts$sub_step -
            ##                  model$hillslope$flux$q_sf_rz +
            ##                  tmp)*model$hillslope$attr$area)
            ##             ))

            ## Step 3: Unsaturated zone
            ## solve ODE
            ##browser()
            tilde_uz <- fode( hillslope$q_rz_uz/ts$sub_step,
                             1 / (hillslope$t_d * hillslope$s_sz),
                             hillslope$s_uz,ts$sub_step )

            hillslope$q_uz_sz <- hillslope$s_uz + hillslope$q_rz_uz - tilde_uz
            hillslope$s_uz <- tilde_uz

            ## print(paste("s_uz",it,inner,
            ##             sum((model0$hillslope$state$s_uz +
            ##                  model$hillslope$flux$q_rz_uz -
            ##                  model$hillslope$flux$q_uz_sz -
            ##                  model$hillslope$state$s_uz)*model0$hillslope$attr$area),
            ##             sum(model0$hillslope$state$s_uz*model0$hillslope$attr$area),
            ##             sum((model$hillslope$state$s_uz -
            ##                  model$hillslope$flux$q_rz_uz +
            ##                  model$hillslope$flux$q_uz_sz)*model$hillslope$attr$area)
            ##             ))

            ## Step 4: Solve saturated zone
            ## move current states to values for start of the time step
            hillslope$sum_l_sz_in_t <- hillslope$sum_l_sz_in # total inflow at start of time step
            hillslope$l_sz_t <- hillslope$l_sz # total outflow at start of time step
            channel$sum_l_sz_in_t <- channel$sum_l_sz_in
            ## evaluate the values of the initial components of the kinematic solution
            hillslope$Q_minus_t <- pmin( hillslope$sum_l_sz_in_t, hillslope$l_szmax ) # inflow to saturated zone at start of timestep
            hillslope$Q_plus_t <- pmin( hillslope$l_sz_t, hillslope$l_szmax ) # outflow at start of timestep

            ## update flows by looping through bands
            for(bnd in sqnc$sz_band){
                ## hillslope
                idx <- bnd$hillslope
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        hillslope$sum_l_sz_in[ii] <- sum( hillslope$Fsz[[ii]]$x *
                                                          lateral_flux$sz[ hillslope$Fsz[[ii]]$j ] )
                    }
                    hillslope$sum_l_sz_in[idx] <- hillslope$sum_l_sz_in[idx]/hillslope$area[idx] ## current inflow in m/s

                    ## compute the new state value
                    hillslope$Q_minus_tDt[idx] <- pmin( hillslope$sum_l_sz_in[idx],hillslope$l_szmax[idx] ) # current inflow to saturated zone

                    qbar <- (hillslope$Q_minus_t[idx] + hillslope$Q_minus_tDt[idx] + hillslope$Q_plus_t[idx])/3
                    cbar <- (qbar*hillslope$delta_x[idx])/hillslope$m[idx]

                    lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/hillslope$delta_x[idx]
                    lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/hillslope$delta_x[idx]

                    k <- lambda_prime * hillslope$Q_plus_t[idx] +
                        (1-lambda_prime) * hillslope$Q_minus_t[idx] +
                        cbar*hillslope$q_uz_sz[idx]/hillslope$delta_x[idx]

                    hillslope$l_sz[idx] <- pmin( (k - (1-lambda)*hillslope$Q_minus_tDt[idx])/lambda , hillslope$l_szmax[idx] )
                    lateral_flux$sz[ hillslope$id[idx] ] <- hillslope$l_sz[idx]*hillslope$area[idx]

                }
                ## channel
                idx <- bnd$channel
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        channel$sum_l_sz_in[ii] <- sum( channel$Fsz[[ii]]$x *
                                                        lateral_flux$sz[ channel$Fsz[[ii]]$j ] )
                    }
                    channel$sum_l_sz_in[idx] <- channel$sum_l_sz_in[idx]/channel$area[idx]
                    laterl_flux$sz[ channel$id[idx] ]  <- 0
                }
            }


            ## update volumes in hillslope
            tilde_sz <- hillslope$s_sz +
                ts$sub_step*(hillslope$l_sz_t + hillslope$l_sz)/2 -
                ts$sub_step*(hillslope$sum_l_sz_in_t + hillslope$sum_l_sz_in)/2 -
                hillslope$q_uz_sz

            hillslope$s_sz <- pmax(0,tilde_sz)
            hillslope$q_sz_sf <- hillslope$s_sz - tilde_sz

            ## update volume of inflow to channel
            channel$s_ch <- channel$s_ch +
                ts$sub_step*(channel$sum_l_sz_in_t + channel$sum_l_sz_in)/2



            ## step 5 - correct the stores for saturation flows
            ##browser()
            ## if( any(hillslope$s_sz <= 1e-5) ){
            ##     browser()
            ## }

            saturated_index <- hillslope$s_sz <= 0
            hillslope$q_uz_sf <- hillslope$s_uz*saturated_index
            hillslope$s_uz <- hillslope$s_uz * !saturated_index
            hillslope$s_sf <- hillslope$s_sf +
                hillslope$q_rz_sf +
                hillslope$q_sz_sf + hillslope$q_uz_sf

            ## mass check for iteration
            if( mass_check ){
                mass_errors[ ((it-1)*ts$n_sub_step) + inner,] <-
                    c(it,inner,
                      sum(hs0$s_sf*hs0$area)+
                      sum(ch0$s_ch*ch0$area) +
                      sum(hillslope$q_rz_sf*hillslope$area) +
                      sum(hillslope$q_sz_sf*hillslope$area) -
                      sum(hillslope$q_sf_rz*hillslope$area) -
                      sum(hillslope$s_sf*hillslope$area) -
                      sum(mc_s_ch_sf*channel$area),
                      sum((hs0$s_rz +
                           hillslope$p*ts$sub_step +
                           hillslope$q_sf_rz -
                           hillslope$q_rz_sf -
                           hillslope$q_rz_uz -
                           hillslope$e_t -
                           hillslope$s_rz)*
                          hs0$area),
                      sum((hs0$s_uz +
                           hillslope$q_rz_uz -
                           hillslope$q_uz_sz -
                           hillslope$q_uz_sf -
                           hillslope$s_uz)*hs0$area),
                      -sum(hs0$s_sz*hs0$area)+
                      sum(mc_s_ch_sf*channel$area) +
                      sum(hillslope$q_uz_sz*hillslope$area) -
                      sum(hillslope$q_sz_sf*hillslope$area) -
                      -sum(hillslope$s_sz*hillslope$area) -
                      sum(channel$s_ch*channel$area)
                      )
                #print( mass_errors[ ((it-1)*ts$n_sub_step) + inner,] )
            }
            ##browser()
            ## step 6 - channel inflow - at the moment a volume / area
            channel_inflow[it,] <- channel$area *
                ( channel$s_ch + channel$p*ts$step ) /
                (ts$step)

            ## handle returning states
            if( return_states$flag ){
                return_states[[index(obs_data)[it]]] <- get_states(hillslope,"hillslope")
            }


        }

    } ## end of timestep loop

    ## copy states back into model
    tmp <- get_states(hillslope,"hillslope")
    model$hillslope <- merge(model$hillslope,get_states(hillslope,"hillslope"),by="id",all=c(FALSE,TRUE))
    tmp <- get_states(channel,"channel")
    model$channel <- merge(model$channel,get_states(channel,"channel"),by="id",all=c(FALSE,TRUE))


    out <- list(model = model,
                channel_input = channel_inflow)

    if(mass_check){
        #browser()
        mass_errors <- as.data.frame(mass_errors)
        mass_errors[,'DateTime'] <- index(obs)[mass_errors[,'DateTime']]
        out[['mass_errors']] <- mass_errors
    }
    if( length(return_states)>0 ){
        return_states <- return_states[setdiff(names(return_states),c("idx","flag"))]
        out[['return_states']] <- return_states
    }
    return( out )
}

#' Run dynamic topmodel
#' @param model A Dynamic TOPMODEL object (see vignette)
#' @param obs_data an xts object containing equally spaced time series of observed data.
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/s
#' @param sim_time_step simulation timestep in hours, default value of NULL results in data time step
#' @param use_states should the states in the model be used (default FALSE)
#' @param return_states a vector of POSIXct objects (e.g. from xts) giving the time stamp at which the states should be returned
#' @param sz_opt a named list e.g. list(omega=1,theta=1) of parameters controlling the kinematic wave solution for the saturated zone.
#' @param mass_check return time series of mass balance errors
#'
#' @details use_states, currently does not impose any checks on the state values. return_states fives the states at the end of the timestep. The default values of sz_opt (omega=1,theta=1) ensure there are no negative fluxes. Other values may produce negative values in whiich case theese are set to 0 and a warning issued.
#'
#' @export
dynatop_sim <- function(model,obs_data,obs_index,ts,
                        mass_check=FALSE,
                        return_states=NULL,
                        sz_opt=list(omega=1,
                                theta=1)){


    ## convert model to variables in function
    list2env(dynatop::convert_form(model),as.environment(-1))
    
    #browser()
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
    channel_inflow <- matrix(NA,nrow(obs_data),length(model$channel$id))
    colnames(channel_inflow) <- model$channel$id

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
        idx <- obs_index %in% return_states
        return_states <- rep(list(NULL),sum(idx))
        names(return_states) <- obs_index[idx]
        return_states[['idx']] <- idx
        return_states$flag <- TRUE
    }else{
        return_states <- list(flag=FALSE)
    }


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")


    for(it in 1:nrow(obs_data)){
        ##print(it)

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

            ## remove fluxes from previous time step since not needed
            ## should be over written but...
            lateral_flux$sf[] <- 0
            lateral_flux$sz[] <- 0
            
            ## Step 1: Distribute any surface storage downslope
            for(idx in sqnc$sf){ ## loop all bands of surface in hillslopes
                ## compute surface flux

                hillslope$l_sf[idx] <- lateral_flux$sf[ hillslope$id[idx] ] / hillslope$area[idx]
                ## compute the new state value
                tilde_sf <- fode( hillslope$l_sf[idx]/ts$sub_step, 1/hillslope$t_sf[idx],
                                 hillslope$s_sf[idx],ts$sub_step)
                ## work out out flow
                hillslope$l_sf[idx] <- hillslope$s_sf[idx] + hillslope$l_sf[idx] - tilde_sf
                hillslope$s_sf[idx] <- tilde_sf

                for(ii in idx){
                    lateral_flux$sf[ hillslope$sf_dir[[ii]]$idx ]  <-
                        lateral_flux$sf[ hillslope$sf_dir[[ii]]$idx ] + 
                        hillslope$sf_dir[[ii]]$frc * hillslope$l_sf[ii] * hillslope$area[ii]
                }
            }

            

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

            ## Step 3: Unsaturated zone
            ## solve ODE
            ##browser()
            tilde_uz <- fode( hillslope$q_rz_uz/ts$sub_step,
                             1 / (hillslope$t_d * hillslope$s_sz),
                             hillslope$s_uz,ts$sub_step )

            hillslope$q_uz_sz <- hillslope$s_uz + hillslope$q_rz_uz - tilde_uz
            hillslope$s_uz <- tilde_uz

            ## Step 4: Solve saturated zone
            ## if mass check compute theinitial mass

            if( mass_check ){
                mass_s_sz <-  -sum(hillslope$s_sz*hillslope$area)
            }

            ## move current states to values for start of the time step
            hillslope$sum_l_sz_in_t <- hillslope$sum_l_sz_in # total inflow at start of time step
            hillslope$l_sz_t <- hillslope$l_sz # total outflow at start of time step
            channel$sum_l_sz_in_t <- channel$sum_l_sz_in
            ## evaluate the values of the initial components of the kinematic solution
            hillslope$Q_minus_t <- pmin( hillslope$sum_l_sz_in_t, hillslope$l_szmax ) # inflow to saturated zone at start of timestep
            hillslope$Q_plus_t <- pmin( hillslope$l_sz_t, hillslope$l_szmax ) # outflow at start of timestep

            ## compute velocity estimate and kinematic parameters
            cbar <- (hillslope$l_szmax/hillslope$m)*
                exp(- hillslope$s_sz / hillslope$m)
            lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/hillslope$delta_x
            lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/hillslope$delta_x
            
            ## update flows by looping through bands
            for(idx in sqnc$sz){
                hillslope$sum_l_sz_in[idx] <- lateral_flux$sz[ hillslope$id[idx] ] / hillslope$area[idx] ## current inflow in m/s
                
                ## compute the new state value
                hillslope$Q_minus_tDt[idx] <- pmin( hillslope$sum_l_sz_in[idx],hillslope$l_szmax[idx] ) # current inflow to saturated zone
                
                k <- lambda_prime[idx] * hillslope$Q_plus_t[idx] +
                    (1-lambda_prime[idx]) * hillslope$Q_minus_t[idx] +
                    cbar[idx]*hillslope$q_uz_sz[idx]/hillslope$delta_x[idx]
                
                hillslope$l_sz[idx] <- pmin( (k - (1-lambda[idx])*hillslope$Q_minus_tDt[idx])/lambda[idx] , hillslope$l_szmax[idx] )
                
                if( any(hillslope$l_sz[idx]<0) ){
                    warning("Negative flow in kinematic solutions, consider revising weights")
                    hillslope$l_sz[idx] <- pmax(hillslope$l_sz[idx],0)
                }
                
                for(ii in idx){
                    
                    lateral_flux$sz[ hillslope$sz_dir[[ii]]$idx ] <-
                        lateral_flux$sz[ hillslope$sz_dir[[ii]]$idx ] +
                        hillslope$sz_dir[[ii]]$frc * hillslope$l_sz[ii] * hillslope$area[ii]
                    
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
            channel$sum_l_sz_in <- lateral_flux$sz[channel$id] / channel$area
            channel$s_ch <-  channel$s_ch +
                (lateral_flux$sf[channel$id] / channel$area ) +
                ts$sub_step*(channel$sum_l_sz_in_t + channel$sum_l_sz_in)/2
            
            ## step 5 - correct the stores for saturation flows
            saturated_index <- hillslope$s_sz <= 0
            hillslope$q_uz_sf <- hillslope$s_uz*saturated_index
            hillslope$s_uz <- hillslope$s_uz * !saturated_index
            hillslope$s_sf <- hillslope$s_sf +
                hillslope$q_rz_sf +
                hillslope$q_sz_sf + hillslope$q_uz_sf
            
            ## mass check for iteration
            if( mass_check ){
                vol_ch_sf <- lateral_flux$sf[ channel$id ]
                vol_ch_sz <- ( (channel$s_ch-ch0$s_ch)*channel$area ) - vol_ch_sf
                browser()
                mass_errors[ ((it-1)*ts$n_sub_step) + inner,] <-
                    c(it,inner,
                      sum(hs0$s_sf*hs0$area)+
                      sum(ch0$s_ch*ch0$area) +
                      sum(hillslope$q_rz_sf*hillslope$area) +
                      sum(hillslope$q_sz_sf*hillslope$area) -
                      sum(hillslope$q_sf_rz*hillslope$area) -
                      sum(hillslope$s_sf*hillslope$area) -
                      sum(vol_ch_sf),
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
                      mass_s_sz +
                      sum(hillslope$q_uz_sz*hillslope$area) -
                      sum(hillslope$q_sz_sf*hillslope$area) -
                      - sum(hillslope$s_sz*hillslope$area) -
                      sum(vol_ch_sz)
                      )
                ##print( mass_errors[ ((it-1)*ts$n_sub_step) + inner,] )
            }
        } ## end of sub_step loop
        

        ## step 6 - channel inflow - at the moment a volume / area
        channel_inflow[it,] <- channel$area *
            ( channel$s_ch + channel$p*ts$step ) /
            (ts$step)
        
        ## handle returning states
        if( return_states$flag ){
            return_states[[index(obs_data)[it]]] <- get_states(hillslope,"dynatop","hillslope")
        }
        


    } ## end of timestep loop


    ## copy states back into model
    tmp <- get_states(hillslope,"dynatop","hillslope")
    nm <- c("id",setdiff( names(model$hillslope),names(tmp)))
    model$hillslope <- merge(model$hillslope[,nm],tmp,by="id",all=TRUE)
    tmp <- get_states(channel,"dynatop","channel")
    nm <- c("id",setdiff( names(model$channel),names(tmp)))
    model$channel <- merge(model$channel[,nm],tmp,by="id",all=TRUE)



    out <- list(model = model,
                channel_input = channel_inflow)

    if(mass_check){
        #browser()
        mass_errors <- as.data.frame(mass_errors)
        mass_errors[,'DateTime'] <- obs_index[mass_errors[,'DateTime']]
        out[['mass_errors']] <- mass_errors
    }
    if( return_states$flag ){
        return_states <- return_states[setdiff(names(return_states),c("idx","flag"))]
        out[['return_states']] <- return_states
    }
    return( out )
}

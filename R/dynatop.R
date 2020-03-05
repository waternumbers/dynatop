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
    ## this creates hillslope, channel and sqnc
    list2env(convert_form(model),as.environment(-1))
    
    ## work out band sequences   
    sqnc <- list(sf_band=list(),sz_band=list())
    for(ii in names(sqnc)){      
        tmp <- sort(unique(c(model$hillslope$attr[[ii]],
                             model$channel$attr[[ii]]))) # sorted list of unique bands
        for(jj in 1:length(tmp)){
            sqnc[[ii]][[jj]] <- list(
                hillslope = which(model$hillslope$attr[[ii]] == tmp[jj]),
                channel = which(model$channel$attr[[ii]] == tmp[jj]))
        }
    }
    
    ## storage for lateral fluxes stored as volumes
    lvol <- list(sf = rep(0,max(c(model$hillslope$attr$id,model$channel$attr$id))),
                 sz = rep(0,max(c(model$hillslope$attr$id,model$channel$attr$id)))) 
    
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
    }
    
                              
    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")
    
    for(it in 1:nrow(obs_data)){
        ## set the inputs to the hillslope and channel
        ## set as rate m/s
        for(ii in c("hillslope","channel")){
            model[[ii]]$input$p <- obs_data[it,model[[ii]]$attr$precip]/ts$step
            model[[ii]]$input$e_t <- obs_data[it,model[[ii]]$attr$pet]/ts$step
        }
        
        model$channel$flux$s_ch[] <- 0
       
        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## stor previous version of model if need mass check
            if( mass_check ){ model0 <- model }
                        
            ## Step 1: Distribute any surface storage downslope            
            lvol$sf[] <- 0 # remove fluxes from previous time step since not needed - should be over written but...
            
            for(bnd in sqnc$sf_band){ ## loop all bands of surface 
                ## hillslope
                idx <- bnd$hillslope ## index of hillslope HSUs
                if(any(idx)){
                    for(ii in idx){
                        model$hillslope$input$l_sf[ii] <- sum( model$sf[[ model$hillslope$attr$id[ii] ]]$x *
                                                               lvol$sf[ model$sf[[ model$hillslope$attr$id[ii] ]]$j])
                    }
                    
                    model$hillslope$input$l_sf[idx] <- model$hillslope$input$l_sf[idx]/model$hillslope$attr$area[idx] # inflow in m depth accrued over sub step
                    
                    ## compute the new state value
                    tilde_sf <- fode( model$hillslope$input$l_sf[idx]/ts$sub_step, 1/model$hillslope$param$t_sf[idx],
                                     model$hillslope$state$s_sf[idx],ts$sub_step)
                    ## work out out flow
                    model$hillslope$flux$l_sf[idx] <- model$hillslope$state$s_sf[idx] + model$hillslope$input$l_sf[idx] - tilde_sf
                    model$hillslope$state$s_sf[idx] <- tilde_sf
                    lvol$sf[ model$hillslope$attr$id[idx] ]  <- model$hillslope$flux$l_sf[idx]*model$hillslope$attr$area[idx]
                }
                

                      
                ## channel
                idx <- bnd$channel #which(model$channel$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    ##model$channel$input$lsf[idx] <- bnd$Fch%*%lvol$sf
                    for(ii in idx){
                        model$channel$input$l_sf[ii] <- sum( model$sf[[ model$channel$attr$id[ii] ]]$x *
                                                             lvol$sf[ model$sf[[ model$channel$attr$id[ii] ]]$j])
                    }
                    model$channel$input$l_sf[idx] <- model$channel$input$l_sf[idx]/model$channel$attr$area[idx]
                    ## compute the new store of channel input in depth
                    model$channel$flux$s_ch[idx] <- model$channel$flux$s_ch[idx] + model$channel$input$l_sf[idx]
                    lvol$sf[ model$channel$attr$id[idx] ]  <- 0
                }
            }
            
            if( mass_check ){ mc_s_ch_sf <- model$channel$flux$s_ch }
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
            model$hillslope$flux$q_sf_rz <- pmin( model$hillslope$param$q_sfmax*ts$sub_step,model$hillslope$state$s_sf )
            model$hillslope$state$s_sf <- model$hillslope$state$s_sf - model$hillslope$flux$q_sf_rz
            
            ## solve ODE
            
            tilde_rz <- fode( model$hillslope$input$p + (model$hillslope$flux$q_sf_rz/ts$sub_step),
                             model$hillslope$input$e_t/model$hillslope$param$s_rzmax,
                             model$hillslope$state$s_rz,ts$sub_step )

            ## work out actual evapotranspiration by mass balance
            model$hillslope$flux$e_t <- model$hillslope$state$s_rz + model$hillslope$input$p*ts$sub_step - tilde_rz
            ## new storage value
            model$hillslope$state$s_rz <- pmin(tilde_rz,model$hillslope$param$s_rzmax)
            
            ## split root zone flow
            tmp <- tilde_rz - model$hillslope$state$s_rz
            saturated_index <- model$hillslope$state$s_sz <= 0 # which areas are saturated
            model$hillslope$flux$q_rz_sf <- tmp * saturated_index
            model$hillslope$flux$q_rz_uz <- tmp * !saturated_index
            
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
            tilde_uz <- fode( model$hillslope$flux$q_rz_uz/ts$sub_step,
                             1 / (model$hillslope$param$t_d * model$hillslope$state$s_sz),
                             model$hillslope$state$s_uz,ts$sub_step )
            
            model$hillslope$flux$q_uz_sz <- model$hillslope$state$s_uz + model$hillslope$flux$q_rz_uz - tilde_uz
            model$hillslope$state$s_uz <- tilde_uz
            
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
            ## compute some required values for the hillslope
            model$hillslope$flux$l_sz_in_t <- model$hillslope$state$l_sz_in # total inflow at start of time step
            model$hillslope$flux$Q_minus_t <- pmin( model$hillslope$state$l_sz_in,model$hillslope$state$l_szmax ) # inflow to saturated zone at start of timestep
            model$hillslope$flux$Q_plus_t <- pmin( model$hillslope$state$l_sz,model$hillslope$state$l_szmax ) # outflow at start of timestep
            ## compute some required values for the channel
            model$channel$flux$l_sz_in_t <- model$channel$state$l_sz_in # total inflow at start of time step
            
            ## update flows by looping through bands
            for(bnd in sqnc$sz_band){
                ## hillslope
                idx <- bnd$hillslope
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        model$hillslope$state$l_sz_in[ii] <- sum(model$sz[[ model$hillslope$attr$id[ii] ]]$x *
                                                                 lvol$sz[ model$sz[[ model$hillslope$attr$id[ii] ]]$j ])
                    }
                    model$hillslope$state$l_sz_in[idx] <- model$hillslope$state$l_sz_in[idx]/model$hillslope$attr$area[idx] ## current inflow in m/s

                    ## compute the new state value
                    Q_minus_tDt <- pmin( model$hillslope$state$l_sz_in[idx],model$hillslope$state$l_szmax[idx] ) # current inflow to saturated zone
                    
                    qbar <- (model$hillslope$flux$Q_minus_t[idx] + Q_minus_tDt + model$hillslope$flux$Q_plus_t[idx])/3
                    cbar <- (qbar*model$hillslope$attr$delta_x[idx])/model$hillslope$param$m[idx]
                    
                    lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/model$hillslope$attr$delta_x[idx]
                    lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/model$hillslope$attr$delta_x[idx]

                    k <- lambda_prime * model$hillslope$flux$Q_plus_t[idx] +
                        (1-lambda_prime) * model$hillslope$flux$Q_minus_t[idx] +
                        cbar*model$hillslope$flux$q_uz_sz[idx]/model$hillslope$attr$delta_x[idx]

                    model$hillslope$state$l_sz[idx] <- pmin( (k - (1-lambda)*Q_minus_tDt)/lambda , model$hillslope$state$l_szmax[idx] )
                    lvol$sz[ model$hillslope$attr$id[idx] ] <- model$hillslope$state$l_sz[idx]*model$hillslope$attr$area[idx]

                }
                ## channel
                idx <- bnd$channel
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        model$channel$state$l_sz_in[ii] <- sum(model$sz[[ model$channel$attr$id[ii] ]]$x *
                                                               lvol$sz[model$sz[[ model$channel$attr$id[ii] ]]$j])
                    }
                    model$channel$state$l_sz_in[idx] <- model$channel$state$l_sz_in[idx]/model$channel$attr$area[idx]
                    lvol$sz[ model$channel$attr$id[idx] ]  <- 0
                }
            }

            
            ## update volumes in hillslope
            tilde_sz <- model$hillslope$state$s_sz +
                ts$sub_step*(model$hillslope$flux$Q_plus_t + model$hillslope$state$l_sz)/2 -
                ts$sub_step*(model$hillslope$flux$l_sz_in_t + model$hillslope$state$l_sz_in)/2 -
                model$hillslope$flux$q_uz_sz
            
            model$hillslope$state$s_sz <- pmax(0,tilde_sz)
            model$hillslope$flux$q_sz_sf <- model$hillslope$state$s_sz - tilde_sz

            ## update volume of inflow to channel
            model$channel$flux$s_ch <- model$channel$flux$s_ch +
                ts$sub_step*(model$channel$flux$l_sz_in_t + model$channel$state$l_sz_in)/2
            

            
            ## step 5 - correct the stores for saturation flows
            ##browser()
            ## if( any(model$hillslope$state$s_sz <= 1e-5) ){
            ##     browser()
            ## }
            
            saturated_index <- model$hillslope$state$s_sz <= 0
            model$hillslope$flux$q_uz_sf <- model$hillslope$state$s_uz*saturated_index
            model$hillslope$state$s_uz <- model$hillslope$state$s_uz * !saturated_index
            model$hillslope$state$s_sf <- model$hillslope$state$s_sf +
                model$hillslope$flux$q_rz_sf +
                model$hillslope$flux$q_sz_sf + model$hillslope$flux$q_uz_sf

            ## mass check for iteration
            if( mass_check ){
                mass_errors[ ((it-1)*ts$n_sub_step) + inner,] <-
                    c(it,inner,
                      sum(model0$hillslope$state$s_sf*model0$hillslope$attr$area)+
                      sum(model0$channel$flux$s_ch*model0$channel$attr$area) +
                      sum(model$hillslope$flux$q_rz_sf*model$hillslope$attr$area) +
                      sum(model$hillslope$flux$q_sz_sf*model$hillslope$attr$area) -
                      sum(model$hillslope$flux$q_sf_rz*model$hillslope$attr$area) -
                      sum(model$hillslope$state$s_sf*model$hillslope$attr$area) -
                      sum(mc_s_ch_sf*model$channel$attr$area),
                      sum((model0$hillslope$state$s_rz +
                           model$hillslope$input$p*ts$sub_step +
                           model$hillslope$flux$q_sf_rz -
                           model$hillslope$flux$q_rz_sf -
                           model$hillslope$flux$q_rz_uz -
                           model$hillslope$flux$e_t -
                           model$hillslope$state$s_rz)*
                          model0$hillslope$attr$area),
                      sum((model0$hillslope$state$s_uz +
                           model$hillslope$flux$q_rz_uz -
                           model$hillslope$flux$q_uz_sz -
                           model$hillslope$flux$q_uz_sf -
                           model$hillslope$state$s_uz)*model0$hillslope$attr$area),
                      -sum(model0$hillslope$state$s_sz*model0$hillslope$attr$area)+
                      sum(mc_s_ch_sf*model$channel$attr$area) +
                      sum(model$hillslope$flux$q_uz_sz*model$hillslope$attr$area) -
                      sum(model$hillslope$flux$q_sz_sf*model$hillslope$attr$area) -
                      -sum(model$hillslope$state$s_sz*model$hillslope$attr$area) -
                      sum(model$channel$flux$s_ch*model$channel$attr$area)
                      )
                #print( mass_errors[ ((it-1)*ts$n_sub_step) + inner,] )
            }     
            ##browser()
            ## step 6 - channel inflow - at the moment a volume / area
            channel_inflow[it,] <- model$channel$attr$area *
                ( model$channel$flux$s_ch + model$channel$input$p*ts$step ) /
                (ts$step)

            ## handle returning states
            if( length(return_states)>0 && return_states$idx[it] ){
                return_states[[index(obs_data)[it]]] <- sim_to_store(model,with_matrix=FALSE)
            }
            
            
        }
        
    } ## end of timestep loop
    
    
    out <- list(model = sim_to_store(model),
                channel_input = channel_inflow)
    
    if(mass_check){
        #browser()
        mass_errors <- as.data.frame(mass_errors)
        mass_errors[,'DateTime'] <- index(obs)[mass_errors[,'DateTime']]
        out[['mass_errors']] <- mass_errors
    }
    if( length(return_states)>0 ){
        return_states <- return_states[setdiff(names(return_states),"idx")]
        out[['return_states']] <- return_states
    }
    return( out )
}

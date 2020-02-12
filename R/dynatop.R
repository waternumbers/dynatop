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
                                theta=0.7)){

    #browser()
    ## check the model
    input_series <- check_model(model)
    ## check input and get model timestep
    ts <- check_obs(obs_data,input_series,
                    sim_time_step)

    if(use_states){ # then just take states from the model object
        mdl <- model$hsu
    }else{ # initialise the model
        mdl <- initialise(model,initial_recharge)
    }

    ## work out bands
    bands <- list(sf = sort(unique(c(model$hillslope$attr$sf_band,
                                     model$channel$attr$sf_band))),
                  sz = sort(unique(c(model$hillslope$attr$sz_band,
                                     model$channel$attr$sz_band)))
                  )

    ## storage for lateral fluxes stored as volumes
    lvol <- list(sf = rep(0,max(c(model$hillslope$attr$id,model$channel$attr$id))),
                 sz = rep(0,max(c(model$hillslope$attr$id,model$channel$attr$id)))) 

    ## TODO populate lval$sz with in initial output$lsz values
    
    ## simple function for solving ODE
    fode <- function(a,b,x0,t){
        b <- pmax(b,1e-10)
        unname( x0*exp(-b*t) + (a/b)*(1-exp(-b*t)) )
    }

    ## initialise the output
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(channel$id)), match.to=obs_data )
    names(channel_inflow) <- channel$id

    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){
        #browser()
        print(it)
        
        ## set the inputs to the hillslope and channel
        for(ii in c("hillslope","channel")){
            mdl[[ii]]$input$precip <- obs_data[it,mdl[[ii]]$series$precip]/ts$step
            mdl[[ii]]$input$pet <- obs_data[it,mdl[[ii]]$series$pet]/ts$step
        }
        
        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## Step 1: Distribute any surface storage downslope
            for(bnd in bands$sf){
                ## hillslope
                idx <- which(mdl$hillslope$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        mdl$hillslope$input$lsf[ii] <- sum(lvol$sf[mdl$sf[[ii]]$j]) #TODO wrong index?
                    }
                    mdl$hillslope$input$lsf[idx] <- mdl$hillslope$input$lsf[idx]/mdl$hillslope$attr$area[idx]
                    ## compute the new state value
                    tilde_sf <- fode( mdl$hillslope$input$lsf[idx]/ts$sub_step, 1/mdl$hillslope$param$tsf[idx],
                                     mdl$hillslope$state$sf[idx],ts$sub_step)
                    mdl$hillslope$output$lsf[idx] <- mdl$hillslope$state$sf[idx] + mdl$hillslope$input$lsf[idx] - tilde_sf
                    mdl$hillslope$state$sf[idx] <- tilde_sf
                    lvol$sf[ mdl$hillslope$attr$id[idx] ]  <- mdl$hillslope$output$lsf[idx]*mdl$hillslope$attr$area[idx]
                }
                ## channel
                idx <- which(mdl$channel$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        mdl$channel$input$lsf[ii] <- sum(lvol$sf[mdl$sf[[ii]]$j])
                    }
                    mdl$channel$input$lsf[idx] <- mdl$channel$input$lsf[idx]/mdl$channel$attr$area[idx]
                    ## compute the new state value
                    mdl$channel$state$sch[idx] <- dl$channel$input$lsf[idx]
                    lvol$sf[ mdl$hillslope$attr$id[idx] ]  <- 0
                }
            }

            ## Step 2: solve the root zone for hillslope elements
            
            ## evaluate max integral of flow to rootzone
            mdl$hillslope$flux$qsf_rz <- pmin( mdl$hillslope$param$qsf_max*ts$sub_step,mdl$hillslope$state$ssf )

            ## solve ODE
            tilde_rz <- fode( mdl$hillslope$input$precip + (mdl$hillslope$flux$qsf_rz/ts$sub_step),
                             mdl$hillslope$input$pet/mdl$hillslope$param$srz_max,
                             mdl$hillslope$state$rz,ts$sub_step )
            ## new storage value
            mdl$hillslope$state$rz <- pmin(tilde_rz,mdl$hillslope$param$srz_max)
            ## split root zone flow
            tmp <- tilde_rz - mdl$hillslope$state$rz
            saturated_index <- mdl$hillslope$state$sz <= 0 # which areas are saturated
            mdl$hillslope$flux$qrz_ex <- tmp * saturated_index
            mdl$hillslope$flux$qrz_uz <- tmp * !saturated_index

            ## Step 3: Unsaturated zone

            ## solve ODE
            tilde_uz <- fode( mdl$hillslope$flux$qrz_uz/ts$sub_step,
                             1 / (mdl$hillslope$param$td * mdl$hillslope$state$sz),
                             mdl$hillslope$state$uz,ts$sub_step )
            mdl$hillslope$flux$quz_sz <- mdl$hillslope$state$uz + mdl$hillslope$flux$qrz_uz - tilde_uz
            mdl$hillslope$state$uz <- tilde_uz


            ## Step 4: Solve saturated zone
            for(bnd in bands$sz){
                ## hillslope
                idx <- which(mdl$hillslope$attr$sz_band == bnd)
                if(any(idx)){
                    ## update the input
                    Q_minus_t <- pmin( mdl$hillslope$input$lsz[idx],mdl$hillslope$state$lsz_max[idx] )
                    for(ii in idx){
                        mdl$hillslope$input$lsz[ii] <- sum(lvol$sz[mdl$sz[[ii]]$j]) #TODO wrong idenx?
                    }
                    mdl$hillslope$input$lsz[idx] <- mdl$hillslope$input$lsz[idx]/mdl$hillslope$attr$area[idx]
                    ## compute the new state value
                    Q_minus_tDt <- pmin( mdl$hillslope$input$lsz[idx],mdl$hillslope$state$lsz_max[idx] )
                    Q_plus_t <- mdl$hillslope$output$lsz[idx]

                    qbar <- (Q_minus_t + Q_minus_tDt + Q_plus_t)/3
                    cbar <- (qbar*mdl$hillslope$attr$delta_x[idx])/mdl$hillslope$param$m[idx]
                    
                    lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/mdl$hillslope$attr$delta_x
                    lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/mdl$hillslope$attr$delta_x

                    k <- lambda_prime*Q_plus_t +
                        (1-lambda_prime)*Q_minus_t +
                        cbar*mdl$hillslope$flux$quz_sz[idx]/mdl$hillslope$attr$delta_x[idx]

                    ## GOT HERE
                    
                    sz$lsz[ii] <- min( (k - (1-lambda)*lsz_in_sat)/lambda , mdl$hillslope$lsz_max[ii] )
                    
                    tilde_sf <- fode( mdl$hillslope$input$lsf[idx]/ts$sub_step, 1/mdl$hillslope$param$tsf[idx],
                                     mdl$hillslope$state$sf[idx],ts$sub_step)
                    mdl$hillslope$output$lsf[idx] <- mdl$hillslope$state$sf[idx] + mdl$hillslope$input$lsf[idx] - tilde_sf
                    mdl$hillslope$state$sf[idx] <- tilde_sf
                    lvol$sf[ mdl$hillslope$attr$id[idx] ]  <- mdl$hillslope$output$lsf[idx]*mdl$hillslope$attr$area[idx]
                }
                ## channel
                idx <- which(mdl$channel$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        mdl$channel$input$lsf[ii] <- sum(lvol$sf[mdl$sf[[ii]]$j])
                    }
                    mdl$channel$input$lsf[idx] <- mdl$channel$input$lsf[idx]/mdl$channel$attr$area[idx]
                    ## compute the new state value
                    mdl$channel$state$sch[idx] <- dl$channel$input$lsf[idx]
                    lvol$sf[ mdl$hillslope$attr$id[idx] ]  <- 0
                }
            }


            
            for(ii in sz$order){
                ## inflow from upstream

                #idx <- sz$hs_parent[[ii]]$j #[,1]
                #w <- sz$hs_parent[[ii]]$x #[,2]

                #sz$lsz_in[ii] <- sum(w * sz$lsz[idx]) # latest inflow
                sz$lsz_in[ii] <- sum(sz$hs_parent[[ii]]$x * sz$lsz[sz$hs_parent[[ii]]$j]) # latest inflow

                lsz_in_sat <- max(sz$lsz_in[ii],mdl$hs$lsz_max[ii]) # constrained to saturated zone max flow

                qbar <- (mdl$hs$lsz[ii]+lsz_in_sat_prev[ii]+lsz_in_sat + max(lsz_in_sat,mdl$hs$lsz[ii]))/4

                # TO DO check the use or not of area in the following
                cbar <- (qbar*mdl$hs$delta_x[ii])/mdl$hs$m[ii]
                lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/mdl$hs$delta_x
                lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/mdl$hs$delta_x

                # not ii enough

                k <- lambda_prime*mdl$hs$lsz[ii] +
                    (1-lambda_prime)*lsz_in_sat_prev[ii] +
                    cbar*mdl$hs$quz_sz[ii]/mdl$hs$delta_x[ii]
                sz$lsz[ii] <- min( (k - (1-lambda)*lsz_in_sat)/lambda , mdl$hs$lsz_max[ii] )
            }

            #browser()
            ## work out integral fluz
            mdl$hs$qsz <- ts$sub_step*(sz$lsz + mdl$hs$lsz)/2
            tilde_sz <- mdl$hs$sz + mdl$hs$qsz - ts$sub_step*(sz$lsz_in + mdl$hs$lsz_in)/2 - mdl$hs$quz_sz
            mdl$hs$sz <- pmax(0,tilde_sz)
            mdl$hs$qsz_ex <- mdl$hs$sz - tilde_sz

            ## update states
            mdl$hs$lsz <- sz$lsz
            mdl$hs$lsz_in <- sz$lsz_in

            ## step 5 - correct the stores
            saturated_index <- mdl$hs$sz <= 0
            mdl$hs$ex <- mdl$hs$ex + mdl$hs$qrz_ex +
                (mdl$hs$qsz_ex + mdl$hs$uz)*saturated_index
            mdl$hs$uz <- mdl$hs$uz * !saturated_index

            ## step 6 - channel inflow - at the moment a volume / area
            channel_inflow[it,] <- channel_inflow[it,] + ex$dex[channel$id] + obs_data[it,channel$precip_series]*ts$sub_step
            for(ii in 1:length(channel$id)){
                idx <- sz$ch_parent[[ii]][,1]
                w <- sz$ch_parent[[ii]][,2]
                channel_inflow[it,ii]  <-  channel_inflow[it,ii] +sum( w*mdl$hs$area[idx]*mdl$hs$qsz[idx]/channel$area[ii] )
            }


        }

        ## step 7: Done through the steps above
        ## put channel inflow in output matric
        channel_inflow[it,] <- channel_inflow[it,]*channel$area / (3600*ts$step)

    } ## end of timestep loop

    model$states <- list(hillslope=mdl$hs,
                         channel=channel,
                         ex=ex,
                         sz=sz)

    return( list(model=model,
                 channel_input = channel_inflow ) )
}

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
    bands <- list(sf=list(),sz=list())
    tmp <- sort(unique(c(mdl$hillslope$attr$sf_band,
                         mdl$channel$attr$sf_band)))
    for(ii in 1:length(tmp)){
        bands$sf[[ii]] <- list(hillslope = which(mdl$hillslope$attr$sf_band == tmp[ii]),
                               channel = which(mdl$channel$attr$sf_band == tmp[ii]))
    }
    tmp <- sort(unique(c(mdl$hillslope$attr$sz_band,
                         mdl$channel$attr$sz_band)))
    for(ii in 1:length(tmp)){
        bands$sz[[ii]] <- list(hillslope = which(mdl$hillslope$attr$sz_band == tmp[ii]),
                               channel = which(mdl$channel$attr$sz_band == tmp[ii]))
    }
    
    ## storage for lateral fluxes stored as volumes
    lvol <- list(sf = rep(0,max(c(mdl$hillslope$attr$id,mdl$channel$attr$id))),
                 sz = rep(0,max(c(mdl$hillslope$attr$id,mdl$channel$attr$id)))) 

    ## TODO populate lval$sz with in initial output$lsz values
    
    ## simple function for solving ODE
    fode <- function(a,b,x0,t){
        b <- pmax(b,1e-10)
        unname( x0*exp(-b*t) + (a/b)*(1-exp(-b*t)) )
    }

    ## initialise the output
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(mdl$channel$attr$id)), match.to=obs_data )
    names(channel_inflow) <- mdl$channel$attr$id

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
            #print("Step 1")
            lvol$sf[] <- 0
            for(bnd in bands$sf){
                ## hillslope
                idx <- bnd$hillslope #which(mdl$hillslope$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        ##jj <- mdl$hillslope$attr$id[ii]
                        ##mdl$hillslope$input$lsf[ii] <- sum( mdl$sf[[ jj]]$x * lvol$sf[ mdl$sf[[jj]]$j])
                        mdl$hillslope$input$lsf[ii] <- sum( mdl$sf[[ mdl$hillslope$attr$id[ii]]]$x *
                                                            lvol$sf[ mdl$sf[[mdl$hillslope$attr$id[ii]]]$j])
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
                idx <- bnd$channel #which(mdl$channel$attr$sf_band == bnd)
                if(any(idx)){
                    ## update the input
                    for(ii in idx){
                        ##jj <-  mdl$channel$attr$id[ii]
                        ##mdl$channel$input$lsf[ii] <- sum(mdl$sf[[jj]]$x * lvol$sf[mdl$sf[[jj]]$j])
                        mdl$channel$input$lsf[ii] <- sum( mdl$sf[[ mdl$channel$attr$id[ii]]]$x *
                                                            lvol$sf[ mdl$sf[[mdl$channel$attr$id[ii]]]$j])
                    }
                    mdl$channel$input$lsf[idx] <- mdl$channel$input$lsf[idx]/mdl$channel$attr$area[idx]
                    ## compute the new state value
                    mdl$channel$state$sch[idx] <- mdl$channel$input$lsf[idx]
                    lvol$sf[ mdl$hillslope$attr$id[idx] ]  <- 0
                }
            }

            ## Step 2: solve the root zone for hillslope elements
            #print("Step 2")
            ## evaluate max integral of flow to rootzone
            mdl$hillslope$flux$qsf_rz <- pmin( mdl$hillslope$param$qsf_max*ts$sub_step,mdl$hillslope$state$sf )

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
            #print("Step 3")
            ## solve ODE
            tilde_uz <- fode( mdl$hillslope$flux$qrz_uz/ts$sub_step,
                             1 / (mdl$hillslope$param$td * mdl$hillslope$state$sz),
                             mdl$hillslope$state$uz,ts$sub_step )
            mdl$hillslope$flux$quz_sz <- mdl$hillslope$state$uz + mdl$hillslope$flux$qrz_uz - tilde_uz
            mdl$hillslope$state$uz <- tilde_uz


            ## Step 4: Solve saturated zone
            #print("Step 4")
            mdl$hillslope$flux$lsz_t <- mdl$hillslope$input$lsz
            mdl$hillslope$flux$Q_minus_t <- pmin( mdl$hillslope$input$lsz,mdl$hillslope$state$lsz_max )
            mdl$hillslope$flux$Q_plus_t <- mdl$hillslope$output$lsz
            mdl$channel$flux$lsz_t <- mdl$channel$input$lsz
            for(bnd in bands$sz){
                ## hillslope
                #print("hillslope")
                #browser()
                idx <- bnd$hillslope#which(mdl$hillslope$attr$sz_band == bnd)
                if(any(idx)){
                    #print(idx)
                    ## update the input
                    for(ii in idx){
                        mdl$hillslope$input$lsz[ii] <- sum(mdl$sf[[mdl$hillslope$attr$id[ii]]]$x *
                                                           lvol$sz[mdl$sz[[mdl$hillslope$attr$id[ii]]]$j])
                    }
                    mdl$hillslope$input$lsz[idx] <- mdl$hillslope$input$lsz[idx]/mdl$hillslope$attr$area[idx]

                    ## compute the new state value
                    mdl$hillslope$flux$Q_minus_tDt[idx] <- pmin( mdl$hillslope$input$lsz[idx],mdl$hillslope$state$lsz_max[idx] )

                    qbar <- (mdl$hillslope$flux$Q_minus_t[idx] + mdl$hillslope$flux$Q_minus_tDt[idx] + mdl$hillslope$flux$Q_plus_t[idx])/3
                    cbar <- (qbar*mdl$hillslope$attr$delta_x[idx])/mdl$hillslope$param$m[idx]
                    
                    lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/mdl$hillslope$attr$delta_x[idx]
                    lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/mdl$hillslope$attr$delta_x[idx]

                    k <- lambda_prime*mdl$hillslope$flux$Q_plus_t[idx] +
                        (1-lambda_prime)*mdl$hillslope$flux$Q_minus_t[idx] +
                        cbar*mdl$hillslope$flux$quz_sz[idx]/mdl$hillslope$attr$delta_x[idx]

                    mdl$hillslope$output$lsz[idx] <- pmin( (k - (1-lambda)*mdl$hillslope$flux$Q_minus_tDt[idx])/lambda , mdl$hillslope$state$lsz_max[idx] )
                    lvol$sz[ mdl$hillslope$attr$id[idx] ] <- mdl$hillslope$output$lsz[idx]*mdl$hillslope$attr$area[idx]

                }
                ## channel
                #print("channel")
                idx <- bnd$channel #which(mdl$channel$attr$sz_band == bnd)
                if(any(idx)){
                    #browser()
                    ## update the input
                    
                    for(ii in idx){
                        ## jj <-  mdl$hillslope$attr$id[ii]
                        
                        mdl$channel$input$lsz[ii] <- sum(mdl$sz[[ mdl$hillslope$attr$id[ii] ]]$x *
                                                         lvol$sz[mdl$sz[[ mdl$hillslope$attr$id[ii] ]]$j])
                    }
                    mdl$channel$input$lsz[idx] <- mdl$channel$input$lsz[idx]/mdl$channel$attr$area[idx]
                    ## compute the new state value
                    #browser()
                    mdl$channel$state$sch[idx] <- mdl$channel$state$sch[idx] +
                        ts$sub_step*(mdl$channel$flux$lsz_t[idx] + mdl$channel$input$lsz[idx])/2
                    lvol$sf[ mdl$hillslope$attr$id[idx] ]  <- 0
                }
            }

            tilde_sz <- mdl$hillslope$state$sz +
                ts$sub_step*(mdl$hillslope$flux$Q_plus_t + mdl$hillslope$output$lsz)/2 -
                ts$sub_step*(mdl$hillslope$flux$lsz_t + mdl$hillslope$input$lsz)/2 -
                mdl$hillslope$flux$quz_sz
            
            mdl$hillslope$state$sz <- pmax(0,tilde_sz)
            mdl$hillslope$flux$qsz_ex <- mdl$hillslope$state$sz - tilde_sz


            ## step 5 - correct the stores
            #print("Step 5")
            saturated_index <- mdl$hillslope$state$sz <= 0
            mdl$hillslope$state$sf <- mdl$hillslope$state$sf +
                mdl$hillslope$flux$qrz_ex +
                (mdl$hillslope$flux$qsz_ex + mdl$hillslope$state$uz)*saturated_index
            mdl$hillslope$state$uz <- mdl$hillslope$state$uz * !saturated_index
        }
        
        ## step 6 - channel inflow - at the moment a volume / area
        channel_inflow[it,] <- mdl$channel$attr$area *
            ( mdl$channel$state$sch + mdl$channel$input$precip*ts$step ) /
            (3600*ts$step)

    } ## end of timestep loop

    model$hsu <- mdl

    return( list(model=model,
                 channel_input = channel_inflow ) )
}

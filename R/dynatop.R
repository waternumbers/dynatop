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
dynatop <- function(model,obs_data,initial_recharge=NA,sim_time_step=NULL, use_states=FALSE){
    # check input and get model timestep
    ts <- check_obs(obs_data,
                    unique( c(unlist(model$hillslope[,c("precip_input","pet_input")]),
                              unlist(model$channel[,c("precip_input","pet_input")]))),
                           sim_time_step)

    if(use_states){ # then just take states from the model object
        hillslope <- model$states$hillslope
        channel <- model$states$channel
    }else{ # initialise the model
        ## check model is valid - fails if not
        check_model(model)

        ## initialise the properties and states of the hillslope hru - include check on initial recharge
        hillslope <- initialise_hillslope(model,initial_recharge)

        ## initialise the properties and states of the channel hru
        channel <- initialise_channel(model)
    }
    

    ## create the common parts of the surface excess solution
    K_ex <- diag(1/c(hillslope$area,channel$area)) %*%
        cbind( rbind( (model$Wex-diag(nrow(model$Wex))), model$Fex),
              matrix(0,length(hillslope$id)+length(channel$id),length(channel$id) ) # no flow from channel
              ) %*% diag(c(hillslope$area,channel$area)) %*%
        diag(c(1/hillslope$tex,rep(0,length(channel$id))))
    ex_eigen <- eigen_routing_setup(K_ex)
    ex_in <- ex_out <- rep(0,nrow(K_ex)) # preassign vector for initial condition storages
    ex_idx <- 1:length(hillslope$id) # index of hillslope elements in the solution

    ## Common parts for the saturated routing
    K_sz <- diag(1/hillslope$area) %*% (model$Wsat-diag(nrow(model$Wsat))) %*%
               diag(hillslope$area)
    FA_sz <- model$Fsat%*%diag(hillslope$area)
    #WA_sz <- model$Wsat%*%diag(hillslope$area)

    ## initialise the output
    qchannel_output <- reclass( matrix(NA,nrow(obs_data),length(channel$id)), match.to=obs_data )
    names(qchannel_output) <- channel$id

    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){

        ## Step 1: Initialise channel stores with precipitation
        ## channel$store[] <- 0
        channel$store <- obs_data[it,channel$precip_input]*ts$step

        ## Step 2: Distribute any surface storage downslope for next time step
  	if(any(hillslope$ex > 0)){
            ## solve eigen routing
            ex_in[ex_idx] <- hillslope$ex
            ex_out <- eigen_routing_step(ex_in,ex_eigen,ts$step)
            hillslope$ex <- ex_out[ex_idx]
            #browser()
            channel$store <- channel$store + ex_out[-ex_idx]
        }

        ## set the inputs to the hillslope
        precip <- obs_data[it,hillslope$precip_input]
        pet <- obs_data[it,hillslope$pet_input]

        ## inner loop to update the flows and storages
        for(inner in 1:ts$n_sub_step){
            ## Step 3: solve the root zone for hillslope elements
            b <- pmax( pet/hillslope$srz_max , 1e-10 ) # stops round errors in explicit solution
            a <- precip

            tilde_srz = hillslope$srz*exp(-b*ts$sub_step) +
                (a/b)*(1-exp(-b*ts$sub_step)) # analytical solution to root zone storage
            hillslope$srz <- pmin(tilde_srz,hillslope$srz_max)
            integral_qrz <- tilde_srz - hillslope$srz

            ## split root zone flow
            saturated_index <- hillslope$ssz <= 0 # which areas are saturated
            hillslope$ex[saturated_index] <- hillslope$ex[saturated_index] +
                integral_qrz[saturated_index] # in saturated zone goes to surface storage
            hillslope$suz[!saturated_index] <- hillslope$suz[!saturated_index] +
                integral_qrz[!saturated_index] # if not satureated goes to unsaturated zone

            ## Step 4: Unsaturated zone
            ## recharge rate through unsaturated drainage into saturated zone
            ## if(use_cpp){
            ##     hillslope$quz <- funcpp_uz(hillslope$suz, hillslope$ssz, hillslope$td, ts$sub_step)
            ## }else{
            hillslope$quz <- pmin( hillslope$suz / (hillslope$td * hillslope$ssz) , hillslope$suz/ts$sub_step )
            hillslope$quz[ hillslope$ssz==0 ] <- 0
            ## }
            

            ## reduce storage by drainage out of zone over time step - limited in above call
            hillslope$suz <- hillslope$suz - hillslope$quz*ts$sub_step

            ## Step 5: Solve saturated zone
            ## solve the saturated zone to distribute baseflows downslope through areas
            ## using precalculated inter-group splits

            ## Solve the ODE for discharge (ignores limit in flow due to saturation
            # browser()
            res <- deSolve::ode(y=hillslope$lsz,
                                times=seq(0, ts$sub_step, length.out=2),
                                func=fun_dlex_dt,
                                parms=list(Wdash=K_sz,
                                           m=hillslope$m,
                                           lsz_max=hillslope$lsz_max,
                                           quz=hillslope$quz))

            res <- res[-1,-1]; names(res) <- NULL # trim to get only final values and not time
            hillslope$lsz <- res # not due to solution above the equality lsz<=lszmax should hold
            #browser()
            ## work out is there is flow to excess due to limit on lsz
            qsz <- pmax(0, hillslope$quz + K_sz %*% hillslope$lsz) * (hillslope$lsz >= hillslope$lsz_max)
            

            ## compute the gradient of the storage change
            grad_ssz <- -hillslope$quz - K_sz %*% hillslope$lsz + qsz

            ## solve for new storage
            tilde_ssz <- hillslope$ssz + grad_ssz*ts$sub_step
            hillslope_ssz <- pmax(0,tilde_ssz)

            ## add extra to surface excess - include component from saturating storage
            hillslope$ex <- hillslope$ex + (hillslope_ssz - tilde_ssz) + qsz*ts$sub_step

            ## add flow to channel
            channel$store <- channel$store + (FA_sz %*% hillslope$lsz)*ts$sub_step/channel$area

            ## ## get inflows to each hillslope HRU
            ## qin <- vect_matrix_mult_cpp(res*all_area,Wsat) / all_area # distribute downstream and convert back to specific input

            ## ## add specific input to chanel and trim off
            ## channel$store <- channel$store + qin[channel$index]*ts$sub_step
            ## qin <- qin[hillslope$index] # trim to just hillslope inflows

            ## ## excess over capacity generates return flow
            ## ## rate of generation of return flow is excess over maximum of net filling of area
            ## ret_fl <- pmax(qin - hillslope$qsz - hillslope$qsz_max, 0)

            ## ## remove it from the net input as has appeared on the surface!
            ## qin <- qin - ret_fl

            ## ## excess surface storage is generated in this time step.
            ## hillslope$ex <- hillslope$ex + ret_fl*ts$sub_step

            ## ## update storage deficit
            ## hillslope$sd <- hillslope$sd + (hillslope$qsz - qin - hillslope$quz)*ts$sub_step

            ## ## handle excess storage
            ## hillslope$ex <- hillslope$ex - pmin(hillslope$sd , 0)
            ## hillslope$sd <- pmax(hillslope$sd , 0)
        } ## end of subseting loop

        ## copy channel data to output
        ## channel$store is the amount of water accumulated over the timestep
        qchannel_output[it,] <- (channel$store*channel$area)/(ts$step*60*60)

    } ## end of timestep loop

    model$states <- list(hillslope=hillslope,
                         channel=channel)
    
    return( list(model=model,
                 channel_input = qchannel_output ) )
}

#' Run dynamic topmodel
#' @param hru TODO
#' @param param TODO
#' @param Wsat TODO
#' @param Wsurf TODO
#' @param obs_data TODO
#' @param init_discharge initial dischage in m^3/s
#' @param sim_time_step simulation timestep in hours
#' @export
dynatop <- function(model,param,obs_data,init_discharge,sim_time_step){
    ## check inputs are valid and return a 'sim' object
    if(!check_input(model,param,obs_data,init_discharge,sim_time_step)){
        stop("Input check failed")
    }
    ##sim <- model
    ## compute the time step and number of sum timesteps
    ts <- list()
    ts$step <- diff(as.numeric(index(obs_data)[1:2]))/(60*60) # hours
    ts$sub_step <- sim_time_step # hours
    ts$n_sub_step <- floor(ts$step/sim_time_step) # dimensionless

    ## initialise the properties and states of the hillslope hru
    hillslope <- initialise_properties(model,param,'hillslope')
    hillslope <- initialise_hillslope(hillslope,model,init_discharge)

    ## initialise the properties and states of the channel hru
    channel <- initialise_properties(model,param,'channel')
    channel <- initialise_channel(channel)

    ## create the common parts of the surface excess solution
    K_ex <- diag(1/c(hillslope$area,channel$area)) %*%
        cbind( rbind( (model$Wsurf-diag(nrow(model$Wsurf))), model$Fsurf),
              matrix(0,length(hillslope$id)+length(channel$id),length(channel$id) ) # no flow from channel
              ) %*% diag(c(hillslope$area,channel$area)) %*%
        diag(c(1/hillslope$Tex,rep(0,length(channel$id))))
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


    message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){
        ##print(it)

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
            ## solve the root zone for hillslope elements
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

            ## recharge rate through unsaturated drainage into saturated zone
            hillslope$quz <- funcpp_uz(hillslope$suz, hillslope$ssz, hillslope$td, ts$sub_step)

            ## reduce storage by drainage out of zone over time step - limited in above call
            hillslope$suz <- hillslope$suz - hillslope$quz*ts$sub_step

            ## solve the saturated zone to distribute baseflows downslope through areas
            ## using precalculated inter-group splits

            ## Solve the ODE for discharge (ignores limit in flow due to saturation
            # browser()
            res <- deSolve::ode(y=hillslope$lsz,
                                times=seq(0, ts$sub_step, length.out=2),
                                func=funR_dqdt,
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
            channel$store <- channel$store + (FA_sz %*% hillslope$lsz)/channel$area

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

    return( qchannel_output )
}

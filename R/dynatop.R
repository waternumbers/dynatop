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
        check_model(model,check_channel=FALSE,verbose=FALSE)
        tmp <- initialise_dynatop(model,initial_recharge)
        hillslope <- tmp$hillslope
        channel <- tmp$channel
    }

    ## create the common parts of the surface excess solution
    K_ex <- diag(1/c(hillslope$area,channel$area)) %*%
        cbind( rbind( (model$Wex-diag(nrow(model$Wex))), model$Fex),
              matrix(0,length(hillslope$id)+length(channel$id),length(channel$id) ) # no flow from channel
              ) %*% diag(c(hillslope$area,channel$area)) %*%
        diag(c(1/hillslope$tex,rep(0,length(channel$id))))
    ex <- list(eigen = eigen_routing_setup(K_ex),
               input = rep(0,nrow(K_ex)),# preassign vector for initial condition storages
               output = rep(0,nrow(K_ex)),# preassign vector for final condition storages
               idx = 1:length(hillslope$id)# index of hillslope elements in the solution
               )
    rm(K_ex)

    ## initialise the vertical flux stores
    integral_q <- list(ex_rz = rep(0,length(hillslope$id)),
                       rz_ex = rep(0,length(hillslope$id)),
                       rz_uz = rep(0,length(hillslope$id)),
                       uz_sz = rep(0,length(hillslope$id)),
                       sz_ex = rep(0,length(hillslope$id)),
                       )
                       
    ## Common parts for the saturated routing


    ## simple function for solving ODE
    fode <- function(a,b,x0,t){
        b <- pmax(b,1e-10)
        unname( x0*exp(-b*t) + (a/b)*(1-exp(-b*t)) )
    }
    fodeint <- function(a,b,x0,t){
        b <- pmax(b,1e-10)
        a <- pmax(a,1e-10)
        unname( (x0/a)*(1-exp(-b*t)) + (a/b)*t + (a/b^2)*(1-exp(-b*t)) )
    }
    
    
        
    ## initialise the output
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(channel$id)), match.to=obs_data )
    names(channel_inflow) <- channel$id

    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){

        ## set the inputs to the hillslope
        precip <- obs_data[it,hillslope$precip_input]
        pet <- obs_data[it,hillslope$pet_input]

        ## initialise the inflow to be added to as the solution progresses
        channel$inflow <- ( channel$area*obs_data[it,channel$precip_input]*ts$step ) / (3600*ts$step) # rainfall input in cumecs

        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## clear all vertical fluxes
            for(ii in names(integral_q)){
                integral_q[[ii]][] <- 0
            }
            
            ## Step 1: Distribute any surface storage downslope
            if(any(hillslope$ex > 0)){
                ## set input
                ex$input[] <- 0
                ex$input[ex$idx] <- hillslope$ex
                ## solve eigen routing
                ex$output <- eigen_routing_step(ex$input,ex$eigen,ts$sub_step)
                ex$output <- pmax(0,ex$output) ## shouldn't be any negative values but...
                ## evaluate integral of low to rootzone
                integral_q$ex_rz <- pmin( hillslope$qex_max*ts$sub_step,ex$output )
                ## assign storage to hillslope
                tilde_ex <- ex$output - integral_q$ex_rz # this is tilde s_ex
                ## add inflow to channel
                channel$inflow <- channel$inflow +
                    ( (channel$area*ex$output[-ex$idx]) / (3600*ts$sub_step) )
            }

            ## Step 2: solve the root zone for hillslope elements

            ## solve ODE
            tilde_rz <- fode( precip + (integral_q$ex_rz/ts$sub_step),
                              pet/hillslope$srz_max,
                              hillslope$rz,ts$sub_step )
            ## new storage value
            hillslope$rz <- pmin(tilde_rz,hillslope$srz_max)
            ## split root zone flow
            tmp <- tilde_rz - hillslope$rz
            saturated_index <- which(hillslope$sz <= 0) # which areas are saturated
            integral_q$rz_ex[saturated_index] <- tmp[saturated_index]
            integral_q$rz_uz[-saturated_index] <- tmp[-saturated_index]

            ## Step 3: Unsaturated zone

            ## solve ODE
            tilde_uz <- fode( integral_q$rz_uz/ts$sub_step,
                              1 / (hillslope$td * hillslope$sz),
                              hillslope$suz,ts$sub_step )
            integral_q$uz_sz <- hillslope$uz + integral_q$rz_uz - tilde_uz

            ## Step 4: Solve saturated zone
        
            ## initialise the value of Q
            Q <- integral_q$uz_sz

            ## initial estimate of tau, beta and upper limit on alpha
            tau <- hillslope$lsz/hillslope$m
            beta <- tau*(1-diag(model$Wsz))
            ebt <- list(exp(-beta*ts$sub_step), 1-exp(-beta*ts$sub_step))
            amax <- pmax( (hillslope$sz*beta)/ebt[[2]] + (hillslope$lsz*beta/tau),
                         beta*(hillslope$lsz_max - hillslope$lsz*ebt[[1]]) / (tau*ebt[[2]])
                         )

            ## loop HRUs to solve
            tilde_sz <- tilde_lsz <- rep(0,length(hillslope$id))
            Q <- integral_q$uz_sz
            for(ii in 1:length(hillslope$id)){
                ## find alpha
                alpha <- min( amax[ii], Q[ii]/ts$sub_step )
                ## evaluate integral of lateral flux
                tmp <- (hillslope$lsz[ii]/beta[ii])*ebt[[2]][ii] +
                    ( alpha*tau[ii]*ts$sub_step / beta[ii] ) -
                    ( alpha*tau[ii]*ebt[[2]][ii]/(beta[ii]^2) )
                ## add to downslope HRUs
                Q <- Q + model$Wsz[ii,]*tmp
                ## solve for final values
                tilde_lsz[ii] <- hillslope$lsz[ii]*ebt[[1]][ii] + (alpha*tau[ii]/beta[ii])*ebt[[2]][ii]
                tilde_sz[ii] <- hillslope$sz[ii] + ((hillslope$lsz[ii]/tau[ii]) - (alpha/beta[ii]))*ebt[[2]][ii]
            }

            ## revise tau, beta and upper limit on alpha
            tau <- (tau + tilde_lsz/hillslope$m)/2
            beta <- tau*(1-diag(model$Wsz))
            ebt <- list(exp(-beta*ts$sub_step), 1-exp(-beta*ts$sub_step))
            amax <- pmax( (hillslope$sz*beta)/ebt[[2]] + (hillslope$lsz*beta/tau),
                         beta*(hillslope$lsz_max - hillslope$lsz*ebt[[1]]) / (tau*ebt[[2]])
                         )

            ## solve again
            ## loop HRUs to solve
            Q <- integral_q$uz_sz
            for(ii in 1:length(hillslope$id)){
                ## find alpha
                alpha <- min( amax[ii], Q[ii]/ts$sub_step )
                ## evaluate integral of lateral flux
                tmp <- (hillslope$lsz[ii]/beta[ii])*ebt[[2]][ii] +
                    ( alpha*tau[ii]*ts$sub_step / beta[ii] ) -
                    ( alpha*tau[ii]*ebt[[2]][ii]/(beta[ii]^2) )
                ## add to downslope HRUs
                Q <- Q + model$Wsz[ii,]*tmp
                ## add to channel HRUs
                channel$inflow <- channel$inflow + model$F[ii,]*hillslope$area*tmp/(3600*ts$time_step)
                ## solve for final values
                tilde_sz[ii] <- hillslope$sz[ii] + ((hillslope$lsz[ii]/tau[ii]) - (alpha/beta[ii]))*ebt[[2]][ii]
                hillslope$lsz[ii] <- hillslope$lsz[ii]*ebt[[1]][ii] + (alpha*tau[ii]/beta[ii])*ebt[[2]][ii]
            }

            ## step 6: Correct for returned vertical fluxes
            hillslope$ex <- hillslope$ex + integral_q$rz_ex
            saturated_index <- which(hillslope$ssz <= 0) # which areas are saturated
            hillslope$ex[saturated_index] <- hillslope$ex[saturated_index] + hillsope$uz[saturated_index]
            hillsope$uz[saturated_index] <- 0

            ## step 7: Done through the steps above
        } # end of inner time step loop
        
        ## put channel inflow in output matric
        channel_inflow[it,] <- channel$inflow
        
    } ## end of timestep loop

    model$states <- list(hillslope=hillslope,
                         channel=channel)
    
    return( list(model=model,
                 channel_input = channel_inflow ) )
}

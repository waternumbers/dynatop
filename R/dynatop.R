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
dynatop <- function(model,obs_data,initial_recharge=NA,sim_time_step=NULL,use_states=FALSE){

    ## check input and get model timestep
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
    browser()
    
    ## create the common parts of the surface excess solution
    K_ex <- Diagonal(1/c(hillslope$area,channel$area)) %*%
        cbind( rbind( (model$Wex-diag(nrow(model$Wex))), model$Fex),
              matrix(0,length(hillslope$id)+length(channel$id),length(channel$id) ) # no flow from channel
              ) %*% diag(c(hillslope$area,channel$area)) %*%
        diag(c(1/hillslope$tex,rep(0,length(channel$id))))
    ex <- list(expAdt = expm_setup(K_ex,ts$sub_step),
               s_0 = rep(0,nrow(K_ex)),# preassign vector for initial condition storages
               s_dt = rep(0,nrow(K_ex)),# preassign vector for final condition storages
               idx = 1:length(hillslope$id)# index of hillslope elements in the solution
               )
    rm(K_ex)
    
    ## initialise the vertical flux stores
    q_vol <- list(ex_rz = rep(0,length(hillslope$id)),
                  rz_ex = rep(0,length(hillslope$id)),
                  rz_uz = rep(0,length(hillslope$id)),
                  uz_sz = rep(0,length(hillslope$id)),
                  sz_ex = rep(0,length(hillslope$id)),
                  sz = rep(0,length(hillslope$id)) # integral of l_ex
                  )
    
    ## Common parts for the saturated routing
    Dsz <- Diagonal(1/hillslope$area) %*% (Diagonal(ncol(model$Wsz)) - model$Wsz) %*% Diagonal(hillslope$area)

    Esz <- model$Fsz %*% Diagonal(hillslope$area) 
    
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
        
        ## set the inputs to the hillslope
        precip <- obs_data[it,hillslope$precip_input]
        pet <- obs_data[it,hillslope$pet_input]
        
        ## set accumulated inflow to channel to be 0
        channel$ves <- channel$vsz <- channel$vp <- rep(0,length(channel$id))
        # volume from precipitation
        channel$vp <- channel$area * obs_data[it,channel$precip_input] * ts$step # rainfall input in m^3
        
        
        ## loop sub steps
        for(inner in 1:ts$n_sub_step){
            
            ## clear all vertical fluxes
            for(ii in names(integral_q)){
                integral_q[[ii]][] <- 0
            }
            
            ## Step 1: Distribute any surface storage downslope
            if(any(hillslope$ex > 0)){
                ## set input
                ex$s_0[] <- 0
                ex$s_0[ex$idx] <- hillslope$ex
                ## solve eigen routing
                ex$s_dt <- pmax( ex$expAdt%*%ex$s_0 ,0 ) ## shouldn't be any negative values but...
                ## evaluate integral of flow to rootzone
                integral_q$ex_rz <- pmin( hillslope$qex_max*ts$sub_step,ex$s_dt[ex$idx] )
                ## inflow to channel
                channel$vex <- channel$vex + channel$area*ex$s_dt[-ex$idx]

                ## set hillslope$ex - corrected in step 5
                hillslope$ex <- ex$s_dt[ex$idx]
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
            sz0 <- hillslope$sz
            lsz0 <- asdf #initial velocity

            ## inital gradient
            gsz0 <- as.numeric(Dsz%*%lsz0) - (int_q$uz_sz/ts$sub_step)
            gsz0[saturated_index] <- pmax(gsz0[saturated_index],0)

            ## intermediate estimate and gradient
            szp <- pmax(sz0 + ts$sub_step*gsz0,0)
            lszp <- asdf
            gszp <- as.numeric(Dsz%*%lszp) - (int_q$uz_sz/ts$sub_step)
            tmp <- which(szp <= 0)
            gszp[tmp] <- pmax(gszp[tmp],0)

            ## final estimate
            hillslope$sz <- pmax(sz0 + ts$sub_step*(gsz0+gszp)/2,0)
            lsz <- asdf
            int_q$sz <- ts$sub_step*(lsz0+lsz)/2
            int_q$sz_ex <- pmax(hillslope$sz - sz0 - Dsz%*%int_q$sz + int_q$uz_sz,0)
            ## step 5 - correct the stores
            saturated_index <- which(hillslope$sz <= 0)
            hillslope$sz <- hillslope$sz - int_q$ex_rz + int_q$sz_ex
            hillslope$sz[saturated_index] <- hillslope$sz[saturated_index] +
                hillslope$uz[saturated_index]
            hillslope$uz[saturated_index] <- 0

            ## step 6 - channel inflow
            channel$vsz <- channel$vsz + Esz %*% int_q$sz
        }
        
                         
        ##     browser()
            
        ##     ## initial estimate of tau, beta and upper limit on alpha
        ##     tau <- hillslope$lsz/hillslope$m
        ##     beta <- tau*(1-diag(model$Wsz))
        ##     ebt <- list(exp(-beta*ts$sub_step), 1-exp(-beta*ts$sub_step))
        ##     amax <- pmax( (hillslope$sz*beta)/ebt[[2]] + (hillslope$lsz*beta/tau),
        ##                  beta*(hillslope$lsz_max - hillslope$lsz*ebt[[1]]) / (tau*ebt[[2]])
        ##                  )
            
        ##     ## loop HRUs to solve
        ##     alpha <- rep(0,length(hillslope$id))
        ##     Q <- integral_q$uz_sz
        ##     for(ii in 1:length(hillslope$id)){
        ##         ## find alpha
        ##         alpha[ii] <- min( amax[ii], Q[ii]/ts$sub_step )
        ##         ## evaluate integral of lateral flux
        ##         tmp <- (hillslope$lsz[ii]/beta[ii])*ebt[[2]][ii] +
        ##             ( alpha*tau[ii]*ts$sub_step / beta[ii] ) -
        ##             ( alpha*tau[ii]*ebt[[2]][ii]/(beta[ii]^2) )
        ##         ## add to downslope HRUs
        ##         Q <- Q + (model$Wsz[ii,]/hillslope$area)*hillslope$area[ii]*tmp
        ##     }
            
        ##     ## solve for intermediate values
        ##     tilde_lsz <- hillslope$lsz*ebt[[1]] + (alpha*tau/beta)*ebt[[2]]
        ##     tilde_sz <- hillslope$sz + ((hillslope$lsz/tau) - (alpha/beta))*ebt[[2]]
            
        ##     ## revise tau, beta and upper limit on alpha
        ##     tau <- (tau + tilde_lsz/hillslope$m)/2
        ##     beta <- tau*(1-diag(model$Wsz))
        ##     ebt <- list(exp(-beta*ts$sub_step), 1-exp(-beta*ts$sub_step))
        ##     amax <- pmax( (hillslope$sz*beta)/ebt[[2]] + (hillslope$lsz*beta/tau),
        ##                  beta*(hillslope$lsz_max - hillslope$lsz*ebt[[1]]) / (tau*ebt[[2]])
        ##                  )
            
        ##     ## solve again
        ##     ## loop HRUs to solve
        ##     alpha <- rep(0,length(hillslope$id))
        ##     Q <- integral_q$uz_sz
        ##     for(ii in 1:length(hillslope$id)){
        ##         ## find alpha
        ##         alpha <- min( amax[ii], Q[ii]/ts$sub_step )
        ##         ## evaluate integral of lateral flux
        ##         tmp <- (hillslope$lsz[ii]/beta[ii])*ebt[[2]][ii] +
        ##             ( alpha*tau[ii]*ts$sub_step / beta[ii] ) -
        ##             ( alpha*tau[ii]*ebt[[2]][ii]/(beta[ii]^2) )
        ##         ## add to downslope HRUs
        ##         Q <- Q + (model$Wsz[ii,]/hillslope$area)*hillslope$area[ii]*tmp
        ##         ## add to channel HRUs
        ##         channel$vsz <- channel$vsz + model$F[ii,]*hillslope$area*tmp
        ##     }
            
        ##     ## solve for final values
        ##     hillslope$sz <- hillslope$sz + ((hillslope$lsz/tau) - (alpha/beta))*ebt[[2]]
        ##     hillslope$lsz <- hillslope$lsz*ebt[[1]] + (alpha*tau/beta)*ebt[[2]]
        ##     integral_q$sz_rz <- Q - alpha*ts$sub_step
            
        ##     ## step 6: Correct for returned vertical fluxes
        ##     hillslope$ex <- ex$s_dt[ex$idx] - integral_q$ex_rz + integral_q$rz_ex + integral_q$sz_ex 
        ##     saturated_index <- which(hillslope$ssz <= 0) # which areas are saturated
        ##     hillslope$ex[saturated_index] <- hillslope$ex[saturated_index] + hillslope$uz[saturated_index]
        ##     hillslope$uz[saturated_index] <- 0
            
        ## } # end of inner time step loop
        
        ## step 7: Done through the steps above
        ## put channel inflow in output matric
        channel_inflow[it,] <- (channel$vex + channel$vsz +  channel$vp) / (3600*ts$step)
        
    } ## end of timestep loop
    
    model$states <- list(hillslope=hillslope,
                         channel=channel)
    
    return( list(model=model,
                 channel_input = channel_inflow ) )
}

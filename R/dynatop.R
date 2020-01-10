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
        hillslope <- model$states$hillslope
        channel <- model$states$channel
    }else{ # initialise the model
        tmp <- initialise_dynatop(model,initial_recharge)
        hillslope <- tmp$hillslope
        channel <- tmp$channel
    }
    
    ## create the common parts of the surface excess solution
    # TO DO - fix this!
    ## tmp_area <- rep(NA,ncol(mdl$Dex))
    ## tmp_area[hillslope$id] <- hillslope$area
    ## tmp_area[channel$id] <- channel$area
    ## tmp_inv_tex <- rep(0,ncol(mdl$Dex))
    ## tmp_inv_tex[hillslope$id] <- 1/hillslope$tex
    ## K_ex <- Diagonal(length(tmp_area),1/tmp_area) %*% mdl$Dex %*%
    ##     Diagonal(length(tmp_area),tmp_area*tmp_inv_tex)
    ## ## TO DO fix matrix exp
    ## ex <- list(expAdt = K_ex, #expm_setup(K_ex,ts$sub_step),
    ##            s_0 = rep(0,nrow(K_ex)),# preassign vector for initial condition storages
    ##            s_dt = rep(0,nrow(K_ex)) # pre assign vector for result
    ##            )
    ## rm(K_ex,tmp_inv_tex,tmp_area)
    
    ## initialise the vertical flux stores
    q_vol <- list(ex_rz = rep(0,length(hillslope$id)),
                  rz_ex = rep(0,length(hillslope$id)),
                  rz_uz = rep(0,length(hillslope$id)),
                  uz_sz = rep(0,length(hillslope$id)),
                  sz_ex = rep(0,length(hillslope$id)),
                  sz_in = rep(0,length(hillslope$id)), # integral of l_ex inflow
                  sz = rep(0,length(hillslope$id)) # integral of l_ex outflow
                  )
    
    ## Common parts for the saturated routing
    Wsz <- mdl$Dsz[hillslope$id,hillslope$id]
    Fsz <- mdl$Dsz[channel$id,hillslope$id]
    A <- Diagonal(length(hillslope$area),hillslope$area)
    invA <- Diagonal(length(hillslope$area),1/hillslope$area)
    invAWA = invA %*% Wsz %*% A
    sz <- list(band=list(),
               FszA = Fsz %*% A)
    for(ii in sort(unique(hillslope$band))){
        sz$band[[ii]] <- list()
        idx <- which(hillslope$band==ii)
        X <- invAWA[,idx,drop=FALSE]
        jdx <- which(rowSums(X)>0)
        
        sz$band[[ii]]$idx <- idx
        sz$band[[ii]]$jdx <- jdx
        sz$band[[ii]]$X <- X[jdx,]
    }
    
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
        channel$vex <- channel$vsz <- channel$vp <- rep(0,length(channel$id))
        # volume from precipitation
        channel$vp <- unname( channel$area * obs_data[it,channel$precip_input] * ts$step) # rainfall input in m^3
        
        
        ## loop sub steps
        for(inner in 1:ts$n_sub_step){
            #print(paste(it,inner))
            ## clear all vertical fluxes
            for(ii in names(q_vol)){
                q_vol[[ii]][] <- 0
            }

            
            ## Step 1: Distribute any surface storage downslope
            ## TO DO reinstate
            if(any(hillslope$ex > Inf)){
                ## set input
                ex$s_0[] <- 0
                ex$s_0[hillslope$id] <- hillslope$ex
                ## solve eigen routing
                ex$s_dt <- pmax( as.numeric(ex$expAdt%*%ex$s_0) ,0 ) ## shouldn't be any negative values but...
                ## revise hillslope stroage
                hillslope$ex <- ex$s_dt[hillslope$id] # this is tilde in document
                ## inflow to channel
                channel$vex <- channel$vex + channel$area*ex$s_dt[channel$id]    
            }
            ## evaluate integral of flow to rootzone
            q_vol$ex_rz <- pmin( hillslope$qex_max*ts$sub_step,hillslope$ex )

            ## Step 2: solve the root zone for hillslope elements
            
            ## solve ODE
            tilde_rz <- fode( precip + (q_vol$ex_rz/ts$sub_step),
                             pet/hillslope$srz_max,
                             hillslope$rz,ts$sub_step )
            ## new storage value
            hillslope$rz <- pmin(tilde_rz,hillslope$srz_max)
            ## split root zone flow
            tmp <- tilde_rz - hillslope$rz
            saturated_index <- hillslope$sz <= 0 # which areas are saturated
            q_vol$rz_ex <- tmp * saturated_index
            q_vol$rz_uz <- tmp * !saturated_index
            
            ## Step 3: Unsaturated zone
            
            ## solve ODE
            tilde_uz <- fode( q_vol$rz_uz/ts$sub_step,
                             1 / (hillslope$td * hillslope$sz),
                             hillslope$uz,ts$sub_step )
            q_vol$uz_sz <- hillslope$uz + q_vol$rz_uz - tilde_uz
            hillslope$uz <- tilde_uz

            
            ## Step 4: Solve saturated zone

            ## compute initial cbar based on previous flows
            Qbar <- pmin(0.5*(hillslope$lsz + hillslope$lsz_in),hillslope$lsz_max)
            c_bar <- (Qbar*hillslope$delta_x)/hillslope$m

            ## compute parameters
            lambda <- sz_opt$omega + sz_opt$theta*c_bar*ts$sub_step/hillslope$delta_x
            lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*c_bar*ts$sub_step/hillslope$delta_x
            ## compute constant
            k <- lambda_prime*hillslope$lsz +
                (1-lambda_prime)*sz$qin +
                c_bar*q_vol$uz_sz/hillslope$delta_x
            
            ## pass to compute estimate
            k <- k/lambda
            lambda_prime <- lambda_prime/lambda
            browser()
            sz$lsz <- hillslope$lsz # old value for computing integral
            sz$lsz_in <- hillslope$lsz_in
            hillslope$lsz_in[] <- 0
            hillslope$lsz[] <- 0
            for(ii in sort(unique(hillslope$band))){
                tmp <- sz$band[[ii]]
                lsz <- k[tmp$idx] - lambda_prime[tmp$idx]*hillslope$lsz_in[tmp$idx]
                lsz <- pmax(lsz,hillslope$lsz_max[tmp$idx])
                hillslope$lsz[tmp$idx] <- lsz
                hillslope$lsz_in[tmp$jdx] <- hillslope$lsz_in[tmp$jdx] + tmp$X%*%lsz
            }
            
                
            ##     hillslope$lsz_in <- hillslope$lsz_in
            ## for(bnd in sort(unique(hillslope$band))){
            ##     idx <- which(hillslope$band==bnd)
            ##     ldz <- hillslope$lsz[idx]
            ##     k <- k[idx]
            ##     lsz_in <- hillslope$lsz_in[idx]
            ##     lp <- lambda_prime[idx]
            ##     X <- sz$invAWA[,idx,drop=FALSE]
            ##     lsz_max <- hillslope$lsz_max[idx]
            
            ##     lsz <- k - lp*lsz_in
            ##     lsz <- pmax(lsz,lsz_max)

            ##     y <- X%*%lsz
            ##     hillslope$lsz_in <- hillslope$lsz_in + y #X%*%lsz
            ##     ##hillslope$lsz[idx] <- k[idx] - lambda_prime[idx]*hillslope$lsz_in[idx]
            ##     ##hillslope$lsz[idx] <- pmax(hillslope$lsz[idx],hillslope$lsz_max[idx])
            ##     ##X <- sz$invAWA[,idx,drop=FALSE]
            ##     ##hillslope$lsz_in <- hillslope$lsz_in + X%*%hillslope$lsz[idx]
            ## }

            
            ## bnd <- 1
            ## for(ii in 1:length(hillslope$id)){
            ##     ## if we have moved down a band recompute the fluxes
            ##     #browser()
            ##     if(hillslope$band[ii]!=bnd){
            ##         hillslope$lsz_in <- sz$invAWA%*%hillslope$lsz
            ##         bnd <- hillslope$band[ii]
            ##         #print(bnd)
            ##     }
                
            ##     ## compute new flow
            ##     hillslope$lsz[ii] <- k[ii] - lambda_prime[ii]*hillslope$lsz_in[ii]
            ##     hillslope$lsz[ii] <- max(hillslope$lsz[ii],hillslope$lsz_max[ii])
             
                
            ## }

            ## work out integral fluz
            #browser()
            q_vol$sz_in <- ts$sub_step*( hillslope$lsz_in + sz$lsz_in )/2
            q_vol$sz <- ts$sub_step*( hillslope$lsz + sz$lsz )/2
            tilde_sz <- hillslope$sz + q_vol$sz - q_vol$sz_in - q_vol$uz_sz
            hillslope$sz <- pmax(0,tilde_sz)
            q_vol$sz_ez <- hillslope$sz - tilde_sz
            
            ## step 5 - correct the stores
            saturated_index <- hillslope$sz <= 0
            hillslope$ex <- hillslope$ex + q_vol$rz_ex +
                (q_vol$sz_ex + hillslope$uz)*saturated_index
            hillslope$uz <- hillslope$uz * !saturated_index

            ## step 6 - channel inflow
            channel$vsz <- channel$vsz + as.numeric(sz$FszA %*% q_vol$sz)
        }
        
        ## step 7: Done through the steps above
        ## put channel inflow in output matric
        channel_inflow[it,] <- (channel$vex + channel$vsz +  channel$vp) / (3600*ts$step)
        
    } ## end of timestep loop
    
    model$states <- list(hillslope=hillslope,
                         channel=channel)
    
    return( list(model=model,
                 channel_input = channel_inflow ) )
}

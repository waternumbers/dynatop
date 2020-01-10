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
dynatop_con <- function(model,obs_data,
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
    tmp_area <- rep(NA,ncol(mdl$Dex))
    tmp_area[hillslope$id] <- hillslope$area
    tmp_area[channel$id] <- channel$area
    tmp_inv_tex <- rep(0,ncol(mdl$Dex))
    tmp_inv_tex[hillslope$id] <- 1/hillslope$tex
    K_ex <- Diagonal(length(tmp_area),1/tmp_area) %*% mdl$Dex %*%
        Diagonal(length(tmp_area),tmp_area*tmp_inv_tex)
    ex <- list(expAdt = expm_setup(K_ex,ts$sub_step),
               s_0 = rep(0,nrow(K_ex)),# preassign vector for initial condition storages
               s_dt = rep(0,nrow(K_ex)) # pre assign vector for result
               )
    rm(K_ex,tmp_inv_tex,tmp_area)
    
    ## initialise the vertical flux stores
    q_vol <- list(ex_rz = rep(0,length(hillslope$id)),
                  rz_ex = rep(0,length(hillslope$id)),
                  rz_uz = rep(0,length(hillslope$id)),
                  uz_sz = rep(0,length(hillslope$id)),
                  sz_ex = rep(0,length(hillslope$id)),
                  sz = rep(0,length(hillslope$id)) # integral of l_ex
                  )
    
    ## Common parts for the saturated routing
    sz <- list(Wsz = mdl$Dsz[hillslope$id,hillslope$id],
               Fsz = mdl$Dsz[channel$id,hillslope$id],
               Diag_lambda = Diagonal(length(hillslope$id)),
               Diag_mlambda = Diagonal(length(hillslope$id)),
               X = Diagonal(length(hillslope$id))
               )
    sz$invAWA = Diagonal(length(hillslope$area),1/hillslope$area) %*% sz$Wsz %*% Diagonal(length(hillslope$area),hillslope$area)
    sz$FszA = sz$Fsz %*% Diagonal(length(hillslope$area),hillslope$area)
    
    
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
            if(any(hillslope$ex > 0)){
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

            ## compute initial cbar
            Qin <- as.numeric( sz$invAWA %*% pmin(hillslope$lsz,hillslope$lsz_max) )
            Qinn <- pmax(Qin,hillslope$lsz_max)
            Qbar <- pmin(0.5*(hillslope$lsz+Qinn),hillslope$lsz_max)
            c_bar <- (Qbar*hillslope$delta_x)/hillslope$m

            ## estimate \check{lsz}
            lambda <- sz_opt$omega + sz_opt$theta*c_bar*ts$sub_step/hillslope$delta_x
            lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*c_bar*ts$sub_step/hillslope$delta_x
            diag(sz$Diag_lambda) <- lambda
            ##Diag_lambda <- Diagonal(length(lambda),lambda)
            k <- lambda_prime*hillslope$lsz +
                (1-lambda_prime)*Qinn +
                c_bar*q_vol$uz_sz/hillslope$delta_x
            
            ##X <- Diagonal(length(lambda),lambda) + Diagonal(length(lambda),1-lambda)%*%sz$invAWA
            ##X <- Diag_lambda + sz$invAWA - Diag_lambda%*%sz$invAWA
            diag(sz$Diag_mlambda) <- 1-lambda
            ##sz$X <- sz$Diag_lambda + sz$Diag_mlambda%*%sz$invAWA
            sz$X <- sz$Diag_mlambda%*%sz$invAWA
            diag(sz$X) <- diag(sz$X) + lambda
            
            #browser()
            lsz <- as.numeric( solve(sz$X,k) )

            ## ## recompute cbar and values for initial condition
            ## Qbar <- pmin( (hillslope$lsz + Qin + lsz +
            ##                as.numeric(sz$invAWA %*% pmin(lsz,hillslope$lsz_max)))/4 ,
            ##              hillslope$lsz_max)
            ## c_bar <- (Qbar*hillslope$delta_x)/hillslope$m
            ## lambda <- sz_opt$omega + sz_opt$theta*c_bar*ts$sub_step/hillslope$delta_x
            ## lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*c_bar*ts$sub_step/hillslope$delta_x
            ## k <- lambda_prime*hillslope$lsz +
            ##     (1-lambda_prime)*Qin +
            ##     c_bar*q_vol$uz_sz/hillslope$delta_x
            ## X <- Diagonal(length(lambda),lambda) + Diagonal(length(lambda),1-lambda)%*%sz$invAWA
            ## lsz <- as.numeric( solve(X,k) )

            #lsz <- pmax(lsz,0)
            
            ## start iterating
            #print(range(hillslope$sz))
            #print(range(q_vol$uz_sz))
            #print(range(hillslope$lsz))
            #print(range(Qbar))
            #print(range(c_bar))
            #X <- Diagonal(length(lambda),1-lambda) %*% sz$invAWA
            #X <- sz$invAWA - Diag_lambda %*% sz$invAWA
            sz$X <- sz$Diag_mlambda %*% sz$invAWA
            flg <- TRUE
            cnt <- 0
            while(flg){
                Qinn <- pmax(as.numeric(sz$X %*% pmin(lsz,hillslope$lsz_max)),hillslope$lsz_max)
                lsz_new <-  (k - Qinn)/lambda
                
                ##lsz_new <-  (k - as.numeric(X %*% pmin(lsz,hillslope$lsz_max)))/lambda
                #print(range(lsz_new))
                cnt <- cnt + 1
                tol <- max(abs(lsz_new-lsz))
                if(tol < sz_opt$tol){
                    flg <- FALSE
                }
                if(cnt > sz_opt$max_iter){
                    warning(paste("Maximum number of iterations exceeded in saturated zone solution on step",it,". Current tol:",tol))
                    flg <- FALSE
                }
                lsz <- lsz_new
                #lsz <- pmax(lsz_new,0)
            }
            ##print(paste(cnt,tol))
            hillslope$lsz <- pmin(lsz,hillslope$lsz_max) 
            hillslope$lsz <- lsz
            #browser()
            ## work out integral fluz
            q_vol$sz <- ts$sub_step*(pmin(lsz,hillslope$lsz_max) +
                         pmin(hillslope$lsz,hillslope$lsz_max))/2
            tilde_sz <- hillslope$sz + q_vol$sz - as.numeric(sz$invAWA%*%q_vol$sz) - q_vol$uz_sz
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

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
        hillslope <- model$states$hillslope
        channel <- model$states$channel
        ex <- model$states$ex
        sz <- model$states$sz
    }else{ # initialise the model
        hillslope <- initialise_hillslope(model,initial_recharge)
        channel <- initialise_channel(model)

        ## create the common parts of the surface excess solution
        tmp <- list(area = rep(NA,ncol(mdl$Dex)),
                    tex = rep(0,ncol(mdl$Dex)))
        tmp$area[c(hillslope$id,channel$id)] <- c(hillslope$area,channel$area)
        tmp$tex[hillslope$id] <- hillslope$tex
        tmp$invAWA <- Diagonal(length(tmp$area),1/tmp$area) %*% mdl$Dex %*%
            Diagonal(length(tmp$area),tmp$area)
        tmp$B1 <- Diagonal(length(tmp$area),exp(-ts$sub_step/tmp$tex))
        tmp$B2 <- Diagonal(length(tmp$area),(1-exp(-ts$sub_step/tmp$tex))*(tmp$tex/ts$sub_step))
        tmp$I <- Diagonal(length(tmp$area))
        ex <- list()
        ex$invAWA_I = tmp$invAWA - tmp$I
        ## TO DO find a solution methiod
        ex$D <- ex$invAWA_I #solve(ex$invAWA_I - tmp$B2%*%tmp$invAWA , tmp$B1-tmp$I)

        ex$K <- tmp$I + (ex$invAWA_I %*% ex$D)
        ex$dex <- rep(NA,ncol(model$Dex)) # depth in store
        rm(tmp)

        ## create the common parts of the saturated zone routing
        #browser()
        tmp <- list()
        tmp$Wsz <- as.matrix(summary(model$Dsz[hillslope$id,hillslope$id]))
        tmp$Wsplit <- by(tmp$Wsz[,2:3],tmp$Wsz[,1],unique,simplify=FALSE)
        ## tmp$Wsplit <- lapply(tmp$Wsplit,function(x){list(j=as.vector(x[,1]),
        ##                                                  x=as.vector(x[,2]))})
        tmp$Fsz <- as.matrix(summary(model$Dsz[channel$id,hillslope$id]))
        tmp$Fsplit <- by(tmp$Fsz[,2:3],tmp$Fsz[,1],unique,simplify=FALSE)
        sz <- list(hs_parent=rep(list(NULL),length(hillslope$id)),
                   ch_parent=rep(list(NULL),length(channel$id)))
        sz$hs_parent[unique(tmp$Wsz[,1])] <- tmp$Wsplit
        sz$ch_parent[unique(tmp$Fsz[,1])] <- tmp$Fsplit
        sz$order <- order(hillslope$band)
        sz$lsz <- rep(NA,length(hillslope$id))
        sz$lsz_in <- rep(NA,length(hillslope$id))
        rm(tmp)
    }


    ## initialise the stores of integral flux
    q_vol <- list(
        ex_rz = rep(0,length(hillslope$id)),
        rz_ex = rep(0,length(hillslope$id)),
        rz_uz = rep(0,length(hillslope$id)),
        uz_sz = rep(0,length(hillslope$id)),
        sz_ex = rep(0,length(hillslope$id)),
        sz = rep(0,length(hillslope$id)) # integral of l_ex
    )

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
        ## set the inputs to the hillslope
        precip <- obs_data[it,hillslope$precip_series]
        pet <- obs_data[it,hillslope$pet_series]


        ## loop sub steps
        for(inner in 1:ts$n_sub_step){
            #print(paste(it,inner))
            ## clear all vertical fluxes
            for(ii in names(q_vol)){
                q_vol[[ii]][] <- 0
            }


            ## Step 1: Distribute any surface storage downslope
            if(any(hillslope$ex > 0)){
                print("ex zone")
                ## set input
                ex$dex[] <- 0
                ex$dex[hillslope$id] <- hillslope$ex
                ## evolve dex
                ex$dex <- pmax(0,ex$K %*% ex$dex)
                ## put back into hillslope
                hillslope$ex <- ex$dex[hillslope$id]
            }

            ## evaluate max integral of flow to rootzone
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
            sz$lsz[] <- 0
            sz$lsz_in[] <- 0
            ## K <- ex$invAWA_I + Diagonal(ncol(ex$invAWA))
            ## K <- K[hillslope$id,hillslope$id]
            ## for(ord in sort(unique(hillslope$band))){
            ##     #print(ord)
            ##     #browser()
            ##     idx <- hillslope$band==ord

            ##     KK <- K[idx,]

            ##     sz$lsz_in[idx] <- KK %*% sz$lsz # latest inflow
            ##     lsz_in_sat <- pmax(sz$lsz_in[idx],hillslope$lsz_max[idx]) # constrained to saturated zone max flow
            ##     lsz_in_sat_prev <- pmax(hillslope$lsz_in[idx],hillslope$lsz_max[idx]) # constrained to saturated zone max flow

            ##     qbar <- (hillslope$lsz[idx]+lsz_in_sat_prev+lsz_in_sat + pmax(lsz_in_sat,hillslope$lsz[idx]))/4

            ##     # TO DO check the use or not of area in the following
            ##     cbar <- (qbar*hillslope$delta_x[idx])/hillslope$m[idx]
            ##     lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/hillslope$delta_x[idx]
            ##     lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/hillslope$delta_x[idx]

            ##     k <- lambda_prime*hillslope$lsz[idx] +
            ##         (1-lambda_prime)*lsz_in_sat_prev +
            ##         cbar*q_vol$uz_sz[idx]/hillslope$delta_x[idx]
            ##     sz$lsz[idx] <- min( (k - (1-lambda)*lsz_in_sat)/lambda , hillslope$lsz_max[idx] )
            ## }

            lsz_in_sat_prev <- pmin(hillslope$lsz_in,hillslope$lsz_max) # previsou inflow limite to saturation
            for(ii in sz$order){
                ## inflow from upstream

                #idx <- sz$hs_parent[[ii]]$j #[,1]
                #w <- sz$hs_parent[[ii]]$x #[,2]

                #sz$lsz_in[ii] <- sum(w * sz$lsz[idx]) # latest inflow
                sz$lsz_in[ii] <- sum(sz$hs_parent[[ii]]$x * sz$lsz[sz$hs_parent[[ii]]$j]) # latest inflow

                lsz_in_sat <- max(sz$lsz_in[ii],hillslope$lsz_max[ii]) # constrained to saturated zone max flow

                qbar <- (hillslope$lsz[ii]+lsz_in_sat_prev[ii]+lsz_in_sat + max(lsz_in_sat,hillslope$lsz[ii]))/4

                # TO DO check the use or not of area in the following
                cbar <- (qbar*hillslope$delta_x[ii])/hillslope$m[ii]
                lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/hillslope$delta_x
                lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/hillslope$delta_x

                # not ii enough

                k <- lambda_prime*hillslope$lsz[ii] +
                    (1-lambda_prime)*lsz_in_sat_prev[ii] +
                    cbar*q_vol$uz_sz[ii]/hillslope$delta_x[ii]
                sz$lsz[ii] <- min( (k - (1-lambda)*lsz_in_sat)/lambda , hillslope$lsz_max[ii] )
            }

            #browser()
            ## work out integral fluz
            q_vol$sz <- ts$sub_step*(sz$lsz + hillslope$lsz)/2
            tilde_sz <- hillslope$sz + q_vol$sz - ts$sub_step*(sz$lsz_in + hillslope$lsz_in)/2 - q_vol$uz_sz
            hillslope$sz <- pmax(0,tilde_sz)
            q_vol$sz_ex <- hillslope$sz - tilde_sz

            ## update states
            hillslope$lsz <- sz$lsz
            hillslope$lsz_in <- sz$lsz_in

            ## step 5 - correct the stores
            saturated_index <- hillslope$sz <= 0
            hillslope$ex <- hillslope$ex + q_vol$rz_ex +
                (q_vol$sz_ex + hillslope$uz)*saturated_index
            hillslope$uz <- hillslope$uz * !saturated_index

            ## step 6 - channel inflow - at the moment a volume / area
            channel_inflow[it,] <- channel_inflow[it,] + ex$dex[channel$id] + obs_data[it,channel$precip_series]*ts$sub_step
            for(ii in 1:length(channel$id)){
                idx <- sz$ch_parent[[ii]][,1]
                w <- sz$ch_parent[[ii]][,2]
                channel_inflow[it,ii]  <-  channel_inflow[it,ii] +sum( w*hillslope$area[idx]*q_vol$sz[idx]/channel$area[ii] )
            }


        }

        ## step 7: Done through the steps above
        ## put channel inflow in output matric
        channel_inflow[it,] <- channel_inflow[it,]*channel$area / (3600*ts$step)

    } ## end of timestep loop

    model$states <- list(hillslope=hillslope,
                         channel=channel,
                         ex=ex,
                         sz=sz)

    return( list(model=model,
                 channel_input = channel_inflow ) )
}

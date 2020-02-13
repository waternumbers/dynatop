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
dynatop_hsu <- function(model,obs_data,
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
        hsu <- model$hsu
        sq <- model$sq
    }else{ # initialise the model
        tmp <- initialise_hsu(model,initial_recharge)
        hsu <- tmp$hsu
        sq <- tmp$sq
        rm(tmp)
    }

    print("general stores")

    ## some code stores to make passing easier - contain area time lateral flux
    alsz <- vapply(hsu,function(x){x$prop['area']*x$state['lsz']},numeric(1))
    alsf <- vapply(hsu,function(x){x$prop['area']*x$state['lsf']},numeric(1))
    
    ## simple function for solving ODE
    fode <- function(a,b,x0,t){
        b[b<1e-10] <- 1e-10
        #if (b < 1e-10) { b <- 1e-10 }
        #b <- pmax(b,1e-10)
        ##unname( x0*exp(-b*t) + (a/b)*(1-exp(-b*t)) )
        x0*exp(-b*t) + (a/b)*(1-exp(-b*t))
    }

    ## initialise the output
    channel_id <- which( vapply(hsu,function(x){x$type=="channel"},logical(1)) )
    channel_inflow <- reclass( matrix(NA,nrow(obs_data),length(channel_id)), match.to=obs_data )
    names(channel_inflow) <- channel_id

    ## remove time properties from input - stops variable type errors later
    obs_data <- as.matrix(obs_data)


    ## message("Running Dynamic TOPMODEL using ", length(hillslope$id), " hillslope units and ", length(channel$id), " channel units")

    for(it in 1:nrow(obs_data)){
        #browser()
        print(it)
        ## set the inputs to the hillslope
        for(ii in 1:length(hsu)){
            hsu[[ii]]$input['precip'] <- obs_data[it,hsu[[ii]]$series['precip']]/ts$step
            hsu[[ii]]$input['pet'] <- obs_data[it,hsu[[ii]]$series['pet']]/ts$step
        }
        
        ## loop sub steps
        for(inner in 1:ts$n_sub_step){

            ## PASS 1: pass in order of surface redistribution
            print("pass 1")
            for(ii in sq$sf){
                ##browser()
                ##print(ii)
                ## compute the inflow
                hsu[[ii]]$input['lsf'] <- unname( sum( hsu[[ii]]$parent$sf$x * alsf[hsu[[ii]]$parent$sf$j] ) / hsu[[ii]]$prop['area'] )
                if( hsu[[ii]]$type == "hillslope" ){
                    ## solve for surface water flux
                    tilde_sf <- fode( hsu[[ii]]$input['lsf'] / ts$sub_step,
                                     1/hsu[[ii]]$param['tex'],
                                     hsu[[ii]]$state['ssf'],
                                     ts$sub_step)
                    hsu[[ii]]$output['lsf'] <- hsu[[ii]]$state['ssf'] +
                        hsu[[ii]]$input['lsf'] - tilde_sf
                    hsu[[ii]]$state['ssf'] <- tilde_sf
                    alsf[ii] <- hsu[[ii]]$output['lsf']*hsu[[ii]]$prop['area']

                    ## solve for root zone
                    hsu[[ii]]$flux['qsf_rz'] <- min(
                        hsu[[ii]]$state['ssf'],
                        hsu[[ii]]$param['qsf_max']*ts$sub_step
                    )

                    a <- hsu[[ii]]$input['precip'] +
                        hsu[[ii]]$flux['qsf_rz'] / ts$sub_step
                    b <- hsu[[ii]]$input['pet']/hsu[[ii]]$param['srz_max']
                    if( b== 0 ){
                        tilde_rz <- hsu[[ii]]$state['srz']
                    }else{
                        ebt <- exp(-b*ts$sub_step)
                        tilde_rz <- hsu[[ii]]$state['srz']*ebt + (a/b)*(1-ebt)
                    }
                    
                    #b[b<1e-10] <- 1e-10
                    #ebt <- exp(-b*ts$sub_step)
                    #x0*exp(-b*t) + (a/b)*(1-exp(-b*t))
                    #tilde_rz <- hsu[[ii]]$state['srz']*ebt + (a/b)*(1-ebt)
                    
                    ## tilde_rz <- fode( hsu[[ii]]$input['precip'] +
                    ##                   hsu[[ii]]$flux['qsf_rz'] / ts$sub_step,
                    ##                  hsu[[ii]]$input['pet']/hsu[[ii]]$param['srz_max'],
                    ##                  hsu[[ii]]$state['srz'],
                    ##                  ts$sub_step )
                    hsu[[ii]]$state['srz'] <- min(tilde_rz,hsu[[ii]]$param['srz_max'])
                    if( hsu[[ii]]$state['srz'] <= 0 ){
                        hsu[[ii]]$flux[c('qrz_sf','qrz_uz')] <- c(tilde_rz - hsu[[ii]]$state['srz'],0)
                    }else{
                        hsu[[ii]]$flux[c('qrz_sf','qrz_uz')] <- c(0,tilde_rz - hsu[[ii]]$state['srz'])
                    }

                    ## solve for unsaturated zone
                    tilde_uz <- fode( hsu[[ii]]$flux['qrz_uz']/ts$sub_step,
                                     1 / (hsu[[ii]]$param['td'] * hsu[[ii]]$state['ssz']),
                                     hsu[[ii]]$state['suz'],
                                     ts$sub_step)
                    hsu[[ii]]$flux['quz_sz'] <- hsu[[ii]]$state['suz'] + hsu[[ii]]$flux['qrz_uz'] - tilde_uz
                    hsu[[ii]]$state['suz'] <- tilde_uz                    
                   
                }
                if( hsu[[ii]]$type == "channel" ){
                    hsu[[ii]]$state['sch'] <- hsu[[ii]]$input['lsf']
                    alsf[ii] <- 0
                }
            }

            ## PASS 2: in order of saturated zone
            print("pass 2")
            for(ii in sq$sz){
                ## compute the inflows
                lsz_in_prev <- hsu[[ii]]$input['lsz'] ## keep previous input for kinematic solution
                hsu[[ii]]$input['lsz'] <- sum( hsu[[ii]]$parent$sz$x * alsz[hsu[[ii]]$parent$sz$j] ) / hsu[[ii]]$prop['area']
                if( hsu[[ii]]$type == "hillslope" ){
                    lsz_sat <- c(lsz_in_prev,hsu[[ii]]$input['lsz'],
                                 hsu[[ii]]$output['lsz'])
                    lsz_sat[lsz_sat > hsu[[ii]]$state['lsz_max']] <- hsu[[ii]]$state['lsz_max']
                    ##lsz_sat <- pmin( c(lsz_in_prev,hsu[[ii]]$input['lsz'],hsu[[ii]]$output['lsz']),hsu[[ii]]$state['lsz_max'])
                    qbar <- (sum(lsz_sat) + max(lsz_sat))/4
                    cbar <- (qbar*hsu[[ii]]$prop['delta_x'])/ hsu[[ii]]$param['m']
                    lambda <- sz_opt$omega + sz_opt$theta*cbar*ts$sub_step/hsu[[ii]]$prop['delta_x']
                    lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*ts$sub_step/hsu[[ii]]$prop['delta_x']
                    
                    k <- lambda_prime * lsz_sat[3] +
                        (1-lambda_prime)*lsz_sat[1] +
                        cbar*hsu[[ii]]$flux['quz_sz']/hsu[[ii]]$prop['delta_x']
                    lsz <- min( (k - (1-lambda)*lsz_sat[2])/lambda , hsu[[ii]]$state['lsz_max'] ) # new output discharge

                    tilde_sz <- hsu[[ii]]$state['ssz'] +
                        ts$sub_step*(lsz_in_prev + hsu[[ii]]$input['lsz'])/2 -
                        ts$sub_step*(lsz + hsu[[ii]]$output['lsz'])/2 -
                        hsu[[ii]]$flux['quz_sz']
                    hsu[[ii]]$output['lsz'] <- lsz
                    hsu[[ii]]$state['ssz'] <-  max(0, tilde_sz)
                    hsu[[ii]]$flux['qsz_ex'] <- hsu[[ii]]$state['ssz'] - tilde_sz
                    alsz[ii] <- hsu[[ii]]$output['lsz']*hsu[[ii]]$prop['area']
                    
                    ## step 5 - correct the stores
                    hsu[[ii]]$state['ssf'] <- hsu[[ii]]$state['ssf'] +
                        hsu[[ii]]$flux['qrz_sf'] +
                        hsu[[ii]]$flux['qsz_ex']
                    if( hsu[[ii]]$state['ssz'] <=0 ){
                        hsu[[ii]]$state['ssf'] <- hsu[[ii]]$state['ssf'] + hsu[[ii]]$state['suz']
                        hsu[[ii]]$state['suz'] <- 0
                    }
                }
                ## channel
                if( hsu[[ii]]$type == "channel" ){
                    hsu[[ii]]$state['sch'] <- hsu[[ii]]$state['sch'] +
                        ts$sub_step*(lsz_in_prev + hsu[[ii]]$input['lsz'])/2
                    alsf[ii] <- 0
                }
            }

            ## compute channel inflow
            for( ii in 1:length(channel_id) ){
                jj <- channel_id[ii]
                channel_inflow[it,ii] <- hsu[[jj]]$prop['area']*
                    (hsu[[jj]]$state['sch'] + hsu[[jj]]$input['precip'])/ts$step
            }
            
        }
    } ## end of timestep loop

    model$hsu <- hsu

    return( list(model=model,
                 channel_input = channel_inflow ) )
}

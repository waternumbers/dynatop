#' Initialises the hru states for a dynatop simulation
#'
#' @description Initialise the states of the different hru types. these can be time varying such as store of fluxes or fixed properties such as parmeter values and areas
#'
#' @param model a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#' 
#' @name initialise_dynatop
#' @export
initialise_dynatop <- function(model,initial_recharge){

    ## check initial discharge
    if( !is.numeric(initial_recharge) | length(initial_recharge) > 1 | any( initial_recharge < 0 ) ){
        stop("Initial discharge should be a single positive numeric value")
    }
    
    ## function for initialising the hillslope
    fhillslope <- function(model,initial_recharge){
        
        hru_var <- c("id","area","atb_bar","s_bar","precip_input","pet_input")
        param_var <- c("srz_max","srz_0","ln_t0","m","td","tex")
        
        hs <- list()
        for(ii in hru_var){
            hs[[ii]] <- unname( model$hillslope[,ii] )
        }
        for(ii in param_var){
            hs[[ii]] <- unname( model$param[ model$hillslope[,ii] ] )
        }
        
        ## initialise surface storage to 0
        hs$ex <- rep(0,length(hs$id))
        
        ## initialise root zone based on fraction constrained within 0,1
        hs$rz <- pmax( pmin( hs$srz_0, 1) ,0 ) * hs$srz_max
        
        ## initialise based on steady state of saturated zone given constant recharge from unsturated zone
        q_0 <- initial_recharge
        
        ## a least squares solution with total outflow equals total inflow
        ## to saturated zone
        ## Assume all outflow to river?
        X <- rbind( diag(nrow(model$Wsat)) - model$Wsat , colSums(model$Fsat) ) %*% diag(hs$area)
        y <- c( hs$area*q_0, sum(hs$area)*q_0 )
        
        tryCatch({
            hs$lsz <- as.numeric( solve( t(X)%*%X, t(X)%*%y ) )
        },error = function(e){
            warning("Solution for saturated zone initialisation produced negative or zero flows. Using constant recharge across catchment")
            hs$lsz <- rep(q_0,length(hs$id))
        })
        ## catch to keep flows within bounds
        hs$lsz <- pmin( hs$lsz , hs$lsz_max )
        
        ## compute the saturated zone storage deficit based on the flow
        hs$sz <- -hs$m * log(hs$lsz/exp(hs$ln_t0-hs$atb_bar))
        hs$sz <- pmax(  hs$sz, 0 )
                
        ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
        hs$uz <- pmax(0, rep(q_0,length(hs$id)) * hs$td * hs$sz, na.rm=TRUE)
        
        ## take only what we need
        out <- hs[c("ex","rz","uz","sz","lsz",
                    "id","area","precip_input","pet_input",
                    "tex","srz_max","td",
                    "lsz_max","m")]
        return(out)
    }
    
    fchannel <- function(model){
        hru_var <- c("id","area","precip_input","pet_input")
        
        ch <- list()
        for(ii in hru_var){
            ch[[ii]] <- unname( model$channel[,ii] )
        }
        ch$inflow <- rep(0,length(ch$id))
        return(ch)
    }
    
    ## initialise the output
    states <- list(hillslope = fhillslope(model,initial_recharge),
                   channel = fchannel(model))

    return(states)
}


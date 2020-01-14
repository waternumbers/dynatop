#' Functions for a hillslope HRU
#'
#' @description Functions for creating, checking, initialising and evolving a hillslope HRU.
#'
#' @param h a dynamic TOPMODEL model list object
#' @param initial_recharge Initial recharge to the saturated zone in steady state in m/hr
#'
#' @return A list of lists, one list for each HRU type which contains vectors of states and properties for use in computation
#'
#' @rdname hillslope_hru

#' @name hillslope_hru
#' @export
create_hillslope <- function(id){
    flst <- function(x){
        list(id=x,
             type="hillslope",
             external_input=c(precip=0,pet=0),
             internal_input=c(lex=0,lsz=0),
             conectivity=list(lex=list(id=numeric(0),
                                       w=numeric(0)),
                              lex=list(id=numeric(0),
                                       w=numeric(0)))
             output=c(lex=0,lsz=0),
             state=c("ex"=0,
                     "rz"=0,
                     "uz"=0,
                     "sz"=0,
                     "lsz_in"=0,
                     "lsz"=0)
             )
    }
    lapply(id,flst)
}
##     prop_names <- c("area","s_bar","atb_bar","delta_x","band")
##     input_names <- c("precip","pet","lex","lsz")
##     output_names <- c("lex","lsz")
##     param_names <- c("qex_max","srz_max","srz_0","ln_t0","m","td","tex",
##                      "precip_series","pet_series")
##     state_names <- c("ex","rz","uz","sz","lsz_in","lsz")
##     #flux_names <- c("ex_rz","rz_uz","rz_ex","uz_sz","sz_ex")

##     out <- list(id=numeric(0),type="hillslope")
##     for(ii in prop_names){
##         out$prop[[ii]] <- numeric(0)
##     }
##     for(ii in input_names){
##         out$input[[ii]] <- numeric(0)
##     }
##     for(ii in output_names){
##         out$output[[ii]] <- list(id=numeric(0),w=numeric(0),val=numeric(0))
##     }
##     for(ii in param_names){
##         out$param[[ii]] <- numeric(0)
##     }
##     for(ii in state_names){
##         out$state[[ii]] <- numeric(0)
##     }

##     return(out)
## }

#' @name hillslope_hru
#' @export
check_hillslope <- function(h){
    #warning("hillslope check not implimented")
    return(h)
}

#' @name hillslope_hru
#' @export
initialise_hillslope <- function(h,q0){
    print("initialise hillslsope")
    ## check hillslope
    h <- check_hillslope(h)

    ## initialise surface storage to 0
    h$state$ex <- 0

    ## initialise root zone based on fraction constrained within 0,1
    h$state$rz <- max( min( h$param$srz_0, 1) ,0 ) *
        h$param$srz_max

    ## initialise based on steady state of saturated zone given constant recharge from unsturated zone

    ## maximum lateral flow from saturated zone
    #print(h)
    h$state$lsz_max <- exp(h$param$ln_t0)*h$prop$s_bar # for unit width
    h$state$lsz_max <- h$state$lsz_max / h$prop$delta_x # for specified width

    ## initialise
    h$state$lsz <- h$state$lsz_in <- q0
    h$state$lsz <- min(h$state$lsz,h$state$lsz_max)

    ## compute the deficit
    gamma <- h$prop$atb_bar - h$param$ln_t0
    h$state$sz <- h$param$m*(gamma + log(h$state$lsz))
    h$state$sz <- max(0,h$sz)

    ## unsaturated storage by inverse of eqn for q_uz in Beven & Wood 1983
    h$state$uz <- max(0, q0 * h$param$td * h$state$sz, na.rm=TRUE)

    return(h)

}

#' @name hillslope_hru
#' @export
evolve_hillslope <- function(h,delta_t,sz_opt){

    ## he <- h #$state <- list(a="dummy")
    ## for(jj in names(he$output)){
    ##     he$output[[jj]]$val <- rnorm(length(he$output[[jj]]$id))
    ## }
    
    ## return(he)
    ## simple function for solving ODE
    if(h$id==22) browser()
    fode <- function(a,b,x0,t){
        if(b<1e-10){b <- 1e-10}
        x0*exp(-b*t) + (a/b)*(1-exp(-b*t))
    }
    
    ## Step 1: Distribute any surface storage downslope
    ex <- fode(h$input$lex/delta_t,
               1/h$param$tex,
               h$state$ex,delta_t)
    lex <- ex + h$input$lex - h$state$ex
    h$state$ex <- ex
    h$out$lex$val <- lex*h$out$lex$w

    ## Step 2: solve the root zone for hillslope elements
    q_ex_rz <- min(h$state$ex, h$param$qex_max*delta_t)

    ## solve ODE
    tilde_rz <- fode( h$input$precip + (q_ex_rz/delta_t),
                     h$input$pet/h$param$srz_max,
                     h$state$rz,delta_t)
    ## new storage value
    h$state$rz <- min(tilde_rz,h$param$srz_max)
    ## split root zone flow
    tmp <- tilde_rz - h$state$rz
    saturated_index <- h$state$sz <= 0 # which areas are saturated
    q_rz_ex <- tmp * saturated_index
    q_rz_uz <- tmp * !saturated_index

    ## Step 3: Unsaturated zone

    ## solve ODE
    tilde_uz <- fode( q_rz_uz/delta_t,
                     1 / (h$param$td * h$state$sz),
                     h$state$uz,delta_t )
    q_uz_sz <- h$state$uz + q_rz_uz - tilde_uz
    h$state$uz <- tilde_uz


    ## Step 4: Solve saturated zone
<<<<<<< HEAD
    qbar <- (h$state$lsz+h$state$lsz_in+h$input$lex+
             max(h$state$lsz,h$input$lex))/4
=======
    qbar <- (h$state$lsz+h$state$lsz_in+h$input$lsz+
             max(h$state$lsz,h$input$lsz))/4
>>>>>>> d6912ae65ad3ee03ee028b53a85f93525320e34d
    cbar <- (qbar*h$prop$delta_x)/h$param$m

    lambda <- sz_opt$omega + sz_opt$theta*cbar*delta_t/h$prop$delta_x
    lambda_prime <- sz_opt$omega + (1-sz_opt$theta)*cbar*delta_t/h$prop$delta_x

    k <- lambda_prime*h$state$lsz +
        (1-lambda_prime)*min(h$state$lsz_in,h$state$lsz_max) +
        cbar*q_uz_sz/h$prop$delta_x
    
    lsz <- (k - lambda_prime*min(h$input$lsz,h$state$lsz_max))/lambda
    
    if(lsz > h$state$lsz_max){lsz  <- h$state$lsz_max}

    ## do storage calc
<<<<<<< HEAD
    tilde_sz <- h$sz + delta_t*( (h$state$lsz+lsz)/2 -
                                 (h$input$lsz+h$state$lsz_in)/2 ) -
        h$flux$uz_sz

    h$state$sz <- max(0,tilde_sz)
    h$flux$sz_ex <- h$state$sz - tilde_sz
=======
    tilde_sz <- h$state$sz + delta_t*( (h$state$lsz+lsz)/2 -
                                 (h$input$lsz+h$state$lsz_in)/2 ) -
        q_uz_sz
    
    if(tilde_sz < 0){
        h$state$sz <- 0
    }else{
        h$state$sz <- tilde_sz
    }
    
    q_sz_ex <- h$state$sz - tilde_sz
>>>>>>> d6912ae65ad3ee03ee028b53a85f93525320e34d

    ## replace states
    h$state$lsz <- lsz
    h$state$lsz_in <- h$input$lsz

    ## redistribute
    h$output$lsz$val <-  h$output$lsz$w *h$state$lsz

    ## step 5 - correct the stores
    saturated_index <- h$state$sz <= 0
    h$state$ex <- h$state$ex + q_rz_ex +
        (q_sz_ex + h$state$uz)*saturated_index
    h$state$uz <- h$state$uz * !saturated_index

    return(h)
}

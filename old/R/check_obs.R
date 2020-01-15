#' Function to check the observations provided to drive a dynamic TOPMODEL run
#'
#' @description This function makes some basic consistency checks on the observed data and summaries of the model
#'
#' @param obs xts object containing the forcing data
#' @param req_series a character vector of series that should be in obs
#' @param sub_step duration of sub step in model evaluation in hours, default value of NULL results in data time step
#' @return Returns a list containing the time interval for the main and sub step as well and the number of substeps
#'
#' @details The checks performed and basic 'sanity' checks. they do not check for the logic of parameter values nor the consistncy of states and parameters. The number of sub steps is
# #' $$ \max\left(1,floor\left(\frac{data time step}{requested substep}\right)\right) $$
#' 
#' @export
check_obs <- function(obs,req_series,sub_step=NULL){

    ## check types
    if(!is.xts(obs)){ stop("observations should be an xts object") }
    if(!is.vector(req_series) | !all(sapply(req_series,class)=='character') ){ stop("req_series should be a character vector") }
    if( (!is.null(sub_step) & !is.numeric(sub_step)) | length(sub_step)>1 ){ stop("sub_step should be a single numeric value or NULL") }

    ## check we have all the series needed
    if( !all( req_series %in% names(obs) ) ){
        stop("Missing input series:",setdiff( req_series , names(obs) ))
    }

    ## check constant time step
    tmp <- diff(as.numeric(index(obs)))
    if( !all( tmp == tmp[1] ) ){
        stop("Time steps in data are not unique")
    }

    ## check all values are finite
    if( !all(is.finite(obs[,req_series])) ){
        stop("There are non finite values in the required time series")
    }
    
    ## work out time steps for use in simulation
    ts <- list()
    ts$step <- tmp[1]/(60*60) # hours
    if(is.null(sub_step)){sub_step <- ts$step}
    ts$n_sub_step <- max(1,floor(ts$step/sub_step)) # dimensionless
    ts$sub_step <- ts$step / ts$n_sub_step

    return(ts)
}


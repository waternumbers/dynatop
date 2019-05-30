#' Function to check the input data to a dynamic TOPMODEL run
#'
#' @description Describe the checks. Describes where the data comes from. Does basic sanity checks NOT a full check on for example paramter values
#'
#' @param hru details of the hydrological units
#' @param param named parameter vector
#' @param Wsat weighting matrix for saturated zone flow
#' @param Wsurf weighting matric for surface flow
#' @param obs_data xts object of observed data used for inputs
#' @param init_discharge initial discharge for hillslope initialisation
#' @param sim_time_step simulation time step in hours
#' @return Logical if tests passed, otherwise fails
#'
#' @export
check_input <- function(model,param,obs_data,init_discharge,sim_time_step){

    ## check object types
    if(!is.data.frame(model$hru)){stop("hru should be a data.frame")}
    if(!is.vector(param)){stop("param should be a vector")}
    if(!is.matrix(model$Wsurf)){stop("Wsurf should be a matrix")}
    if(!is.matrix(model$Wsat)){stop("Wsat should be a matrix")}
    if(!is.matrix(model$Fsurf)){stop("Fsurf should be a matrix")}
    if(!is.matrix(model$Fsat)){stop("Fsat should be a matrix")}
    
    if(!is.xts(obs_data)){stop("param should be a vector")}
    if(!is.numeric(init_discharge)){stop("initial_discharge should be numeric")}
    if(!is.numeric(sim_time_step)){stop("sim_time_step should be numeric")}

    ## check hru
    ## has required variables
    ## variables are required types
    ## has unique ids
    ## check only recognised types
    ## check only appropriate parameters set
    ## check string not factors


    ## check param
    ## is numeric
    ## has all names in hru$param

    ## check Wsurf
    ## is numeric
    ## has correct id in correct order in column and row
    ## rows sum to 1
    ## check no flow back from channel

    ## check Wsat
    ## is numeric
    ## has correct id in correct order in column and row
    ## rows sum to 1
    ## check no flow back from channel

    ## check obs_data
    ## has unique timestep
    ## has multiple timesteps
    ## has all variables

    ## check init_discharge
    ## check posistive and length 1

    ## check sim_time_step
    ## check multiplies to timestep and <1

    return(TRUE)
}


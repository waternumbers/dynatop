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
dynatop <- R6::R6Class(
    "dynatop",
    public = list(
        initialize = function(hru_list) {
            if( check_hrus(hru_list) ){
                out <- list()
                ## then passed checks and build model
                out[[1]] <- private$hs_init()
                out[[2]] <- private$ch_init()
            }
            print(out)
        }
    ),
    active=list(
        observations = function(value){
            if(missing(value)){
                ##print
            }else{
                ## else check and import
            }
        },
        parameters = function(value){
            if(missing(value)){
                ##print
            }else{
                ## else check and change
            }
        },
        simulate = function(value){
            if(missing(value)){
            }else{
            }
        }    
    ),
    private = list(
        ## inputs to the system
        par = NULL,
        obs = NULL,
        ## states of the system
        hs=list(),
        ch=list(),
        ## initialisation functions
        hs_init = function(){},
        ch_init = function(){},
        ## time step functions
        hs_step = function(ii){}
    )
)

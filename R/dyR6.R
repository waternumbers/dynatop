#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
dynatop <- R6::R6Class(
    "dynatop",
    public = list(
        initialize = function(mdl, initial_recharge, verbose=FALSE,
                              use_states=FALSE,delta = 1e-13){
            ## check the model
            private$info$data_series <- dynatop::check_model(mdl,verbose,use_states,delta)
            
            ## initialise the model
            if( !use_states ){
                mdl <- dynatop::initialise(mdl,initial_recharge)
            }
            
            ## convert to lists for simulation
            ## this creates hillslope, channel, sqnc and lateral_flux
            private$model <- mdl
            
            invisible(self)
        },
        add_data = function(obs_data,sim_time_step=NULL){
            self$clear_data()
            ## check input and get model timestep
            private$info$ts <- dynatop::check_obs(obs_data,private$info$data_series,
                                                  sim_time_step)
            private$time_series$obs_data <- as.matrix(obs_data)
            private$time_series$index <- index(obs_data)
        },
        clear_data = function(){
            
            private$time_series <- list()
            private$info$ts <- list()
        },
        sim = function(mass_check=FALSE,return_states=NULL,
                       sz_opt=list(omega=1,theta=1)){
            browser()
            out <- dynatop::dynatop_sim(private$model,
                                        private$time_series$obs_data,
                                        private$time_series$index,
                                        private$info$ts,
                                        mass_check,
                                        return_states,
                                        sz_opt)
            browser()
            private$model <- out$model
            private$time_series$channel_input <- out$channel_input
            if(mass_check){
                private$time_series$mass_check <- out$mass_errors
            }
            if( "return_states" %in% names(out) ){
                private$return_states <- out$return_states
            }
        }
    ),
    private = list(
        ## stores of data
        version = 0.1,
        time_series = list(),
        info=list(),
        model=NULL
    )
)

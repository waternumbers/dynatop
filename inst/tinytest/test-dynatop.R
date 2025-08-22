library(tinytest)
library(dynatop)

## test baseline simulation
info_str <- "default"
expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux,n_thread=1)
    y_exp <- dt$get_output()
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## test multiple cores on baseline setup
info_str <- "multicore test"
expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux,n_thread=10)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)
expect_equal(dt$get_output(), y_exp, info=info_str)

## check sub stepping
info_str <- "substeping"
expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux,sub_step=300)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## ###############################################
## check saturated zone types

## constant velocity
info_str <- "sz: cnst"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sz$type <- "cnst"; h$sz$parameters <- c("h_szmax" = 0.1, "v_sz" = 0.01); h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## bounded exponential
info_str <- "sz: bexp"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sz$type <- "bexp"; h$sz$parameters["h_szmax"] <- 0.5; h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux)
    error_bexp <- max(abs(dt$get_mass_errors()[,6]))
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## double exponential
info_str <- "sz: dexp"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sz$type <- "dexp"; h$sz$parameters[c("m2","omega")] <- c(0.1,0.5);h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)


## #################################################
## check surface types and rafs

## constant with rafs
info_str <- "sf: cnst with raf"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sf$parameters[c("s_raf","t_raf")] <- c(100,10*60*60);h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    raf_error <- max(abs(dt$get_mass_errors()[,6]))
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## Manning with shallow water approx
info_str <- "sf: Mannings"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sf$type <- "kin";h$sf$parameters["n"] <- 0.03;h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    kin_sf_error <- max(abs(dt$get_mass_errors()[,6]))
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## compound channel with two velocities
info_str <- "sf: compound"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){
                      h$sf$parameters <- c("v_sf_1" =  as.numeric(h$sf$parameters["v_sf"]),
                                           "s_1" = Inf,
                                           "v_sf_2" = 1e-100)
                      h$sf$type <- "comp"
                      h
                  })
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    comp_sf_error <- max(abs(dt$get_mass_errors()[,6]))
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## arbitarty relationship
info_str <- "sf: arbitary relationship"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){
                      h$sf$parameters <- setNames(c(seq(0,1000,length=10),
                                                    seq(0,100,length=10)),
                                                  c(paste0("area_",1:10),
                                                    paste0("flow_",1:10)))
                      h$sf$type <- "arb_kin"
                      h
                  })
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## Power law
info_str <- "sf: powerlaw"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){
                      h$sf$parameters <- c("sc" = 10,
                                           "pwr" = 1.5,
                                           "s_raf" = 0,
                                           "t_raf" = 999)
                      h$sf$type <- "power_law"
                      h
                  })
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## trapezoid
info_str <- "sf: trapezoid"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){
                      h$sf$parameters <- c("n" = 0.03,
                                           "bank_slope" = 1,
                                           "bed_width" = 5)
                      h$sf$type <- "mct"
                      h
                  })
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)

## rectangular and trapezoid
info_str <- "sf: rectangular and trapezoid"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){
                      h$sf$parameters <- c("n" = 0.03,
                                           "b_lower" = 5,
                                           "tan_alpha" = 0.01,
                                           "q_crit" = 20)
                      h$sf$type <- "mct_rect"
                      h
                  })
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux)
}, info=info_str)
expect_true({ all(dt$get_output() > 0) }, info=info_str)
expect_true({ max(abs(dt$get_mass_errors()[,6])) < 1e-6 }, info=info_str)


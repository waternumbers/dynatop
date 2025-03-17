library(tinytest)
library(dynatop)

expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux)
    error_exp <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ error_exp < 1e-6 })

expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux,sub_step=300)
    error_exp_substep <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ error_exp_substep < 1e-6 })


expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sz$type <- "cnst"; h$sz$parameters <- c("h_szmax" = 0.1, "v_sz" = 0.1); h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    error_cnst <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ error_cnst < 1e-6 })

expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sz$type <- "bexp"; h$sz$parameters["h_szmax"] <- 0.1; h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    error_bexp <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ error_bexp < 1e-6 })


expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sz$type <- "dexp"; h$sz$parameters[c("m2","omega")] <- c(0.1,0.5);h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    dexp_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ dexp_error < 1e-6 })

expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sf$parameters[c("s_raf","t_raf")] <- c(100,10*60*60);h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    raf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ raf_error <1e-6 })

## Dynatop mass errors with kinematic surface are <1e-6", {
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sf$type <- "kin";h$sf$parameters["n"] <- 0.03;h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    kin_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ kin_sf_error < 1e-6 })

## Dynatop mass errors with compound surface are <1e-6"
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){
                      h$sf$parameters <- c("v_sf_1" =  as.numeric(h$sf$parameters["c_sf"]),
                                           "d_sf_1" = 0,
                                           "s_1" = Inf,
                                           "v_sf_2" = 0,
                                           "d_sf_2" = 0)
                      h$sf$type <- "comp"
                      h
                  })
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    comp_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ comp_sf_error < 1e-6 })

## ## there are differences in the initialisation which meant he comparision is only valid with s_1 =0
## test_that("Dynatop cnst and compound solutions are consistent without raf to <1e-3", {
##     data(Swindale)
##     mdl_raf <- Swindale$model$hru
##     mdl_cmp <- Swindale$model$hru
##     for(ii in 1:length(mdl_cmp)){
##         mdl_cmp[[ii]]$sf$parameters <- c("v_sf_1" = 999,
##                                          "d_sf_1" = 0,
##                                          "s_1" = 0,
##                                          "v_sf_2" = as.numeric(mdl_cmp[[ii]]$sf$parameters["c_sf"]), ## needs to be positive else get NaN from C++
##                                          "d_sf_2" = 0)
##         mdl_cmp[[ii]]$sf$type <- "comp"
##     }
##     dt_raf <- dynatop$new(mdl_raf)$add_data(Swindale$obs)
##     dt_raf$initialise()
##     dt_cmp <- dynatop$new(mdl_cmp)$add_data(Swindale$obs)
##     dt_cmp$initialise()

##     s_raf <- dt_raf$get_states()
##     s_cmp <- dt_cmp$get_states()
##     e <- s_cmp - s_raf
##     testthat::expect_lt( max(abs(e)), 1e-8 )

##     dt_raf$sim(Swindale$model$output_flux,vtol=1e-8)
##     dt_cmp$sim(Swindale$model$output_flux,vtol=1e-8)

##     tmp <- max(abs(dt_cmp$get_output() - dt_raf$get_output()))
##     testthat::expect_lt( tmp, 1e-8 )
## })




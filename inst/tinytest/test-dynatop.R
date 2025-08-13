library(tinytest)
library(dynatop)

## testing number of threads works - maybe remove for CRAN
expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux,n_thread=1)
    error_exp <- max(abs(dt$get_mass_errors()[,6]))

    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()
    dt$sim(Swindale$model$output_flux,n_thread=10)
    error_exp10 <- max(abs(dt$get_mass_errors()[,6]))


})
expect_true({ error_exp < 1e-6 })

## check sub stepping
expect_silent({
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux,sub_step=300)
    error_exp_substep <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ error_exp_substep < 1e-6 })

## ###############################################
## check saturated zone types
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

## #################################################
## check surface types adn rafs
expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sf$parameters[c("s_raf","t_raf")] <- c(100,10*60*60);h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    raf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ raf_error <1e-6 })

expect_silent({
    data(Swindale)
    mdl <- lapply(Swindale$model$hru,
                  function(h){h$sf$type <- "kin";h$sf$parameters["n"] <- 0.03;h})
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    kin_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ kin_sf_error < 1e-6 })

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
})
expect_true({ comp_sf_error < 1e-6 })

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
    arb_kin_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ arb_kin_sf_error < 1e-6 })

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
    power_law_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ power_law_sf_error < 1e-6 })

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
    mct_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ mct_sf_error < 1e-6 })

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
    mct_rect_sf_error <- max(abs(dt$get_mass_errors()[,6]))
})
expect_true({ mct_rect_sf_error < 1e-6 })


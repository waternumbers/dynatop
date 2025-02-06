test_that("Dynatop mass errors for exponential profile are <1e-6", {
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors are correctly computed for substeps and less then <1e-6", {
    data(Swindale)
    dt <- dynatop$new(Swindale$model$hru)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux,sub_step=300)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors for constant profile are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sz$type <- "cnst"
        mdl[[ii]]$sz$parameters <- c("h_sz_max" = 0.1, "v_sz" = 0.1)
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors for bounded exponential profile are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sz$type <- "bexp"
        mdl[[ii]]$sz$parameters["h_sz_max"] <- 0.1
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors for double exponential transmissivity profile are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sz$type <- "dexp"
        mdl[[ii]]$sz$parameters[c("m2","omega")] <- c(0.1,0.5)
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors with two path cnst surface are  <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sf$parameters["v_sf_1"] <- mdl[[ii]]$sf$parameters["v_sf_2"]
        mdl[[ii]]$sf$parameters["s_1"] <- 1
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors with kinematic surface are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sf$type <- "kin"
        mdl[[ii]]$sf$parameters <- c("n"=0.03,"s_raf"=0,"t_raf"=999)
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors with power law surface are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sf$parameters <- c("pwr_raf" = 1, as.numeric(mdl[[ii]]$sf$parameters["c_sf"]),
                                     "sc_raf" = 1,
                                     "s_raf" = 0,
                                     "pwr" = 0.5,
                                     "sc" = 5)
        mdl[[ii]]$sf$type <- "power_law"
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})


test_that("Dynatop cnst and compound solutions are consistent without raf to <1e-3", {
    data(Swindale)
    mdl_cnst <- Swindale$model$hru
    mdl_pwr <- Swindale$model$hru
    for(ii in 1:length(mdl_cnst)){
        mdl_cnst[[ii]]$sf$parameters["v_sf_1"] <- 0.002
        mdl_cnst[[ii]]$sf$parameters["s_1"] <- 1.2
        mdl_pwr[[ii]]$sf$type <- "power_law"
        mdl_pwr[[ii]]$sf$parameters <- c("pwr_raf" = 1,
                                         "sc_raf" = as.numeric(mdl_cnst[[ii]]$sf$parameters["v_sf_1"] /
                                             mdl_cnst[[ii]]$properties["Dx"]),
                                         "s_raf" =  as.numeric( mdl_cnst[[ii]]$sf$parameters["s_1"] ),
                                         "pwr" = 1,
                                         "sc" = as.numeric(mdl_cnst[[ii]]$sf$parameters["v_sf_2"] /
                                                           mdl_cnst[[ii]]$properties["Dx"]))
    }
    dt_cnst <- dynatop$new(mdl_cnst)$add_data(Swindale$obs)
    dt_cnst$initialise()
    dt_pwr <- dynatop$new(mdl_pwr)$add_data(Swindale$obs)
    dt_pwr$initialise()

    s_cnst <- dt_cnst$get_states()
    s_pwr <- dt_pwr$get_states()
    e <- s_pwr - s_cnst
    testthat::expect_lt( max(abs(e)), 1e-8 )

    dt_cnst$sim(Swindale$model$output_flux,vtol=1e-8)
    dt_pwr$sim(Swindale$model$output_flux,vtol=1e-8)

    tmp <- max(abs(dt_pwr$get_output() - dt_cnst$get_output()))
    testthat::expect_lt( tmp, 1e-8 )
})




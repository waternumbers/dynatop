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

test_that("Dynatop mass errors for double exponential are <1e-6", {
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

test_that("Dynatop mass errors with RAF are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sz$parameters[c("s_raf","t_raf")] <- c(100,10*60*60)
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
        mdl[[ii]]$sf$parameters["n"] <- 0.03
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors with compound surface are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model$hru
    for(ii in 1:length(mdl)){
        mdl[[ii]]$sf$parameters <- c("v_sf_1" =  as.numeric(mdl[[ii]]$sf$parameters["c_sf"]),
                                     "d_sf_1" = 0,
                                     "s_1" = Inf,
                                     "v_sf_2" = 0,
                                     "d_sf_2" = 0)
        mdl[[ii]]$sf$type <- "comp"
    }
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise()$sim(Swindale$model$output_flux)
    tmp <- max(abs(dt$get_mass_errors()[,6]))
    testthat::expect_lt( tmp, 1e-6 )
})

## there are differences in the initialisation which meant he comparision is only valid with s_1 =0
test_that("Dynatop cnst and compound solutions are consistent without raf to <1e-3", {
    data(Swindale)
    mdl_raf <- Swindale$model$hru
    mdl_cmp <- Swindale$model$hru
    for(ii in 1:length(mdl_cmp)){
        mdl_cmp[[ii]]$sf$parameters <- c("v_sf_1" = 999,
                                         "d_sf_1" = 0,
                                         "s_1" = 0,
                                         "v_sf_2" = as.numeric(mdl_cmp[[ii]]$sf$parameters["c_sf"]), ## needs to be positive else get NaN from C++
                                         "d_sf_2" = 0)
        mdl_cmp[[ii]]$sf$type <- "comp"
    }
    dt_raf <- dynatop$new(mdl_raf)$add_data(Swindale$obs)
    dt_raf$initialise()
    dt_cmp <- dynatop$new(mdl_cmp)$add_data(Swindale$obs)
    dt_cmp$initialise()
   
    s_raf <- dt_raf$get_states()
    s_cmp <- dt_cmp$get_states()
    e <- s_cmp - s_raf
    testthat::expect_lt( max(abs(e)), 1e-8 )

    dt_raf$sim(Swindale$model$output_flux,vtol=1e-8)
    dt_cmp$sim(Swindale$model$output_flux,vtol=1e-8)

    tmp <- max(abs(dt_cmp$get_output() - dt_raf$get_output()))
    testthat::expect_lt( tmp, 1e-8 )
})




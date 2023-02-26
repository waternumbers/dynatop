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

## currently this is to slow for CRAN - run manually it passes
## test_that("Dynatop mass errors for constant profile are <1e-6", {
##     data(Swindale)
##     mdl <- Swindale$model
##     mdl$hillslope$opt <- "cnst"
##     mdl$hillslope$c_sz<- 0.5; mdl$hillslope$D <- 5
##     dt <- dynatop$new(mdl)$add_data(Swindale$obs)
##     dt$initialise()$sim_hillslope()
##     tmp <- max(abs(dt$get_mass_errors()[,6]))
##     testthat::expect_lt( tmp, 1e-6 )
## })

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

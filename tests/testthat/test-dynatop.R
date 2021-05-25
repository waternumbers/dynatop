test_that("Dynatop mass errors for exponential profile are <1e-6", {
    data(Swindale)
    dt <- dynatop$new(Swindale$model)$add_data(Swindale$obs)
    dt$initialise(1e-6)$sim_hillslope()
    tmp <- max(rowSums(dt$get_mass_errors()))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors for constant profile are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model
    mdl$hillslope$c_sz <- 0.5; mdl$hillslope$D <- 5
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise(1e-6)$sim_hillslope()
    tmp <- max(rowSums(dt$get_mass_errors()))
    testthat::expect_lt( tmp, 1e-6 )
})

test_that("Dynatop mass errors for bounded exponential profile are <1e-6", {
    data(Swindale)
    mdl <- Swindale$model
    mdl$hillslope$D <- 5
    dt <- dynatop$new(mdl)$add_data(Swindale$obs)
    dt$initialise(1e-6)$sim_hillslope()
    tmp <- max(rowSums(dt$get_mass_errors()))
    testthat::expect_lt( tmp, 1e-6 )
})


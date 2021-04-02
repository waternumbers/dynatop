rm(list=ls())
.libPaths(c("/home/paul/R/x86_64-suse-linux-gnu-library/4.0/",.libPaths()))
devtools::load_all("./dynatop"); data("Swindale");

mdl <- Swindale$model
## mdl$hillslope <- mdl$hillslope[1,,drop=FALSE]
## mdl$hillslope$id <- as.integer(2)
## mdl$channel <- mdl$channel[1,,drop=FALSE]
## mdl$channel$id <- as.integer(1)
## mdl$flow_direction <- mdl$flow_direction[1,,drop=FALSE]
## mdl$flow_direction$from <- as.integer(2)
## mdl$flow_direction$to <- as.integer(1)
## mdl$flow_direction$frc <- 1
## mdl$rainfall_input <- mdl$rainfall_input[1:2,,drop=FALSE]
## mdl$rainfall_input$id <- as.integer(c(1,2))
## mdl$pet_input <- mdl$pet_input[1:2,,drop=FALSE]
## mdl$pet_input$id <- as.integer(c(1,2))

mdl$param <- c(r_sfmax_default=Inf,
               m_default=0.16, ## 0.05 - 0.6 m
               ln_t0_default=0.746, ## 0.1 - 8 m2/h
               s_rz0_default=0.98,
               s_rzmax_default=0.1,
               v_ch_default=0.4, ## 1000 -5000 m/h
               t_d_default=80*60*60, ## 0.1-120 m/h - are these units correct - or is cobceptualisation other way round?
               c_sf_default=0.4
               )
## bits that I shouldn't have to do...
mdl$options=c("transmisivity_profile"="exponential","channel_solver"="histogram")

m1 <- dynatop$new(mdl)
##Swindale$obs[,c("Rainfall","PET")] <- 0
m1$add_data(Swindale$obs)
m1$initialise(1e-6)
m1$sim_hillslope(keep_states=index(Swindale$obs))
tmp <- m1$get_states(TRUE)
v = sapply(tmp,colMeans)
#matplot(t(v))

mb <- m1$get_mass_errors()
mb <- cbind(mb,"final"=NA)
idx <- 2:nrow(mb)
mb[idx-1,"final"] <- -mb[idx,"initial_state"]
plot(rowSums(mb))










## rm(list=ls())
## devtools::load_all("./dynatop"); data("Swindale");

## Swindale$model$param <- c(r_sfmax_default=Inf,
##                           m_default=0.006, ## 0.05 - 0.6 m
##                           ln_t0_default=0.746, ## 0.1 - 8 m2/h
##                           s_rz0_default=0.98,
##                           s_rzmax_default=0.1,
##                           v_ch_default=0.4, ## 1000 -5000 m/h
##                           t_d_default=80*60*60, ## 0.1-120 m/h - are these units correct - or is cobceptualisation other way round?
##                           c_sf_default=0.4
##                           )
## ## bits that I shouldn't have to do...
## Swindale$model$options=c("transmisivity_profile"="exponential","channel_solver"="histogram")

## m1 <- dynatop$new(Swindale$model)
## ##tmp <- m1$get_model(); mtmp <- dynatop$new(tmp)
## Swindale$obs[,c("Rainfall","PET")] <- 0
## m1$add_data(Swindale$obs[1:50,])
## ##tmp <- m1$get_model(); mtmp <- dynatop$new(tmp)
## m1$initialise(1e-6)
## ##tmp <- m1$get_model(); mtmp <- dynatop$new(tmp)

## m1$sim_hillslope(keep_states=index(Swindale$obs))
## tmp <- m1$get_states(TRUE)
## v = sapply(tmp,colMeans)
## matplot(t(v))

## mb <- m1$get_mass_errors()
## mb <- cbind(mb,"final"=NA)
## idx <- 2:nrow(mb)
## mb[idx-1,"final"] <- -mb[idx,"initial_state"]
## plot(rowSums(mb))

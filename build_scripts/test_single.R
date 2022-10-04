## test some class and inheritance through Rcpp
rm(list=ls())
devtools::load_all()

data("Swindale"); obs <- Swindale$obs; names(obs) <- c("Flow","precip","pet")

mdl <- list(list(
    id = as.integer(0),
    states = setNames(as.numeric(rep(NA,6)), c("s_sf","s_rz","s_uz","s_sz","q_sf","q_sz")),
    properties = c(area = as.numeric(100), width=as.numeric(2), Dx = as.numeric(50), gradient= as.numeric(0.01)),
    sf = list( type="cnstCD", parameters = c(c_sf = 1,d_sf =0, s_raf=0.0,t_raf=10)),
    ##sf = list( type="cnstC_raf", parameters = c(c_sf = 0.1,s_raf=0,t_raf=10)),
    rz = list( type="orig", parameters = c(s_rzmax = 0.05)),
    uz = list( type="orig", parameters = c(t_d = 3)),
    sz = list( type="exp", parameters = c(t_0=1,m=0.001,D=5)),
    initialisation = c(s_rz_0 = 0.75,r_uz_sz_0 = 1e-3),
    sf_flow_direction = list(id=integer(0),fraction=numeric(0)),
    sz_flow_direction = list(id=integer(0),fraction=numeric(0)),
    precip = list(name="precip",fraction=as.numeric(1)),
    pet = list(name="pet",fraction=as.numeric(1))
))

output_flux = data.frame(name = c("q_sf_1","q_sz_1","s_sf","s_rz","s_uz","s_sz","r_sf_rz","r_rz_uz","r_uz_sz","precip","pet","aet"),
                         id = as.integer(rep(0, 12)),
                         flux = c("q_sf","q_sz","s_sf","s_rz","s_uz","s_sz","r_sf_rz","r_rz_uz","r_uz_sz","precip","pet","aet"))

## dt <- dynatop$new(mdl)
## #obs[,"precip"] <- 0
## #obs[,"pet"] <- 0
## dt$add_data(obs) # [1:35,])

## dt$initialise()


## is <- dt$get_states()

##                                         #
## dt$sim(output_flux)
## y <- dt$get_output()

## x <- dt$get_mass_errors()
## range(y$s_sf)

obs$precip <- 1 #[] <- 0
dt <- dynatop$new(mdl)$add_data(obs)$initialise()
dt$sim(output_flux)
dt$plot_output(c("q_sf_1","q_sz_1"))
##$sim(output_flux)
##cnst <- dt$get_output()

## mdl[[1]]$sf$type <- "cnstC_raf"
## dt <- dynatop$new(mdl)$add_data(obs)$initialise()$sim(output_flux)
## cnst_raf <- dt$get_output()

## par(mfrow=c(2,1)); plot( merge(cnst$s_sf, cnst_raf$s_sf) );plot( cnst$s_sf - cnst_raf$s_sf );

## range(x$error)

#plot( merge( cnst$q_sf, cnst_raf$q_sf) )
#plot( merge( cnst$r_sf_rz, cnst_raf$r_sf_rz) )

head(dt$get_output(c("q_sf_1","q_sz_1","r_sf_rz","r_rz_uz","r_uz_sz","s_sz","s_uz","s_rz","s_sf")),30)


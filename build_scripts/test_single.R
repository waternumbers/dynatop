## test some class and inheritance through Rcpp
rm(list=ls())
devtools::load_all("../")

data("Swindale"); obs <- Swindale$obs; names(obs) <- c("Flow","precip","pet")

mdl <- list(list(
    id = as.integer(0),
    states = setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz")),
    properties = c(area = as.numeric(100), width=as.numeric(2), gradient= as.numeric(0.01)),
    sf = list( type="cnst", parameters = c(c_sf = 0.1)),
    rz = list( type="orig", parameters = c(s_rzmax = 0.05)),
    uz = list( type="orig", parameters = c(t_d = 3600)),
    sz = list( type="bexp", parameters = c(t_0=1,m=0.01,D=5.0)),
    initialisation = c(s_rz_0 = 0.75,r_uz_sz_0 = 1e-3),
    sf_flow_direction = list(id=integer(0),fraction=numeric(0)),
    sz_flow_direction = list(id=integer(0),fraction=numeric(0)),
    precip = list(name="precip",fraction=as.numeric(1)),
    pet = list(name="pet",fraction=as.numeric(1))
))

output_flux = data.frame(name = c("q_sf_1","q_sz_1","s_sf","
s_rz","s_uz","s_sz","r_sf_rz","r_rz_uz","r_uz_sz","precip","pet","aet"),
                         id = as.integer(rep(0, 12)),
                         flux = c("q_sf","q_sz","s_sf","s_rz","s_uz","s_sz","r_sf_rz","r_rz_uz","r_uz_sz","precip","pet","aet"))

dt <- dynatop$new(mdl)
#obs[,"precip"] <- 0
#obs[,"pet"] <- 0
dt$add_data(obs)#;[1:2,]) #[1:2,,drop=F])
dt$initialise()
dt$get_states()
                                        #
dt$sim(output_flux)
y <- dt$get_output()
x <- dt$get_mass_errors()
range(y$s_sf)

range(x$error)

rm(list=ls())

## ----load_library-------------------------------------------------------------
devtools::load_all() ##library("dynatop")


## ----tempory_dir--------------------------------------------------------------
demo_dir <- tempdir()

## ----initialisation-----------------------------------------------------------
ctch <- dynatopGIS$new(file.path(demo_dir,"example.tif"))

## ----data_files---------------------------------------------------------------
dem_file <- system.file("extdata", "gis","SwindaleDTM40m.tif", package="dynatop", mustWork = TRUE)
channel_file <- system.file("extdata", "gis","SwindaleRiverNetwork.shp", package="dynatop", mustWork = TRUE)


## ----add_catchment------------------------------------------------------------
dem <- terra::rast(dem_file)
dem <- terra::extend(dem,1) ## pad with NA values
catchment_outline <- terra::ifel(is.finite(dem),1,NA)
ctch$add_catchment(catchment_outline)


## ----add_dem------------------------------------------------------------------
ctch$add_dem(dem)


## ----channel_current----------------------------------------------------------
sp_lines <- terra::vect(channel_file)
head(sp_lines)


## ----channel_properties-------------------------------------------------------
property_names <- c(name="identifier",
                    endNode="endNode",
                    startNode="startNode",
                    length="length")
chn <- convert_channel(sp_lines,property_names)


## ----add_channel--------------------------------------------------------------
ctch$add_channel(chn)


## ----list_layers--------------------------------------------------------------
ctch$get_layer()


## ----plot---------------------------------------------------------------------
ctch$plot_layer("dem", add_channel=TRUE)


## ----get_layer----------------------------------------------------------------
ctch$get_layer("dem")


## ----sink_fill----------------------------------------------------------------
ctch$sink_fill()

terra::plot( ctch$get_layer('filled_dem') - ctch$get_layer('dem'),
            main="Changes to height")


ctch$compute_properties()


## ----plot_atb-----------------------------------------------------------------
## plot of topographic index (log(a/tan b))
ctch$plot_layer('atb')


## ----flow_length--------------------------------------------------------------
## ctch$compute_flow_lengths(flow_routing="shortest")


## ----flow_length_plot---------------------------------------------------------
ctch$get_layer()
#ctch$plot_layer("shortest_flow_length")


## ----extract_filled-----------------------------------------------------------
tmp <- ctch$get_layer("filled_dem")


## ----height layer-------------------------------------------------------------
## T
tmp <- terra::ifel(tmp<=500,0,1)


## ----add_height_layer---------------------------------------------------------
ctch$add_layer(tmp, "greater_500")
ctch$get_layer()


## ----atb_split----------------------------------------------------------------
ctch$add_layer( ctch$classify("atb_20","atb",cuts=20) )
ctch$plot_layer("atb_20")

## ----atb_20_band--------------------------------------------------------------
ctch$add_layer( ctch$combine_classes("atb_20_band",c("atb_20","band")) )
ctch$plot_layer("atb_20_band")


## ----atb_20_band_burn---------------------------------------------------------
terra::plot(ctch$combine_classes("atb_20_band_500",pairs=c("atb_20","band"),burns="greater_500"))
##ctch$plot_layer("atb_20_band_500")


## ----see_class----------------------------------------------------------------
##head( ctch$get_method("atb_20_band_500")$groups )


## ----model_atb_split----------------------------------------------------------
ctch$create_model(file.path(demo_dir,"new_model"),"atb_20")
##ctch$create_model(file.path(".","/build_scripts","vignette_debug","new_model"),"atb_20")


## ----model files--------------------------------------------------------------
list.files(demo_dir,pattern="new_model*")


## ----setup--------------------------------------------------------------------
#library(dynatop)
data("Swindale")


## ----data_loaded--------------------------------------------------------------
#names(Swindale)


## ----sep----------------------------------------------------------------------
swindale_model <- readRDS(file.path(demo_dir,"new_model.rds")) #Swindale$model

##rm(list=ls())
##devtools::load_all()
##swindale_model <- readRDS(file.path(".","/build_scripts","vignette_debug","new_model.rds"))


swindale_obs <- Swindale$obs


## ----model_parts--------------------------------------------------------------
names(swindale_model)


## ----set_map------------------------------------------------------------------
swindale_model$map <- file.path(demo_dir,"new_model.tif")
##swindale_model$map <- file.path(".","/build_scripts","vignette_debug","new_model.tif")
##system.file("extdata","mdl","Swindale.tif",package="dynatop",mustWork=TRUE)


## ----obs----------------------------------------------------------------------
head(swindale_obs)


## ----set_obs_names------------------------------------------------------------
head(swindale_model$hru[[1]]$precip)
head(swindale_model$hru[[1]]$pet)


## ----change_param-------------------------------------------------------------
hru <- swindale_model$hru
for(ii in 1:length(hru)){
    if(!("endNode" %in% names(hru[[ii]]$class))){
        ## then HRU is not a channel
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["m"] <- 0.0063
        hru[[ii]]$sz$parameters["t_0"] <- exp(7.46)
        ## unsaturated zone parameters
        hru[[ii]]$uz$parameters["t_d"] <- 8*60*60
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.1
        ## surface parameters
        hru[[ii]]$sf$parameters["c_sf"] <- 0.4
    }else{
        ## then HRU is a channel - set so no subsurface response
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["t_0"] <- 0.000
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.001
        ## surface parameters
        hru[[ii]]$sf$parameters["c_sf"] <- 0.8
    }
    ## initialisation parameters
    hru[[ii]]$initialisation["s_rz_0"] <- 0.98
    hru[[ii]]$initialisation["r_uz_sz_0"] <- 1.755582e-07 ## initial outflow divided by catchment area
}


## ----create_object------------------------------------------------------------
rm(list=ls())
devtools::load_all()
data("Swindale")
swindale_model = Swindale$model
#swindale_model <- readRDS(file.path(".","/build_scripts","vignette_debug","new_model.rds"))
hru <- swindale_model$hru
for(ii in 1:length(hru)){
    hru[[ii]]$sz$type <- "exp"
##    hru[[ii]]$sz$parameters <- hru[[ii]]$sz$parameters,"h_max"=2)
    if(!("endNode" %in% names(hru[[ii]]$class))){
        ## then HRU is not a channel
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["m"] <- 0.04# 0.0063
        hru[[ii]]$sz$parameters["t_0"] <- 0.135 #exp(7.46) ##0.135 #exp(-1) #5) #)
        ## unsaturated zone parameters
        hru[[ii]]$uz$parameters["t_d"] <- 8*60*60
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.1
        ## surface parameters
        hru[[ii]]$sf$parameters["v_sf"] <- 0.4
        ## test of raf
        ##hru[[ii]]$sf$parameters["s_raf"] <- 9000
    }else{
        ## then HRU is a channel - set so no subsurface response
        ## saturated zone parameters
        hru[[ii]]$sz$parameters["t_0"] <- 0.000
        ## root zone parameters
        hru[[ii]]$rz$parameters["s_rzmax"] <- 0.1
        ## surface parameters
        hru[[ii]]$sf$parameters["v_sf"] <- 0.8
    }
    ## initialisation parameters
    hru[[ii]]$initialisation["s_rz_0"] <- 0.98
    hru[[ii]]$initialisation["r_uz_sz_0"] <- 1.755582e-07 ## initial outflow divided by catchment area
}

## ctch_mdl <- dynatop$new(hru) #,map=swindale_model$map)
## ## ----add_data-----------------------------------------------------------------
## #data("Swindale")
swindale_obs <- Swindale$obs
ctch_mdl <- dynatop$new(hru)
ctch_mdl$add_data(swindale_obs)
ctch_mdl$initialise()
outDefn <- swindale_model$output_flux
system.time({ ctch_mdl$sim(outDefn) })
o1 <- ctch_mdl$get_states()
y1 <- ctch_mdl$get_output()



## #swindale_obs <- swindale_obs["2009-11-16::2009-11-19"]
## #swindale_obs$precip <- 0
## ctch_mdl$add_data(swindale_obs)
## ## ----initialise---------------------------------------------------------------
## ctch_mdl$initialise()
## st <- list(initial=ctch_mdl$get_states())
## #

                                        #ctch_mdl$plot_state("s_sz")
outID <- 0:9915 #c(9740,9744,9755,9761,9767,9774)
outDefn <- data.frame( name = c(paste0(outID,"_q_sf"),
                                 paste0(outID,"_q_sf_in"),
                                 paste0(outID,"_q_sz"),
                                 paste0(outID,"_q_sz_in"),
                                 paste0(outID,"_precip"),
                                 paste0(outID,"_pet")),
                      id = rep(outID,6),
                      flux = c(rep("q_sf",length(outID)),
                               rep("q_sf_in",length(outID)),
                               rep("q_sz",length(outID)),
                               rep("q_sz_in",length(outID)),
                               rep("precip",length(outID)),
                               rep("pet",length(outID))),
                      scale=1 )

devtools::load_all() ##library("dynatop")
outDefn <- swindale_model$output_flux
gc()
ctch_mdl <- dynatop$new(hru)
ctch_mdl$add_data(swindale_obs)
system.time({ ctch_mdl$initialise() })
s1 <- ctch_mdl$get_states()
system.time({ ctch_mdl$sim(outDefn,sub_step=900) })
o1 <- ctch_mdl$get_states()
y1 <- ctch_mdl$get_output()
print(max(abs(ctch_mdl$get_mass_errors()[,6])))
x11();plot(y1); lines(Swindale$obs$flow)

## gc()
## n_thread <- 3
## ctch_mdl <- dynatop$new(hru)
## ctch_mdl$add_data(swindale_obs[tdx,])
## system.time({ ctch_mdl$initialise(n_thread=n_thread) })
## s3 <- ctch_mdl$get_states()
## system.time({ ctch_mdl$sim(outDefn, n_thread=n_thread) })
## o3 <- ctch_mdl$get_states() #output()
## y3 <- ctch_mdl$get_output()

## all(s1==s3)
## all(o1==o3)
## all(y1==y3)

## plot(y1);points(y3)

## any(abs(o1-o3)>1e-6)
## which(colSums(abs(o1-o3)>1e-6)>1)




## ctch_mdl$add_data(swindale_obs)
## ## ----initialise---------------------------------------------------------------
## ctch_mdl$initialise()
## st <- list(initial=ctch_mdl$get_states())
## ##9740 9744 9755 9761 9767 9774

## print( system.time({sim1 <- ctch_mdl$sim(swindale_model$output_flux, keep_states= index(swindale_obs)[1:3])$get_output()}) )
## st[["sim1"]] <- ctch_mdl$get_states(rec=TRUE)

## ctch_mdl$initialise()
## print(system.time({sim2 <- ctch_mdl$sim(swindale_model$output_flux, keep_states= index(swindale_obs)[1:3],n_thread=3)$get_output()}))
## st[["sim2"]] <- ctch_mdl$get_states(rec=TRUE) #mass_errors()

## theSame <- rep(NA,length(st[["sim1"]]))
## hasSurface <- rep(NA,length(st[["sim1"]]))
## tp <- 1
## for(ii in 1:length(st[["sim1"]])){
##     theSame[ii] <- ( st[["sim1"]][[tp]][[ii]]$id == st[["sim2"]][[tp]][[ii]]$id ) &
##         all( st[["sim1"]][[tp]][[ii]]$states == st[["sim2"]][[tp]][[ii]]$states )
##     hasSurface[ii] <- (st[["sim1"]][[tp]][[ii]]$states["s_sf"]>0) +
##         2*(st[["sim2"]][[tp]][[ii]]$states["s_sf"]>0)
## }



## ctch_mdl$initialise()
## print(system.time({ sim3 <- ctch_mdl$sim(swindale_model$output_flux,sub_step=60,n_thread=10)$get_output() }))
## st[["sim3"]] <- ctch_mdl$get_states() #mass_errors()

## out <- Reduce(merge,list(swindale_obs,sim1,sim2,sim3))
## names(out) <- c(names(swindale_obs),'sim_1','sim_2',"sim_3")
## plot(out[,c('flow','sim_1','sim_2',"sim_3")], main="Discharge",ylab="m3/s",legend.loc="topright")
## ##plot(out[,c('sim_1','sim_2',"sim_3")], main="Discharge",ylab="m3/s",legend.loc="topright")

## lapply(st,function(s){sapply(s,range)})
## ## ----mass_check---------------------------------------------------------------
## mb <- ctch_mdl$get_mass_errors()
## plot( mb[,6] , main="Mass Error", ylab="[m^3]")


## ## converting to a strin for alternative usage
## str <- supressMessages({jsonlite::toJSON(hsc,pretty=TRUE,keep_vec_names = TRUE)})


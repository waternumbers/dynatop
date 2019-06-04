rm(list=ls())

devtools::load_all('./dynatop/dynatop')

data(brompton)

## Examine the November 2012 event that flooded the village (see Metcalfe et al., 2017)
sel <- "2012-11-23 12:00::2012-12-01"
obs_data <- brompton$obs_data[sel,]
obs_data[!is.finite(obs_data)] <- 0

obs_data[obs_data[,'rain']>0,'pe'] <- 0
obs_data[,'rain'] <- obs_data[,'rain']
## set up a model
model <- list(hru = brompton$hru,
              Wsat = t( brompton$Wsat[-1,-1] ),
              Wsurf = t( brompton$Wsurf[-1,-1] ),
              Fsat = t(brompton$Wsat[-1,1,drop=FALSE]),
              Fsurf = t(brompton$Wsurf[-1,1,drop=FALSE])
              )

## alter discharge to m3/s
obs_data[,'q'] <- obs_data[,'q']*sum(model$hru$area) /3600

## alter to reflect previous example
param <- brompton$param
param['srz_max_default'] <- 0.1
param['td_default'] <- 29
param['srz_0_default'] <- 0.98
param['m_default'] <- 0.0051
param['ln_t0_default'] <- 8.1
param['Tex_default'] <- 0.01

## get rid of unused v_of
param <- param[setdiff(names(param),'vof_default')]
model$hru <- model$hru[,setdiff(names(model$hru),'vof')]
model$hru[['Tex']] <- 'Tex_default'

## set hru tabel to characters not factors
i <- sapply(model$hru, is.factor)
model$hru[i] <- lapply(model$hru[i], as.character)

## convert factors to character vectors
init_discharge <- as.numeric(obs_data[1,'q'])
sim_time_step <- 0.25
out <- dynatop(model, param, obs_data, as.numeric(obs_data[1,'q']), 0.25)
tmp <- merge(obs_data,out)
plot(tmp[,c('q','X1')])


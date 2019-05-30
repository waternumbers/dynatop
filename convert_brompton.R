## convert the old brompton data into the new format
rm(list=ls())
library(xts)

load("./old/brompton_old.rda")
hru <- data.frame(id = as.numeric(brompton$groups[,'id']),
                  label = as.character(brompton$groups[,'id']),
                  type = ifelse(is.finite(brompton$groups[,'chan.no']),'channel','hillslope'),
                  area = as.numeric(brompton$groups[,'area']),
                  atb_bar = as.numeric(brompton$groups[,'atb.bar']),
                  precip_input = rep("rain",nrow(brompton$groups)),
                  pet_input = rep("pe",nrow(brompton$groups)),
                  srz_max = rep("srz_max_default",nrow(brompton$groups)),
                  srz_0 = rep("srz_0_default",nrow(brompton$groups)),
                  ln_t0 = rep("ln_t0_default",nrow(brompton$groups)),
                  m = rep("m_default",nrow(brompton$groups)),
                  td = rep("td_default",nrow(brompton$groups)),
                  sd_max = rep("sd_max_default",nrow(brompton$groups)),
                  vof = rep("vof_default",nrow(brompton$groups))
                  )
param <- c("srz_max_default" = 0.05,
           "srz_0_default" = 0.99,
           "ln_t0_default" = 19,
           "m_default" = 0.004,
           "td_default"=20,
           "vof_default"=100)
Wsat <- brompton$weights
colnames(Wsat) <- rownames(Wsat) <- paste(hru[,'id'])
Wsurf <- Wsat

## set up the data
source('../dynatopmodel/R/xts_util.r')
tmp <- disaggregate_xts(brompton$rain, dt = 15/60)
tmp <- merge(tmp,brompton$pe,all=c(TRUE,FALSE))
tmp <- merge(tmp,brompton$qobs,all=c(TRUE,FALSE))
names(tmp) <- c("rain","pe","q")
obs_data <- tmp

brompton <- list(hru=hru,
                 obs_data=obs_data,
                 Wsurf=Wsurf,
                 Wsat=Wsat,
                 param = param)

save(brompton,file="./data/brompton.rda")

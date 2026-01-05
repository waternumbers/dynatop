## create a single hsu model for testing
rm(list=ls())
library(xts)
devtools::load_all("../../") #library(dynatop)

mdl <- list(
    options = c(channel_solver="histogram"),
    map = list(hillslope="", channel="",channel_id=""),
    point_inflow = data.frame(id = integer(0), name = character(0)),
    diffuse_inflow = data.frame(id = integer(0), name = character(0)),
    gauge = data.frame(name = "outlet", id = as.integer(1)),
    precip_input = data.frame(name="rain",id=1:2,frc=1),
    pet_input = data.frame(name="pet",id=1:2,frc=1),
    channel = data.frame(id=as.integer(1), area=1, length=100, v_ch=1),
    hillslope = data.frame(id = as.integer(2),
                           area = 1, atb_bar = 12, s_bar = 0.001, width = 1,
                           s_sf = 1, s_rz = 0.1, s_uz = 0, s_sz = 0,
                           opt = "exp",
                           r_sfmax = 0, c_sf = 0.1,
                           s_rzmax = 0.1,
                           t_d = 7200,
                           ln_t0 = -2,
                           c_sz = as.numeric(NA), m = 0.04, D = as.numeric(NA), m_2 = as.numeric(NA), omega = as.numeric(NA),
                           s_rz0 = 0.75, r_uz_sz0 = 1e-6,
                           s_raf = 1, t_raf = 10*60*60
                           ),
    flow_direction = data.frame(from = as.integer(2), to = as.integer(1), frc=1)
)

dt <- dynatop$new(mdl,use_states=TRUE)

obs <- xts( matrix(0,100,2), order.by = as.POSIXct("2022-01-01 00:00:00") + (0:99)*900)
names(obs) <- c("rain","pet")

dt$add_data(obs)
dt$sim_hillslope(keep_states=index(obs))
s <- dt$get_states(record=TRUE)
s <- do.call(rbind,s)
s <- xts(as.matrix(s),order.by=index(obs))
head( s[,"s_sf"] )
head(exp( - (1:100)*900/(10*60*60) ))

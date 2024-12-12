rm(list=ls())
source("./R/evap.R")

st <- as.POSIXct("1970-01-02 00:00:00",tz='GMT')
fn <- as.POSIXct("1971-01-01 00:00:00",tz='GMT')
daily_ts <- seq(st,fn,by=24*60*60)
dpet <- evap_est(daily_ts,0,1)

st <- as.POSIXct("1970-01-01 01:00:00",tz='GMT')
fn <- as.POSIXct("1971-01-01 00:00:00",tz='GMT')
hour_ts <- seq(st,fn,by=1*60*60)
hpet <- evap_est(hour_ts,0,1)

## the totals should be the same...
stopifnot(all.equal(sum(hpet), sum(dpet)))
 
## Generating the hourly data with daily input
hpet2 <- evap_est(hour_ts, dailyPET=dpet)

## should be identical
stopifnot(all.equal(hpet, hpet2))


## try updated version
dpet_v <- evap_est_v(daily_ts,0,1)
hpet_v <- evap_est(hour_ts,0,1)
hpet2_v <- evap_est(hour_ts, dailyPET=dpet_v)

stopifnot(all.equal(dpet, dpet_v))
stopifnot(all.equal(hpet, hpet_v))
stopifnot(all.equal(hpet2, hpet2_v))

#' Create sinsuiodal time series of potential evapotranspiration input
#'
#' @description Generate series of potential evapotranspiration
#'
#' @param ts as vector of POSIXct data/times
#' @param eMin Minimum daily PE total (m or mm)
#' @param eMax Maximum daily PE total (m or mm)
#'
#' @details Dynamic TOPMODEL requires a time series of potential
#'   evapotranspiration in order to calculate and remove actual
#'   evapotranspiration from the root zone during a run. Many sophisticated
#'   physical models have been developed for estimating PE and AE, including the
#'   Priestly-Taylor (Priestley and Taylor, 1972) and Penman-Monteith (Montieth,
#'   1965) methods. These, however, require detailed meteorological data such as
#'   radiation input and relative humidities that are, in general, difficult to
#'   obtain. Calder (1983) demonstrated that a simple approximation using a
#'   sinusoidal variation in potential evapotranspiration to be a good
#'   approximation to more complex schemes.
#'
#'   If the insolation is also taken to vary sinusoidally through the daylight
#'   hours then, ignoring diurnal meteorological variations, the potential
#'   evapotranspiration during daylight hours for each year day number can be
#'   calculated (for the catchment's latitude). Integration over the daylight
#'   hours allows the daily maximum to be calculated and thus a sub-daily series
#'   generated.
#' 
#' @return Time series (xts) of potential evapotranspiration totals for the timesteps geven in same units as eMin and eMax
#'
#' @references Beven, K. J. (2012). Rainfall-runoff modelling : the primer. Chichester, UK, Wiley-Blackwell.
#' @references Calder, I. R. (1986). A stochastic model of rainfall interception. Journal of Hydrology, 89(1), 65-71.
#'
#' @examples
#' ## Generating daily PET data for 1970
#' ## the values of eMin and eMax may not by not be realistic
#' st <- as.POSIXct("1970-01-02 00:00:00",tz='GMT')
#' fn <- as.POSIXct("1971-01-01 00:00:00",tz='GMT')
#' daily_ts <- seq(st,fn,by=24*60*60)
#' dpet <- evap_est(daily_ts,0,1)
#'
#' ## create hourly data for the same period
#' st <- as.POSIXct("1970-01-01 01:00:00",tz='GMT')
#' fn <- as.POSIXct("1971-01-01 00:00:00",tz='GMT')
#' hour_ts <- seq(st,fn,by=1*60*60)
#' hpet <- evap_est(hour_ts,0,1)
#'
#' ## the totals should eb the same...
#' stopifnot(sum(hpet)==sum(dpet))
#' @export
evap_est <- function(ts, eMin=0, eMax=0){

    ## checks on input - TODO
    ## fail if more then daily
    dt <- diff(as.numeric(ts))
    if(!all(dt[]==dt[1])){
        stop("Irregularly spaced time series supplied")
    }else{
        dt <- dt[1]
    }

    ## create a series of daily PET values based on eMin and eMax
    ## day 0 is Jan 1, 31 Dec is day 364 or day 365 depending on if leap year
    yday <- 0:365
    fact <- 1+sin(2*pi*yday/365-pi/2)    
    daily_pet <- eMin + 0.5*(eMax-eMin)*fact
    dawn <- (10 - 2.5*fact)*60*60 # in seconds from start of day
    dayLength <- (6 + 4*fact) * 60*60 # in sec from start fo day

    ## make some elements for computation
    sts <- ts-dt # this is the start of the timestep
    dts <- as.POSIXct(3600*24*floor(as.numeric(ts)/(3600*24)),origin="1970-01-01",tz='GMT') # start of the day on which the timestep finishes
    dsts <- as.POSIXct(3600*24*floor(as.numeric(sts)/(3600*24)),origin="1970-01-01",tz='GMT') # start of the day on which the timestep starts

    ##
    pet <- rep(0,length(ts))
    tmp <- dsts
    idx <- tmp<dts
    while( any(idx) ){
        pet[idx] <- pet[idx] + daily_pet[ as.POSIXlt(tmp[idx])$yday + 1]
        tmp <- tmp + 24*60*60
        idx <- tmp<dts
    }
    ## take away the fraction from the start
    sc <- as.numeric(sts) - as.numeric(dsts)
    frc <- (sc - dawn[ as.POSIXlt(sts)$yday +1 ])/ dayLength[ as.POSIXlt(sts)$yday +1]
    frc <- pmin(1,pmax(0,frc))    
    pet <- pet -  daily_pet[ as.POSIXlt(sts)$yday + 1]*0.5*(1-cos(frc*pi))
    ## add fraction from start fo currect day
    sc <- as.numeric(ts) - as.numeric(dts)
    frc <- (sc - dawn[ as.POSIXlt(ts)$yday +1 ])/ dayLength[ as.POSIXlt(ts)$yday +1]
    frc <- pmin(1,pmax(0,frc))
    pet <- pet + daily_pet[ as.POSIXlt(ts)$yday + 1]*0.5*(1-cos(frc*pi))

    return(as.xts(pet,order.by=ts))
}


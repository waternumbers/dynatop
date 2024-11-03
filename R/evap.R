#' Create sinusoidal time series of potential evapotranspiration input
#'
#' @description Generate series of potential evapotranspiration
#'
#' @param ts as vector of POSIXct data/times
#' @param eMin Minimum daily PE total (m or mm)
#' @param eMax Maximum daily PE total (m or mm)
#' @param dailyPET A `xts` series of daily PET totals to use instead of eMin and eMax
#'
#' @details Dynamic TOPMODEL requires a time series of potential
#'   evapotranspiration in order to calculate and remove actual
#'   evapotranspiration from the root zone during a run. Many sophisticated
#'   physical models have been developed for estimating potential and actual evapotranspiration, including the
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
#'   As an alternative to the generation of the sinusoidal daily total is to input
#'   an `xts` object of daily PET values on which to apply the sub-daily sinusoid.
#' 
#' @return Time series (xts) of potential evapotranspiration totals for the time steps given in same units as eMin and eMax, or the daily totals
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
#' ## the totals should be the same...
#' stopifnot(all.equal(sum(hpet), sum(dpet)))
#' 
#' ## Generating the hourly data with daily input
#' hpet2 <- evap_est(hour_ts, dailyPET=dpet)
#'
#' ## should be identical
#' stopifnot(all.equal(hpet, hpet2))
#' 
#' @export
evap_est <- function(ts, eMin=0, eMax=0.03, dailyPET=NULL){

    ## Check min and max
    if(!(eMin < eMax)){
        stop("eMin should be less then eMax")
    }
    
    ## Check timestep
    dt <- diff(as.numeric(ts))
    if(!all(dt[]==dt[1])){
        stop("Irregularly spaced time series supplied")
    }else{
        dt <- dt[1]
    }
    if(dt > 24*60*60){ stop("Time step can be no more then a day") }

    ## for each value in ts we need:
    ## the number of seconds since start of day and date the data is recorded as
    ## recall the daily data is recorded at the date of the midnight at the end!
    ts_sec <- as.numeric(ts) %% 86400 ## seconds since start of day
    ts_date <- as.Date(ts)
    ## represent midnight as the end of a day
    idx <- (ts_sec==0)
    ts_sec[idx] <- 86400
    ts_date[idx] <- ts_date[idx]-1
    ## move on to represent cummulative time indexing
    ts_date <- ts_date + 1

    ## work out range of dates needed
    date_range <- range(ts_date)
    if(ts_sec[1]-dt <0){ date_range[1] <- date_range[1] - 1 }
    
    ## Check if daily pet series is valid
    if(!is.null(dailyPET)){
        stopifnot("dailyPET should be an xts object" = xts::is.xts(dailyPET),
                  "dailyPET index should be at 00:00:00 GMT" = all( (as.numeric(zoo::index(dailyPET)) %% 86400)==0 ),
                  "dailyPET index should have a daily timestep" = all(unique(diff(as.numeric(zoo::index(dailyPET))))==86400),
                  "dailyPET should not start after ts" = zoo::index(dailyPET)[1]<= as.POSIXct(date_range[1],tz="GMT"),
                  "dailyPET should not finish before ts" = as.POSIXct(date_range[2],tz="GMT") >= zoo::index(dailyPET)[nrow(dailyPET)]
                  )
        
        yday <- setNames(as.POSIXlt(zoo::index(dailyPET))$yday-1,
                         format(zoo::index(dailyPET)))
        fact <- 1+sin(2*pi*yday/365-pi/2)
        dailyPET <- as.matrix(dailyPET)
    }else{
    
        tmp <- seq(min(ts_date)-1,max(ts_date),by=1)
        yday <- setNames(as.POSIXlt(tmp)$yday-1,
                         format(tmp))
        
        fact <- 1+sin(2*pi*yday/365-pi/2)    
        dailyPET <- matrix(eMin + 0.5*(eMax-eMin)*fact,dimnames=list(format(tmp),"pet"))
    }

    ## create a series of daily PET values based on eMin and eMax
    ## day 0 is Jan 1, 31 Dec is day 364 or day 365 depending on if leap year
    dawn <- (10 - 2.5*fact)*60*60 # in seconds from start of day
    dayLength <- (6 + 4*fact) * 60*60 # in sec from start fo day

    ## work out cumulative pet in 
    dstr <- format(ts_date)
    frc <- (ts_sec - dawn[ dstr ]) / dayLength[ dstr ]
    frc <- pmin(1,pmax(0,frc))
    pet <- dailyPET[ dstr, ,drop=FALSE]*0.5*(1-cos(frc*pi))

    ## adjust back
    ts_sec <- ts_sec - dt
    idx <- ts_sec < 0
    ts_date[idx] <- ts_date[idx]-1
    dstr <- format(ts_date)
    frc <- (ts_sec - dawn[ dstr ]) / dayLength[ dstr ]
    frc <- pmin(1,pmax(0,frc))
    tmp_pet <- dailyPET[ dstr, ,drop=FALSE]*0.5*(1-cos(frc*pi))
    tmp_pet[idx] <- -tmp_pet[idx]
    pet <- pet - tmp_pet
    return( xts::xts(pet,order.by=ts) )
}

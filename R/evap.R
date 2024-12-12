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
#' @references Calder, I. R., Harding, R. J. & Rosier, P. T. W. (1983) An Objective Assessment of soil-moisture deficit models. Journal of Hydrology, 60, 329-355
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
    ## recall that data recorded at midnight come from the day before
    ts_sec <- as.numeric(ts) %% 86400 ## seconds since start of day
    ts_date <- as.Date(ts)
    ## represent midnight as the end of a day
    idx <- (ts_sec==0)
    ts_sec[idx] <- 86400
    ts_date[idx] <- ts_date[idx]-1

    ## work out which sort of daily pet we are using
    if(!is.null(dailyPET)){
        ## Check if daily pet series is valid
        d_ts <- zoo::index(dailyPET)
        stopifnot("dailyPET should be an xts object" = xts::is.xts(dailyPET),
                  "dailyPET index should have a daily timestep" = all( diff(as.numeric(d_ts))==86400 )
                  )
          
        ## convert to calender days
        d_ts_date <- as.Date(d_ts)
        offset <- as.numeric(d_ts[1]) %% 86400
        if( offset != 0 ){
            warning("Approximating date to allow for offset in daily data")
        }
        
        if( offset <= (86400/2) ){
            ##  recorded before midday on given date - take previous date
            d_ts_date <- d_ts_date - 1
        }
        
        stopifnot(
            "dailyPET should not start after first ts date" = d_ts_date[1] <= ts_date[1],
            "dailyPET should not finish before last ts date" = utils::tail(d_ts_date,1) >= utils::tail(ts_date,1)
        )

        dailyPET <- as.matrix(dailyPET)
        rownames(dailyPET) <- format(d_ts_date)
        yday <- setNames(as.POSIXlt(d_ts_date)$yday,
                         format(d_ts_date)
                         )
        fact <- 1+sin(2*pi*yday/365-pi/2)
    }else{
        ## Check min and max
        if(!(eMin < eMax)){
            stop("eMin should be less then eMax")
        }
        
        ## create a series of daily PET values based on eMin and eMax    
        tmp <- seq(min(ts_date),max(ts_date),by=1)
        yday <- setNames(as.POSIXlt(tmp)$yday,
                         format(tmp))
        
        fact <- 1+sin(2*pi*yday/365-pi/2)    
        dailyPET <- matrix(eMin + 0.5*(eMax-eMin)*fact,dimnames=list(format(tmp),"pet"))
    }
    
    ## day 0 is Jan 1, 31 Dec is day 364 or day 365 depending on if leap year
    dawn <- (10 - 2.5*fact)*60*60 # in seconds from start of day
    dayLength <- (6 + 4*fact) * 60*60 # in sec from start fo day

    ## work out cumulaive pet to start of the day
    dailyCumulative <- apply(dailyPET,2,function(x){ utils::head(cumsum(c(0,x)), -1) })
    rownames(dailyCumulative) <- rownames( dailyPET )
    
    ## work out cumulative pet in 
    dstr <- format(ts_date)
    frc <- (ts_sec - dawn[ dstr ]) / dayLength[ dstr ]
    frc <- pmin(1,pmax(0,frc))
    pet <- dailyCumulative[dstr,,drop=FALSE] + dailyPET[ dstr, ,drop=FALSE]*0.5*(1-cos(frc*pi))

    pet <- apply(pet,2,function(x){ c(x[1],diff(x)) })

    ## ## adjust back
    ## ts_sec <- ts_sec - dt
    ## idx <- ts_sec < 0
    ## ts_date[idx] <- ts_date[idx]-1
    ## dstr <- format(ts_date)
    ## frc <- (ts_sec - dawn[ dstr ]) / dayLength[ dstr ]
    ## frc <- pmin(1,pmax(0,frc))
    ## tmp_pet <- dailyPET[ dstr, ,drop=FALSE]*0.5*(1-cos(frc*pi))
    ## tmp_pet[idx] <- -tmp_pet[idx]
    ## pet <- pet - tmp_pet
    return( xts::xts(pet,order.by=ts) )
}





## ## this is the original implimentation for checking
## evap_est_v <- function(ts, eMin=0, eMax=0.03, dailyPET=NULL){

##     ## Check min and max
##     if(!(eMin < eMax)){
##         stop("eMin should be less then eMax")
##     }
    
##     ## Check timestep
##     dt <- diff(as.numeric(ts))
##     if(!all(dt[]==dt[1])){
##         stop("Irregularly spaced time series supplied")
##     }else{
##         dt <- dt[1]
##     }
##     if(dt > 24*60*60){ stop("Time step can be no more then a day") }

##     ## for each value in ts we need:
##     ## the number of seconds since start of day and date the data is recorded as
##     ## recall the daily data is recorded at the date of the midnight at the end!
##     ts_sec <- as.numeric(ts) %% 86400 ## seconds since start of day
##     ts_date <- as.Date(ts)
##     ## represent midnight as the end of a day
##     idx <- (ts_sec==0)
##     ts_sec[idx] <- 86400
##     ts_date[idx] <- ts_date[idx]-1
##     ## move on to represent cummulative time indexing
##     ts_date <- ts_date + 1

##     ## work out range of dates needed
##     date_range <- range(ts_date)
##     if(ts_sec[1]-dt <0){ date_range[1] <- date_range[1] - 1 }
    
##     ## Check if daily pet series is valid
##     if(!is.null(dailyPET)){
##         stopifnot("dailyPET should be an xts object" = xts::is.xts(dailyPET),
##                   "dailyPET index should be at 00:00:00 GMT" = all( (as.numeric(zoo::index(dailyPET)) %% 86400)==0 ),
##                   "dailyPET index should have a daily timestep" = all(unique(diff(as.numeric(zoo::index(dailyPET))))==86400),
##                   "dailyPET should not start after ts" = zoo::index(dailyPET)[1]<= as.POSIXct(date_range[1],tz="GMT"),
##                   "dailyPET should not finish before ts" = as.POSIXct(date_range[2],tz="GMT") >= zoo::index(dailyPET)[nrow(dailyPET)]
##                   )
        
##         yday <- setNames(as.POSIXlt(zoo::index(dailyPET))$yday-1,
##                          format(zoo::index(dailyPET)))
##         fact <- 1+sin(2*pi*yday/365-pi/2)
##         dailyPET <- as.matrix(dailyPET)
##     }else{
##         ## create a series of daily PET values based on eMin and eMax    
##         tmp <- seq(min(ts_date)-1,max(ts_date),by=1)
##         yday <- setNames(as.POSIXlt(tmp)$yday-1,
##                          format(tmp))
        
##         fact <- 1+sin(2*pi*yday/365-pi/2)    
##         dailyPET <- matrix(eMin + 0.5*(eMax-eMin)*fact,dimnames=list(format(tmp),"pet"))
##     }

##     ## day 0 is Jan 1, 31 Dec is day 364 or day 365 depending on if leap year
##     dawn <- (10 - 2.5*fact)*60*60 # in seconds from start of day
##     dayLength <- (6 + 4*fact) * 60*60 # in sec from start fo day

##     ## work out cumulative pet in 
##     dstr <- format(ts_date)
##     frc <- (ts_sec - dawn[ dstr ]) / dayLength[ dstr ]
##     frc <- pmin(1,pmax(0,frc))
##     pet <- dailyPET[ dstr, ,drop=FALSE]*0.5*(1-cos(frc*pi))

##     ## adjust back
##     ts_sec <- ts_sec - dt
##     idx <- ts_sec < 0
##     ts_date[idx] <- ts_date[idx]-1
##     dstr <- format(ts_date)
##     frc <- (ts_sec - dawn[ dstr ]) / dayLength[ dstr ]
##     frc <- pmin(1,pmax(0,frc))
##     tmp_pet <- dailyPET[ dstr, ,drop=FALSE]*0.5*(1-cos(frc*pi))
##     tmp_pet[idx] <- -tmp_pet[idx]
##     pet <- pet - tmp_pet
##     return( xts::xts(pet,order.by=ts) )
## }

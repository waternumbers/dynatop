#' Functions to resample an xts time series
#'
#' @description Takes an xts time series object and resamples then to a new time step.
#' @param obs A times series (xts) object with a POSIXct index.
#' @param dt New time interval in seconds
#' @param is.rate If TRUE then these are rates i.e m/hr. Otherwise they are absolute values accumulated within the preceeding time interval. Values are scaled before returning so resampling is conservative.
#'
#' @return An xts object with the new timestep
#'
#' @details Time series of observation data are often of different temporal resolutions, however the input to most hydrological models, as is the case with the Dynamic TOPMODEL, requires those data at the same interval. This provides a method to resample a collection of such data to a single interval.
#'
#' Because of the methods used the results:
#' 
#' - are not accurate when the input data does not have a constant timestep. The code issued a warnign and proceeds assuming the data are equally spaced with the modal timestep.
#' - do not guarentee the requested time step but returns a series with the timestep computed from an integer rounding the ratio of the current and requested time step.
#'
#' @details
#'
#' @export
#' @examples
#' # Resample Swindale Rainfall to hourly intervals
#' require(dynatop)
#' data("Swindale")
#' obs <- Swindale$obs
#' cobs <- resample_xts(obs, dt=60*60) # hourly data
#' dobs <- resample_xts(cobs,dt=15*60) # back to 15 minute data
#' cdobs <- resample_xts(dobs,dt=60*60) # back to hourly data - checks time stamp conversion
#' obs <- obs[zoo::index(obs)<=max(zoo::index(cobs)),]
#' 
#' # check totals
#' stopifnot( all.equal(sum(obs),sum(cobs)) )
#' stopifnot( all.equal(sum(obs),sum(dobs)) )
#' stopifnot( all.equal(cobs,cdobs) )
resample_xts <- function(obs, dt, is.rate=FALSE){

    ## if the set is NULL then return
    if(is.null(obs)){return(obs)}
    if(!zoo::is.zoo(obs)){stop("Time series required")}
    if(is.null(dt)){stop("Supply new time interval")}
    tms <- as.double(index(obs))
    dt_series <- diff(tms)
    if(!all(dt_series[]==dt_series[1])){
        warning(paste("Irregularly spaced time series supplied to resample_xts -",
                      "proceding with modal timestep, results are questionable"))
        dt_series <- min(Modes(dt_series))
    }else{
        dt_series <- dt_series[1]
    }


    ## if factor greater then 1 take to smaller timestep
    if(dt_series >= dt){
        ## disaggregation factor (nearest integer)
        tryCatch({ fact <- ceiling(dt_series/dt) },
                 error = {fact <- 1})
        ## disaggregating observations from larger to smaller time interval by given factor
        ## recalc interval
        dt <- dt_series/fact
        ## new values duplicate series and rescale if not a rate in terms of a fixed
        ## period e.g m/hr
        vals <- apply(as.matrix(obs), MARGIN=2, FUN=rep, each=fact)   # won't wotk if fact is not an integer
        vals <- matrix(vals, ncol=ncol(vals))
        tms <- seq(index(obs)[1]-((fact-1)*dt), along.with=vals, by=dt)
        ## if the value is a rate then it should be applied to all of the values in
        ## the interval "as is". Otherwise each values needs to be divided across the
        ## smaller time steps so that the total across the original time intervals is the same
        obs_agg <- xts::xts(vals, order.by=tms[1:nrow(vals)])
        if(!is.rate){obs_agg <- obs_agg/fact}
        names(obs_agg) <- names(obs)
    }else{
        ## then we need to aggregate (e.g. quarterly to hourly)
        tryCatch({ fact <- floor(dt/dt_series) },
                 error = {fact <- 1})
        tms <- index(obs)[seq(fact,nrow(obs),by=fact)]
        # tms <- seq(start(obs), end(obs), by=dt)
        idx <- rep(tms, each=fact)
        nn <- min(length(idx),nrow(obs))
        idx <- idx[1:nn]
        obs <- obs[1:nn,]
        ##idx <- idx[1:nrow(obs)]
        if(is.rate){fun=mean}else{fun=sum}
        obs_agg <- as.xts(aggregate(zoo::zoo(obs), by = idx, FUN=fun))
        names(obs_agg) <- names(obs)
        
    }

    return(obs_agg)
}


Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  mtab <- max(tab)
  if( max(tab)==1 ){
      ## issue warning
      warning("dynatop::Modes - All values occur once, Modes are poorly defined")
  }
  ux[tab == max(tab)]
}
